# packJPG test suite

A test suite for packJPG, added where there was none. It is organised in the
three layers requested, in priority order.

## Running

From the `source/` directory:

```
make test        # build + run the C++ suite AND the CLI golden test
make test-cli    # just the CLI golden round-trip against the real binary
make test-asan   # deterministic AddressSanitizer reproduction of the model_b bug
make fixtures    # regenerate the JPEG corpus (needs python3+Pillow, ImageMagick)
```

`make test` exits 0 when the suite is green. The JPEG fixtures are committed, so
running the tests needs only a C++14 compiler; `make fixtures` is only required
when you want to change the corpus.

## Layout

```
tests/
  test_framework.h    tiny single-header framework (TEST_CASE / CHECK / REQUIRE)
  test_main.cpp       defines main() (runs all auto-registered cases)
  test_roundtrip.cpp  Layer 1: in-process pjglib round-trip + error handling
  cli_golden.sh       Layer 1: end-to-end round-trip against the real binary
  test_bitops.cpp     Layer 2: BitReader/BitWriter/Memory{Reader,Writer} units
  test_aricoder.cpp   Layer 2: arithmetic coder encode/decode inverse + guards
  gen_fixtures.py     regenerates fixtures/ deterministically
  fixtures/valid/     libjpeg-encoded JPEGs (grayscale/color/progressive/cmyk/...)
  fixtures/invalid/   malformed inputs (empty, truncated, non-JPEG, corrupt PJG)
```

### Why a home-grown framework instead of googletest/doctest/Catch2?

The brief asked for a lightweight, single-header framework with no new build
system. `test_framework.h` is ~150 lines and needs nothing vendored or fetched
(the build sandbox has no package network to pull a dependency). It deliberately
mirrors the `TEST_CASE` / `CHECK` / `REQUIRE` spelling of doctest, so switching to
doctest later is a near drop-in change. googletest was considered but it requires
building a library and a CMake/Bazel step, which conflicts with the
single-header / no-new-build-system constraint.

## Layer 1 — round-trip / golden tests (highest priority)

The core packJPG guarantee is a lossless round-trip: `JPG -> PJG -> JPG` is
byte-identical (`Readme.txt:85-91`). This is covered two ways:

* **In-process** (`test_roundtrip.cpp`) through the library API
  `pjglib_init_streams` + `pjglib_convert_stream2mem` (`packjpg.cpp:982,903`) in
  memory→memory mode, so there are no temp files. It also checks compression and
  decompression determinism and that every malformed input is rejected with a
  message and without crashing.
* **End-to-end** (`cli_golden.sh`) by driving the real `packJPG` binary over the
  whole corpus and byte-comparing the reconstruction, plus feeding it the
  malformed inputs to confirm it does not crash on them.

## Layer 2 — unit tests for self-contained pieces

`bitops` (`BitReader`/`BitWriter`/`Memory{Reader,Writer}`) and `aricoder`
(`ArithmeticEncoder`/`ArithmeticDecoder` + `model_b`/`model_s`) carry no
file-scope global state, so they are tested in isolation. The properties tested
are inverse/round-trip ones (what a writer writes, a reader reads back; what a
model encodes, an identically-initialised model decodes) plus boundary and
buffer-edge behaviour.

The rest of the engine in `packjpg.cpp` is a wall of `INTERN static` globals
(`packjpg.cpp:295,652` etc.), and its core routines (`read_jpeg`, `decode_jpeg`,
…) are `static`, so they cannot be linked or exercised in isolation without
refactoring. Per the brief this was left alone and covered at the round-trip
boundary instead.

## Findings (bugs discovered while building this)

1. **Heap-buffer-overflow in `model_b::update_model` (`aricoder.cpp:629`) —
   FIXED.** Symptom: during PJG compression of many images, `context->counts[
   symbol]++` indexed far past the model's 2-entry count table (ASan showed a
   read at `symbol == 255`). Usually benign (lands on harmless adjacent heap and
   output is still correct) but it **intermittently segfaulted** depending on
   heap layout — larger inputs / longer working paths made a crash more likely.

   Stack:
   ```
   model_b::update_model(int)  aricoder.cpp:629
   encode_ari                  aricoder.h:345
   pjg_encode_bit              packjpg.cpp:5321
   pack_pjg                    packjpg.cpp:3318   <- pjg_encode_bit(encoder, padbit)
   ```

   **Root cause: `char` signedness.** `padbit` was declared `char`
   (`packjpg.cpp:564`) and uses `-1` as an "unset" sentinel. Where `char` is
   *unsigned* (e.g. ARM/aarch64) the sentinel is stored as `255`, so the guard
   `if (padbit == -1) padbit = 1;` (`pack_pjg`, packjpg.cpp:3307) never matches
   (`255 == -1` is false after integer promotion). An unset `padbit` (255) was
   then passed to `pjg_encode_bit(..., unsigned char bit)` and encoded as a 1-bit
   symbol through a binary model — overflowing the 2-entry count table. It only
   bites JPEGs whose final Huffman scan ends byte-aligned (so `unpad` never
   assigns a real 0/1 pad bit). x86 (`signed char`) silently did the right thing,
   which is why upstream never saw it.

   **Fix:** declare the sentinel `signed char padbit` (`packjpg.cpp:564`), which
   restores the intended cross-platform behaviour. Bonus: the same 255-vs-(-1)
   confusion also made the "inconsistent use of padbits" consistency check
   (packjpg.cpp:2727) misfire on ARM, spuriously rejecting valid files under
   default options; with the fix, all fixtures now compress and round-trip under
   default options too (the in-process suite went from 3/9 to 9/9 accepted).
   `make test-asan` is the regression guard and now PASSES.

2. **`new[]` / `free()` alloc-dealloc mismatch in `BitWriter::get_c_bytes()`
   (`bitops.cpp`) — FIXED.** The Huffman buffer from `get_c_bytes()` was
   allocated with `new unsigned char[]`, but the same `huffdata` pointer is
   released with `free()` (`reset_buffers`, packjpg.cpp:2018) — a mismatched
   allocator (undefined behaviour), surfaced by ASan on the decompress path.
   Fixed by allocating with `malloc`, matching the `free()` convention and the
   sibling `get_c_data()` (which already used `malloc`).

3. **Negative shift UB in `decode_jpeg` (`packjpg.cpp:3962`, `:3980`) — OPEN,
   benign.** UBSan reports "shift exponent -1 is negative" via the `DEVLI` macro
   when a Huffman size code is 0. It does not corrupt output (round-trips stay
   bit-identical) and is left as a documented, separate finding rather than
   risking a behavioural change in the coder. `make test-asan` intentionally runs
   AddressSanitizer only, so this UB does not mask memory-safety regressions.

4. **`packJPG` returns exit code 0 even on error.** The CLI prints an error
   summary but still exits 0 (observed on every malformed input). The reliable
   success signal is therefore "did a correct output file get produced", not
   `$?`. The golden test asserts on output existence + byte-compare accordingly,
   and treats a *crash* (rc ≥ 128) as a hard failure.

5. **`packJPG` is picky about JPEG Huffman pad bits.** Files whose final Huffman
   byte is padded with 0-bits (e.g. Pillow's encoder) are rejected outright with
   "inconsistent use of padbits" (`packjpg.cpp:2726-2732`). libjpeg (ImageMagick
   `convert`) pads per spec and is accepted, so the corpus is libjpeg-encoded.
   `-p` (proceed-on-warnings) makes packJPG use padbit=1 — exactly what libjpeg
   wrote — so round-trips stay bit-identical; the CLI test uses `-p` throughout.

6. **Out-of-bounds heap write in `pjg_unoptimize_header` (`packjpg.cpp`) —
   FIXED.** On the decompress path the reconstructed JFIF header is walked with
   segment lengths and Huffman/quant "skip" counts taken straight from the
   (untrusted) decompressed stream. A malformed DHT/DQT segment could drive
   `hpos` past the end of the `hdrdata` allocation, so the in-place rewrites
   (`hdrdata[hpos+spos] += …`, std-table reinsertion) scribbled over adjacent
   heap and corrupted the allocator — the process then aborted inside `free()`
   in `reset_buffers` (`EXC_BREAKPOINT`/SIGTRAP on macOS, SIGSEGV under a
   different heap layout). The std-table index `i = hdrdata[hpos+1]` was also
   used to read `std_huff_lengths[i]`/`std_huff_tables[i]` (size-4 arrays)
   without bounding `i`. **Fix:** bound every `hdrdata` access against `hdrs`
   (and `i` against `[0,4)`) in both `pjg_unoptimize_header` and its encode-side
   mirror `pjg_optimize_header`; a violation is reported via `pjg_header_error()`
   (sets `errorlevel`, so the pipeline actually halts instead of proceeding into
   `recode_jpeg` with half-initialized state). Regression fixture:
   `fixtures/invalid/pjg_header_oob.pjg`.

## Notes / limitations

* The shipped `docs/sample_images.zip` contains only PNG coefficient dumps, not
  JPEGs, so the corpus is generated from scratch (`gen_fixtures.py`).
* Compressing the same JPEG again *after* a decompress in the same process is not
  a documented invariant and empirically exposes residual global state, so the
  suite tests only the documented single round-trip, not re-compression
  stability.

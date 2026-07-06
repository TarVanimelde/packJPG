// test_roundtrip.cpp - Layer 1: round-trip / golden tests (highest priority).
//
// The core packJPG guarantee is a LOSSLESS round-trip: compressing a JPEG to PJG
// and decompressing back yields a byte-identical JPEG (Readme.txt:85-91). These
// tests exercise that invariant *in-process* through the library API, so there is
// no dependency on temp files or the CLI:
//
//     pjglib_init_streams(in, in_type, in_size, out, out_type)   (packjpg.cpp:982)
//     pjglib_convert_stream2mem(&out, &out_size, msg)            (packjpg.cpp:903)
//
// with in_type/out_type == 1 meaning "memory" (packjpglib.h:19-40). The filetype
// (compress vs decompress) is auto-detected from the magic bytes: 0xFFD8 => JPEG
// (compress), 'J''S' => PJG (decompress) (packjpg.cpp:1083-1100, :711).
//
// Note on the padbit quirk: under DEFAULT options packJPG rejects some libjpeg
// JPEGs with "inconsistent use of padbits" (packjpg.cpp:2726-2732); whether a
// given file trips it is content-dependent. Rather than depend on which files
// pass, this test *probes* acceptance: refused-with-padbits files are recorded
// and skipped, accepted files MUST round-trip bit-identically, and we REQUIRE at
// least one file to actually round-trip so the suite can never pass vacuously.
// The CLI golden test (tests/cli_golden.sh) covers the whole corpus using `-p`.

#include "test_framework.h"

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>

// --- packJPG library API (compiled with -DBUILD_LIB; EXPORT == extern) -------
extern bool pjglib_convert_stream2mem(unsigned char** out_file,
                                      unsigned int* out_size, char* msg);
extern void pjglib_init_streams(void* in_src, int in_type, int in_size,
                                void* out_dest, int out_type);
extern const char* pjglib_version_info(void);

#ifndef FIXTURE_DIR
#define FIXTURE_DIR "tests/fixtures"
#endif

namespace {

const char* kValid[] = {
    "grayscale.jpg", "grayscale_large.jpg", "baseline_rgb.jpg",
    "baseline_444.jpg", "baseline_q20.jpg", "progressive.jpg",
    "cmyk.jpg", "odd_dimensions.jpg", "large_256.jpg",
};

const char* kInvalid[] = {
    "empty.bin", "not_a_jpeg.txt", "bad_after_soi.jpg",
    "truncated.jpg", "soi_eoi_only.jpg", "fake_pjg.pjg",
    // Fuzzed PJG that drove model_s::current_order out of bounds -> segfault in
    // totalize_table before the fix (issue #41).
    "issue41_context_order.pjg",
    // Fuzzed PJG whose decompressed JFIF header has a DQT/DHT segment whose
    // length drives pjg_unoptimize_header past the hdrdata allocation: an
    // out-of-bounds heap write that corrupted the allocator and crashed on the
    // next free() (SIGTRAP/SIGSEGV) before the bounds checks were added.
    "pjg_header_oob.pjg",
    // Valid PJG truncated to 40 bytes: the arithmetic stream runs out before the
    // generic header decoder sees its 256 terminator. read_bit used to feed zero
    // bits forever, so pjg_decode_generic looped (and grew) indefinitely -> a DoS
    // hang. Must now be rejected promptly; if this test hangs, that fix regressed.
    "pjg_truncated_hang.pjg",
    // JPEG whose SOF component quantization-table selector (Tq) is 4 (valid range
    // 0..3). It indexed qtables[4][64] out of bounds and handed back a wild
    // pointer dereferenced in jpg_setup_imginfo -> global-buffer-overflow under
    // ASan (upstream issues #23/#32, follow-on SEGVs #27/#35).
    "jpg_qtable_index_oob.jpg",
    // Fuzzed JPEG that passes read_jpeg but fails jpg_parse_jfif inside
    // decode_jpeg: the early "return false" there used to leak the BitReader
    // (`huffr`, 32 bytes) -> LeakSanitizer report (upstream issue #34). Kept as a
    // rejection test here; `make test-asan` (detect_leaks=1) is the leak guard.
    "jpg_decode_huffr_leak.jpg",
    // Fuzzed JPEG with a DHT whose 16 code-length counts sum past the end of the
    // marker segment: jpg_build_huffcodes walked the code-value bytes off the end
    // of the hdrdata buffer -> heap-buffer-overflow READ under ASan (upstream
    // issues #26/#33). `make test-asan` is the guard.
    "jpg_dht_overflow.jpg",
    // Fuzzed JPEG whose corrupt DC Huffman table decodes an out-of-range size
    // code s (> 15), which drove oversized shifts `1 << s` / read(s) with s > 31
    // -> undefined behavior under UBSan. `make test-ubsan` is the guard.
    "jpg_dc_size_ub.jpg",
};

using Bytes = std::vector<std::uint8_t>;

Bytes read_file(const std::string& path) {
    FILE* f = std::fopen(path.c_str(), "rb");
    if (!f) return Bytes{};
    std::fseek(f, 0, SEEK_END);
    long n = std::ftell(f);
    std::fseek(f, 0, SEEK_SET);
    Bytes v(n > 0 ? static_cast<size_t>(n) : 0);
    if (n > 0 && std::fread(v.data(), 1, v.size(), f) != v.size()) v.clear();
    std::fclose(f);
    return v;
}

std::string valid_path(const char* name) {
    return std::string(FIXTURE_DIR) + "/valid/" + name;
}
std::string invalid_path(const char* name) {
    return std::string(FIXTURE_DIR) + "/invalid/" + name;
}

struct ConvertResult {
    bool ok = false;
    Bytes out;
    std::string msg;
};

// Runs one pjglib memory->memory conversion. `in` is copied so the library may
// read/consume it freely. Frees the library-allocated output buffer.
ConvertResult convert_mem(const Bytes& in) {
    ConvertResult r;
    char msg[128] = {0};  // MSG_SIZE == 128 (packjpg.cpp:327)
    Bytes in_copy = in;   // library reads from this buffer in place
    pjglib_init_streams(in_copy.empty() ? (void*)"" : (void*)in_copy.data(),
                        1, static_cast<int>(in_copy.size()), nullptr, 1);
    unsigned char* out = nullptr;
    unsigned int out_size = 0;
    r.ok = pjglib_convert_stream2mem(&out, &out_size, msg);
    r.msg = msg;
    if (r.ok && out != nullptr && out_size > 0) {
        r.out.assign(out, out + out_size);
    }
    if (out) std::free(out);
    return r;
}

bool is_padbit_refusal(const ConvertResult& r) {
    return !r.ok && r.msg.find("padbit") != std::string::npos;
}

}  // namespace

TEST_CASE("library API is linked and reports a version") {
    const char* v = pjglib_version_info();
    REQUIRE(v != nullptr);
    CHECK(std::strstr(v, "packJPG") != nullptr);
}

TEST_CASE("round-trip: accepted JPEGs decompress bit-identically") {
    int round_tripped = 0, padbit_skipped = 0;
    for (const char* name : kValid) {
        Bytes jpg = read_file(valid_path(name));
        REQUIRE(!jpg.empty());                 // fixture must exist
        CHECK(jpg[0] == 0xFF && jpg[1] == 0xD8);  // is a JPEG (SOI)

        ConvertResult comp = convert_mem(jpg);  // JPG -> PJG
        if (is_padbit_refusal(comp)) {
            ++padbit_skipped;
            std::printf("      (skip %s: default-options padbit refusal)\n", name);
            continue;
        }
        REQUIRE(comp.ok);
        REQUIRE(!comp.out.empty());
        CHECK(comp.out[0] == 'J' && comp.out[1] == 'S');  // PJG magic
        CHECK(comp.out.size() < jpg.size());              // it compresses

        ConvertResult decomp = convert_mem(comp.out);     // PJG -> JPG
        REQUIRE(decomp.ok);
        // The whole point: reconstruction is byte-for-byte identical.
        CHECK_EQ(decomp.out.size(), jpg.size());
        CHECK(decomp.out == jpg);
        if (decomp.out == jpg) ++round_tripped;
        else std::printf("      MISMATCH on %s\n", name);
    }
    std::printf("      round-tripped %d file(s), skipped %d padbit refusal(s)\n",
                round_tripped, padbit_skipped);
    // Guard against a vacuous pass (e.g. every file skipped).
    REQUIRE(round_tripped >= 1);
}

TEST_CASE("compression is deterministic (same input -> same PJG)") {
    // Use a file known to be accepted under default options.
    Bytes jpg = read_file(valid_path("grayscale_large.jpg"));
    REQUIRE(!jpg.empty());
    ConvertResult a = convert_mem(jpg);
    ConvertResult b = convert_mem(jpg);
    REQUIRE(a.ok);
    REQUIRE(b.ok);
    CHECK(a.out == b.out);
}

TEST_CASE("decompression is deterministic (same PJG -> same JPG)") {
    // Decompressing the same PJG twice must reproduce the same JPEG.
    // NOTE: we deliberately do NOT test the compress->decompress->compress chain
    // here. Re-compressing after a decompress in the same process is not a
    // documented invariant and empirically exposes residual global state in the
    // pjglib path (the second compression can differ). The stable, documented
    // guarantee is the single round-trip above.
    Bytes jpg = read_file(valid_path("grayscale_large.jpg"));
    REQUIRE(!jpg.empty());
    ConvertResult pjg = convert_mem(jpg);
    REQUIRE(pjg.ok);
    ConvertResult a = convert_mem(pjg.out);
    ConvertResult b = convert_mem(pjg.out);
    REQUIRE(a.ok);
    REQUIRE(b.ok);
    CHECK(a.out == jpg);
    CHECK(a.out == b.out);
}

TEST_CASE("error handling: malformed inputs are rejected without crashing") {
    for (const char* name : kInvalid) {
        Bytes bad = read_file(invalid_path(name));
        // empty.bin is legitimately zero bytes; others are non-empty.
        ConvertResult r = convert_mem(bad);
        CHECK_FALSE(r.ok);            // must not claim success
        CHECK(!r.msg.empty());        // must report a reason
        std::printf("      %-20s -> rejected: %s\n", name, r.msg.c_str());
    }
}

TEST_CASE("error handling: unknown filetype is reported") {
    Bytes txt = read_file(invalid_path("not_a_jpeg.txt"));
    REQUIRE(!txt.empty());
    ConvertResult r = convert_mem(txt);
    REQUIRE(!r.ok);
    CHECK(r.msg.find("unknown") != std::string::npos ||
          r.msg.find("filetype") != std::string::npos);
}

#!/usr/bin/env bash
# cli_golden.sh - Layer 1 end-to-end golden test against the real packJPG binary.
#
# Complements the in-process pjglib test (tests/run_tests) by driving the actual
# command-line tool over the WHOLE valid/ corpus, using `-p` (proceed-on-warnings)
# so that the content-dependent "inconsistent use of padbits" warning does not
# abort compression. Because libjpeg pads with 1-bits and `-p` makes packJPG use
# padbit=1, reconstruction stays byte-for-byte identical (Readme.txt:73-91).
#
# It also checks that malformed inputs are handled without crashing.
#
# NOTE (documented quirk / finding): packJPG returns exit code 0 even when it
# reports errors, so the reliable success signal is "did a correct output file
# get produced", not $?. We assert on file existence + byte-compare.
#
# A crash (rc >= 128) is always a hard FAIL. This used to XFAIL a known heap-
# buffer-overflow in model_b::update_model (aricoder.cpp:629) that segfaulted on
# byte-aligned JPEGs; that bug is now FIXED (padbit char-signedness, packjpg.cpp:
# 564 -- see tests/README.md "Findings"). The crash check below is retained as a
# regression guard so the bug can never silently return.
#
# Usage: bash cli_golden.sh /path/to/packJPG
set -u

BIN="${1:-./packJPG}"
HERE="$(cd "$(dirname "$0")" && pwd)"
VALID="$HERE/fixtures/valid"
INVALID="$HERE/fixtures/invalid"

if [ ! -x "$BIN" ]; then echo "FATAL: packJPG binary not found/executable at '$BIN'"; exit 2; fi

pass=0; fail=0
red()   { printf '  \033[31mFAIL\033[0m %s\n' "$1"; fail=$((fail+1)); }
green() { printf '  \033[32mok  \033[0m %s\n' "$1"; pass=$((pass+1)); }

echo "== CLI round-trip (valid corpus, -p) =="
for f in "$VALID"/*.jpg; do
    name="$(basename "$f")"
    work="$(mktemp -d)"
    cp "$f" "$work/img.jpg"
    cp "$f" "$work/orig.jpg"                      # keep a pristine copy
    "$BIN" -p -np "$work/img.jpg" >/dev/null 2>&1 # -> img.pjg
    crc=$?
    if [ "$crc" -ge 128 ]; then
        red "$name (compress CRASHED rc=$crc -- regression of the model_b/padbit bug?)"
        rm -rf "$work"; continue
    fi
    if [ ! -f "$work/img.pjg" ]; then red "$name (no .pjg produced, rc=$crc)"; rm -rf "$work"; continue; fi
    rm -f "$work/img.jpg"                          # clear the way for reconstruction
    "$BIN" -p -np "$work/img.pjg" >/dev/null 2>&1  # -> img.jpg
    drc=$?
    if [ "$drc" -ge 128 ]; then
        red "$name (decompress CRASHED rc=$drc)"
    elif [ -f "$work/img.jpg" ] && cmp -s "$work/orig.jpg" "$work/img.jpg"; then
        green "$name  ($(wc -c <"$f") -> $(wc -c <"$work/img.pjg") bytes)"
    else
        red "$name (reconstruction not bit-identical)"
    fi
    rm -rf "$work"
done

echo "== CLI bad-input handling (must not crash) =="
for f in "$INVALID"/*; do
    name="$(basename "$f")"
    work="$(mktemp -d)"
    cp "$f" "$work/in.bin"
    timeout 20 "$BIN" -p -np "$work/in.bin" >/dev/null 2>&1
    rc=$?
    if [ "$rc" -ge 128 ]; then
        red "$name (crashed / killed on malformed input, rc=$rc)"
    else
        green "$name (handled cleanly, rc=$rc)"
    fi
    rm -rf "$work"
done

echo "== CLI output directory (-out) =="
# a representative fixture is enough; the redirect logic is content-independent
sample="$(ls "$VALID"/*.jpg | head -1)"
if [ -n "$sample" ]; then
    work="$(mktemp -d)"
    mkdir -p "$work/src" "$work/out" "$work/out2"
    cp "$sample" "$work/src/img.jpg"; cp "$sample" "$work/orig.jpg"
    # (1) compress into out/ -> output must land there, NOT next to the input
    "$BIN" -p -np -o -out "$work/out" "$work/src/img.jpg" >/dev/null 2>&1
    if [ -f "$work/out/img.pjg" ] && [ ! -f "$work/src/img.pjg" ]; then
        green "-out redirects compressed output into the target directory"
    else
        red "-out did not place img.pjg in the target directory"
    fi
    # (2) decompress with -out into out2/ and byte-compare to the original
    "$BIN" -p -np -o -out "$work/out2" "$work/out/img.pjg" >/dev/null 2>&1
    if [ -f "$work/out2/img.jpg" ] && cmp -s "$work/orig.jpg" "$work/out2/img.jpg"; then
        green "-out round-trip is bit-identical"
    else
        red "-out round-trip not bit-identical"
    fi
    # (3) a missing output directory must be reported clearly and produce no file
    out="$("$BIN" -p -np -o -out "$work/nope" "$work/src/img.jpg" 2>&1)"
    if printf '%s' "$out" | grep -qi "output directory does not exist" && [ ! -f "$work/nope/img.pjg" ]; then
        green "-out with a missing directory errors clearly and writes nothing"
    else
        red "-out with a missing directory did not report the expected error"
    fi
    rm -rf "$work"
fi

echo
echo "CLI golden: $pass passed, $fail failed."
[ "$fail" -eq 0 ]

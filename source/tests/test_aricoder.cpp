// test_aricoder.cpp - Layer 2 unit tests for the arithmetic coder (aricoder.h).
//
// The coder pairs an ArithmeticEncoder(Writer&) with an ArithmeticDecoder(Reader&)
// (aricoder.h:176-221) driven by a statistical model (model_b binary, model_s
// multi-symbol). The property under test is the fundamental one: whatever a model
// encodes, an identically-initialised model must decode back exactly. We route
// bytes through Memory{Writer,Reader} so the whole thing stays in memory.
//
// The encode/decode helpers encode_ari()/decode_ari() are the same generic
// drivers the engine uses (aricoder.h:304-362); the model constructors mirror the
// engine's INIT_MODEL_B/INIT_MODEL_S macros (packjpg.cpp:297-298).

#include "test_framework.h"
#include "../aricoder.h"
#include "../bitops.h"

#include <cstdint>
#include <vector>

namespace {

// Encode a binary sequence with a single-context binary model (max_order 0, so
// no shift_context needed -- like `INIT_MODEL_B(1, 0)` in packjpg.cpp:4862).
std::vector<std::uint8_t> encode_bits(const std::vector<int>& bits) {
    MemoryWriter mw;
    ArithmeticEncoder enc(mw);
    model_b m(1, 0, 255);
    for (int b : bits) encode_ari(&enc, &m, b);
    enc.finalize();          // flush the coder into the writer
    return mw.get_data();
}

std::vector<int> decode_bits(const std::vector<std::uint8_t>& bytes, std::size_t n) {
    MemoryReader mr(bytes);
    ArithmeticDecoder dec(mr);
    model_b m(1, 0, 255);
    std::vector<int> out;
    out.reserve(n);
    for (std::size_t i = 0; i < n; i++) out.push_back(decode_ari(&dec, &m));
    return out;
}

}  // namespace

TEST_CASE("binary model: encode/decode is a round-trip inverse") {
    std::vector<int> bits;
    // A deterministic, skewed pattern so the adaptive model actually adapts.
    for (int i = 0; i < 500; i++) bits.push_back(((i * 7 + 3) % 11) < 4 ? 1 : 0);
    std::vector<std::uint8_t> enc = encode_bits(bits);
    REQUIRE(!enc.empty());
    std::vector<int> dec = decode_bits(enc, bits.size());
    CHECK(dec == bits);
}

TEST_CASE("binary model: all-zeros and all-ones compress and invert") {
    std::vector<int> zeros(256, 0), ones(256, 1);
    CHECK(decode_bits(encode_bits(zeros), zeros.size()) == zeros);
    CHECK(decode_bits(encode_bits(ones), ones.size()) == ones);
    // A highly skewed stream should code to far fewer than 256 bytes.
    CHECK(encode_bits(zeros).size() < 64u);
}

TEST_CASE("binary model: single-symbol streams invert") {
    CHECK(decode_bits(encode_bits({0}), 1) == std::vector<int>{0});
    CHECK(decode_bits(encode_bits({1}), 1) == std::vector<int>{1});
}

TEST_CASE("multi-symbol model: encode/decode is a round-trip inverse") {
    const int kAlphabet = 16;
    std::vector<int> syms;
    for (int i = 0; i < 400; i++) syms.push_back((i * 5 + 1) % kAlphabet);

    // Single context (max_context 1, max_order 0) -> no shift_context needed.
    std::vector<std::uint8_t> bytes;
    {
        MemoryWriter mw;
        ArithmeticEncoder enc(mw);
        model_s m(kAlphabet, 1, 0, 255);
        for (int c : syms) encode_ari(&enc, &m, c);
        enc.finalize();
        bytes = mw.get_data();
    }
    REQUIRE(!bytes.empty());

    std::vector<int> out;
    {
        MemoryReader mr(bytes);
        ArithmeticDecoder dec(mr);
        model_s m(kAlphabet, 1, 0, 255);
        for (std::size_t i = 0; i < syms.size(); i++) out.push_back(decode_ari(&dec, &m));
    }
    CHECK(out == syms);
}

TEST_CASE("multi-symbol model: context shifting round-trips (regression guard)") {
    // Exercises shift_context on both sides -- the shape of real engine usage
    // (e.g. packjpg.cpp:4758-4760). Encoder and decoder must shift identically.
    const int kAlphabet = 8;
    const int kContexts = 4;
    std::vector<int> syms;
    for (int i = 0; i < 300; i++) syms.push_back((i * 3) % kAlphabet);
    auto ctx_of = [](std::size_t i) { return static_cast<int>(i % 4); };

    std::vector<std::uint8_t> bytes;
    {
        MemoryWriter mw;
        ArithmeticEncoder enc(mw);
        model_s m(kAlphabet, kContexts, 1, 255);
        for (std::size_t i = 0; i < syms.size(); i++) {
            m.shift_context(ctx_of(i));
            encode_ari(&enc, &m, syms[i]);
        }
        enc.finalize();
        bytes = mw.get_data();
    }

    std::vector<int> out;
    {
        MemoryReader mr(bytes);
        ArithmeticDecoder dec(mr);
        model_s m(kAlphabet, kContexts, 1, 255);
        for (std::size_t i = 0; i < syms.size(); i++) {
            m.shift_context(ctx_of(i));
            out.push_back(decode_ari(&dec, &m));
        }
    }
    CHECK(out == syms);
}

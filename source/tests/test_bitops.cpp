// test_bitops.cpp - Layer 2 unit tests for the bit/byte I/O layer (bitops.h).
//
// These classes are the most self-contained pieces in packJPG: BitReader/
// BitWriter and the Memory{Reader,Writer} carry no file-scope global state, so
// they can be exercised in isolation. Both bit classes are MSB-first
// (bitops.cpp: BitReader::read / BitWriter::write_u16), which the round-trip
// tests below rely on.

#include "test_framework.h"
#include "../bitops.h"

#include <cstdint>
#include <vector>

namespace {
struct Field { std::uint16_t val; std::size_t bits; };
}

TEST_CASE("BitWriter/BitReader are inverses for multi-width fields") {
    const Field fields[] = {
        {0b101, 3}, {0xA, 4}, {0x1FF, 9}, {0, 1}, {1, 1}, {0x7F, 7}, {0xFFFF, 16},
    };
    BitWriter w(0);
    for (const auto& f : fields) w.write_u16(f.val, f.bits);
    std::vector<std::uint8_t> bytes = w.get_bytes();  // pads the final byte

    BitReader r(bytes.data(), static_cast<int>(bytes.size()));
    for (const auto& f : fields) {
        unsigned int got = r.read(static_cast<int>(f.bits));
        CHECK_EQ((int)got, (int)f.val);
    }
    CHECK_FALSE(r.eof());  // exactly the written bits consumed (+ pad remains)
}

TEST_CASE("write_bit/read_bit preserve MSB-first bit order") {
    const std::uint8_t pattern[] = {1,0,1,1,0,0,0,1, 1,1,0,1,0};
    BitWriter w(0);
    for (std::uint8_t b : pattern) w.write_bit(b);
    std::vector<std::uint8_t> bytes = w.get_bytes();

    BitReader r(bytes.data(), static_cast<int>(bytes.size()));
    for (std::uint8_t b : pattern) CHECK_EQ((int)r.read_bit(), (int)b);
}

TEST_CASE("num_bytes_written counts only fully written bytes") {
    BitWriter w(0);
    CHECK_EQ((int)w.num_bytes_written(), 0);
    w.write_u16(0xFF, 8);
    CHECK_EQ((int)w.num_bytes_written(), 1);
    w.write_bit(1);                       // partial second byte
    CHECK_EQ((int)w.num_bytes_written(), 1);
    std::vector<std::uint8_t> bytes = w.get_bytes();  // pad completes byte 2
    CHECK_EQ((int)bytes.size(), 2);
}

TEST_CASE("pad fills the final byte with the pad bit") {
    // Write a single 1 bit, pad with 1s -> byte should be 0b1111'1111 = 0xFF.
    BitWriter w1(1);
    w1.write_bit(1);
    std::vector<std::uint8_t> b1 = w1.get_bytes();
    REQUIRE(b1.size() == 1);
    CHECK_EQ((int)b1[0], 0xFF);

    // Same but pad with 0s -> 0b1000'0000 = 0x80.
    BitWriter w0(0);
    w0.write_bit(1);
    std::vector<std::uint8_t> b0 = w0.get_bytes();
    REQUIRE(b0.size() == 1);
    CHECK_EQ((int)b0[0], 0x80);
}

TEST_CASE("BitReader signals EOF and counts overshoot via peof") {
    std::uint8_t data[1] = {0xAB};
    BitReader r(data, 1);
    CHECK_EQ((int)r.read(8), 0xAB);
    CHECK(r.eof());
    // Reading past the end returns 0 and accumulates the overshoot in peof.
    CHECK_EQ((int)r.read(5), 0);
    CHECK_EQ(r.peof(), 5);
}

TEST_CASE("BitReader::rewind_bits re-reads bits") {
    std::uint8_t data[2] = {0b10110010, 0b11000000};
    BitReader r(data, 2);
    unsigned int first = r.read(4);       // 0b1011
    CHECK_EQ((int)first, 0b1011);
    r.rewind_bits(4);                     // go back
    CHECK_EQ((int)r.read(4), 0b1011);     // same bits again
}

TEST_CASE("MemoryWriter accumulates bytes and reports counts") {
    MemoryWriter w;
    CHECK_EQ((int)w.num_bytes_written(), 0);
    const std::vector<std::uint8_t> chunk = {1, 2, 3, 4};
    w.write(chunk);
    w.write_byte(5);
    CHECK_EQ((int)w.num_bytes_written(), 5);
    std::vector<std::uint8_t> got = w.get_data();
    const std::vector<std::uint8_t> want = {1, 2, 3, 4, 5};
    CHECK(got == want);
    w.reset();
    CHECK_EQ((int)w.num_bytes_written(), 0);
}

TEST_CASE("MemoryReader reads back what MemoryWriter wrote") {
    const std::vector<std::uint8_t> payload = {9, 8, 7, 6, 5, 4};
    MemoryReader r(payload);
    CHECK_EQ((int)r.get_size(), (int)payload.size());
    for (std::uint8_t b : payload) CHECK_EQ((int)r.read_byte(), (int)b);
    CHECK(r.end_of_reader());
}

TEST_CASE("MemoryReader::read_byte throws at end of data") {
    const std::vector<std::uint8_t> one = {42};
    MemoryReader r(one);
    CHECK_EQ((int)r.read_byte(), 42);
    bool threw = false;
    try { r.read_byte(); } catch (const std::runtime_error&) { threw = true; }
    CHECK(threw);
}

TEST_CASE("MemoryReader skip and rewind_bytes move the cursor") {
    const std::vector<std::uint8_t> payload = {0, 1, 2, 3, 4, 5};
    MemoryReader r(payload);
    r.skip(2);
    CHECK_EQ((int)r.read_byte(), 2);
    CHECK_EQ((int)r.num_bytes_read(), 3);
    r.rewind_bytes(3);
    CHECK_EQ((int)r.read_byte(), 0);
    r.rewind();
    CHECK_EQ((int)r.num_bytes_read(), 0);
}

TEST_CASE("MemoryReader::read into pointer clamps to available bytes") {
    const std::vector<std::uint8_t> payload = {10, 20, 30};
    MemoryReader r(payload);
    std::uint8_t buf[8] = {0};
    std::size_t got = r.read(buf, 8);     // only 3 available
    CHECK_EQ((int)got, 3);
    CHECK_EQ((int)buf[0], 10);
    CHECK_EQ((int)buf[2], 30);
}

#ifndef JPGENCODER_H
#define JPGENCODER_H

#include <cstdint>
#include <map>
#include <memory>
#include <vector>

#include "bitwriter.h"

#include "codingstatus.h"
#include "frameinfo.h"
#include "huffcodes.h"
#include "scaninfo.h"
#include "segment.h"

class JpgEncoder {
public:
	JpgEncoder(FrameInfo& frame_info, const std::vector<Segment>& segments, std::uint8_t padbit);
	// JPEG encoding routine.
	void encode();

	std::vector<std::uint8_t> get_huffman_data() const;
	std::vector<std::size_t> get_restart_marker_pos() const;
	std::vector<std::size_t> get_scan_pos() const;

private:
	// encoding for interleaved data.
	CodingStatus encode_interleaved(int rsti, int& cmp, int& dpos, int& rstw, int& csc, int& mcu, int& sub);
	// encoding for non interleaved data.
	CodingStatus encode_noninterleaved(int rsti, int cmp, int& dpos, int& rstw);

	CodingStatus encode_sequential_noninterleaved(const Component& component, int cmp, int rsti, int& dpos, int& rstw);
	CodingStatus encode_progressive_noninterleaved_dc(const Component& component, int cmp, int rsti, int& dpos, int& rstw);
	CodingStatus encode_progressive_noninterleaved_ac(const Component& component, int rsti, int& dpos, int& rstw);

	CodingStatus encode_sequential_interleaved(int rsti, int& cmp, int& dpos, int& rstw, int& csc, int& mcu, int& sub);
	// Progressive interleaved DC encoding.
	CodingStatus encode_progressive_interleaved_dc(int rsti, int& cmp, int& dpos, int& rstw, int& csc, int& mcu, int& sub);

	// Sequential block encoding routine.
	void block_seq(const HuffCodes& dc_table, const HuffCodes& ac_table);
	// Progressive DC encoding routine.
	void dc_prg_fs(const HuffCodes& dc_table);
	// Progressive AC encoding routine.
	void ac_prg_fs(const HuffCodes& ac_table, int& eobrun);
	// Progressive DC SA encoding routine.
	void dc_prg_sa();
	// Progressive AC SA encoding routine.
	void ac_prg_sa(const HuffCodes& ac_table, int& eobrun);
	// Run of EOB encoding routine.
	void eobrun(const HuffCodes& ac_table, int& eobrun);
	// Correction bits encoding routine.
	void write_correction_bits();

	static constexpr std::int16_t fdiv2(std::int16_t v, int p) {
		return (v < 0) ? -((-v) >> p) : (v >> p);
	}

	static constexpr int envli(int s, int v) {
		return (v > 0) ? v : v - 1 + (1 << s);
	}

	static constexpr int e_envli(int s, int v) {
		return v - (1 << s);
	}

	static constexpr int bitlen(int v) {
		int length = 0;
		while ((v >> length) != 0) {
			length++;
		}
		return length;
	}

	void encode_scan(int restart_interval, int& restart_markers);

	void copy_colldata_to_block_in_scan(const Component& component, int dpos);

	void encode_sequential(const Component& component, int cmp, int dpos);

	// DC successive approximation first stage.
	void dc_succ_approx_first_stage(const Component& component, int cmp, int dpos);
	// DC successive approximation later stage.
	void dc_succ_approx_later_stage(const Component& component, int dpos);

	FrameInfo& frame_info_;

	std::unique_ptr<BitWriter> huffman_writer_; // Bitwise writer for image data.
	std::vector<std::uint8_t> correction_bits_; // Store for correction bits.
	std::array<std::int16_t, 64> block_{};  // Store block for coefficientss.
	ScanInfo scan_info_;

	std::map<int, std::unique_ptr<HuffCodes>> dc_tables_;
	std::map<int, std::unique_ptr<HuffCodes>> ac_tables_;

	std::array<std::int16_t, 4> lastdc_{};  // last dc for each component (used for diff coding)

	const std::vector<Segment>& segments_;

	std::vector<std::uint8_t> huffman_data_;

	std::vector<std::size_t> restart_marker_pos_; // restart markers positions in huffdata
	std::vector<std::size_t> scan_pos_; // scan start positions in huffdata
};

#endif
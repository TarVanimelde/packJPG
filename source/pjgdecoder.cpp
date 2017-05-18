#include "pjgdecoder.h"

#include <algorithm>
#include <string>

#include "bitops.h"
#include "dct8x8.h"
#include "jfif.h"
#include "pjgcontext.h"
#include "pjpgtbl.h"
#include "programinfo.h"

PjgDecoder::PjgDecoder(Reader& decoding_stream) {
	const auto image_pjg_version = decoding_stream.read_byte();
	if (image_pjg_version >= 0x14) {
		if (image_pjg_version != program_info::appversion) {
			throw std::runtime_error("Incompatible file, use " + program_info::appname
				+ " v" + std::to_string(image_pjg_version / 10) + "." + std::to_string(image_pjg_version % 10));
		}
	} else {
		throw std::runtime_error("Unknown pjg version, use newer version of " + program_info::appname);
	}
	decoder_ = std::make_unique<ArithmeticDecoder>(decoding_stream);
}

void PjgDecoder::decode() {
	segments_ = Segment::parse_segments(this->generic());
	padbit_ = this->bit();
	const bool rst_err_used = this->bit() == 1;
	if (rst_err_used) {
		// Decode the number of false set RST markers per scan only if available:
		rst_err_ = this->generic();
	}

	for (auto& segment : segments_) {
		segment.undo_optimize();
	}
	frame_info_ = jfif::get_frame_info(segments_);

	for (auto& component : frame_info_->components) {
		component.freqscan = this->decode_zero_sorted_scan();
		this->zdst_high(component);
		this->ac_high(component);
		this->zdst_low(component);
		this->ac_low(component);
		this->dc(component);
	}

	const bool garbage_exists = this->bit() == 1;
	if (garbage_exists) {
		garbage_data_ = this->generic();
	}
}

std::unique_ptr<FrameInfo> PjgDecoder::get_frame_info() {
	return std::move(frame_info_);
}

std::vector<Segment> PjgDecoder::get_segments() {
	return segments_;
}

std::uint8_t PjgDecoder::get_padbit() {
	return padbit_;
}

std::vector<std::uint8_t> PjgDecoder::get_rst_err() {
	return rst_err_;
}

std::vector<std::uint8_t> PjgDecoder::get_garbage_data() {
	return garbage_data_;
}

std::array<std::uint8_t, 64> PjgDecoder::decode_zero_sorted_scan() {
	std::array<std::uint8_t, 64> zero_sorted_scan{};
	// Skip the first (DC) element, since it is always 0 in the zero-sorted scan order.
	std::vector<std::uint8_t> standard_scan(std::begin(pjg::stdscan) + 1, std::end(pjg::stdscan));

	auto model = std::make_unique<UniversalModel>(64, 64, 1);

	// Decode the zero-sorted scan order:
	for (int i = 1; i < zero_sorted_scan.size(); i++) {
		model->exclude_symbols_above(64 - i);

		int coded_pos = decoder_->decode(*model);
		model->shift_context(coded_pos);

		if (coded_pos == 0) {
			// The remainder of the zero-sorted scan is identical to the standard scan:
			std::copy(std::begin(standard_scan), std::end(standard_scan), std::begin(zero_sorted_scan) + i);
			break;
		}
		coded_pos--;
		zero_sorted_scan[i] = standard_scan[coded_pos];
		standard_scan.erase(std::begin(standard_scan) + coded_pos);
	}

	return zero_sorted_scan;
}

void PjgDecoder::zdst_high(Component& component) {
	// init model, constants
	auto model = std::make_unique<UniversalModel>(49 + 1, 25 + 1, 1);
	auto& zdstls = component.zdstdata;
	const int w = component.bch;
	const int bc = component.bc;

	// arithmetic decode zero-distribution-list
	for (int dpos = 0; dpos < bc; dpos++) {
		// context modelling - use average of above and left as context		
		auto coords = PjgContext::get_context_nnb(dpos, w);
		coords.first = (coords.first >= 0) ? zdstls[coords.first] : 0;
		coords.second = (coords.second >= 0) ? zdstls[coords.second] : 0;
		// shift context
		model->shift_context((coords.first + coords.second + 2) / 4);
		// decode symbol
		zdstls[dpos] = decoder_->decode(*model);
	}
}

void PjgDecoder::zdst_low(Component& component) {
	// init model, constants
	auto model = std::make_unique<UniversalModel>(8, 8, 2);

	auto& zdstls_x = component.zdstxlow;
	auto& zdstls_y = component.zdstylow;

	const auto& ctx_eobx = component.eobxhigh;
	const auto& ctx_eoby = component.eobyhigh;
	const auto& ctx_zdst = component.zdstdata;
	const int bc = component.bc;

	// arithmetic encode zero-distribution-list (first row)
	for (int dpos = 0; dpos < bc; dpos++) {
		model->shift_context((ctx_zdst[dpos] + 3) / 7); // shift context
		model->shift_context(ctx_eobx[dpos]); // shift context
		zdstls_x[dpos] = decoder_->decode(*model); // decode symbol
	}
	// arithmetic encode zero-distribution-list (first collumn)
	for (int dpos = 0; dpos < bc; dpos++) {
		model->shift_context((ctx_zdst[dpos] + 3) / 7); // shift context
		model->shift_context(ctx_eoby[dpos]); // shift context
		zdstls_y[dpos] = decoder_->decode(*model); // decode symbol
	}
}

void PjgDecoder::dc(Component& component) {
	// decide segmentation setting
	const auto& segm_tab = pjg::segm_tables[component.segm_cnt - 1];

	// get max absolute value/bit length
	const int max_val = component.max_v(0); // Max value.
	const int max_len = pjg::bitlen1024p(max_val); // Max bitlength.

	// init models for bitlenghts and -patterns
	auto mod_len = std::make_unique<UniversalModel>(max_len + 1, std::max(component.segm_cnt, max_len + 1), 2);
	auto mod_res = std::make_unique<BinaryModel>(std::max(component.segm_cnt, 16), 2);
	auto mod_sgn = std::make_unique<BinaryModel>(1, 0);

	// set width/height of each band
	const int bc = component.bc;
	const int w = component.bch;

	PjgContext context(component);

	auto& coeffs = component.colldata[0];
	const auto& zero_dist_list = component.zdstdata;

	// arithmetic compression loop
	for (int dpos = 0; dpos < bc; dpos++) {
		//calculate x/y positions in band
		const int p_y = dpos / w;
		// r_y = h - ( p_y + 1 );
		const int p_x = dpos % w;
		const int r_x = w - (p_x + 1);

		// get segment-number from zero distribution list and segmentation set
		const int snum = segm_tab[zero_dist_list[dpos]];
		// calculate contexts (for bit length)
		const int ctx_avr = context.aavrg_context(dpos, p_y, p_x, r_x); // Average context
		const int ctx_len = pjg::bitlen1024p(ctx_avr); // Bitlength context				
		// shift context / do context modelling (segmentation is done per context)
		mod_len->shift_model(ctx_len, snum);
		// decode bit length of current coefficient
		const int clen = decoder_->decode(*mod_len);

		// simple treatment if coefficient is zero
		if (clen == 0) {
			// coeffs[ dpos ] = 0;
		} else {
			// decoding of residual
			int absv = 1;
			// first set bit must be 1, so we start at clen - 2
			for (int bp = clen - 2; bp >= 0; bp--) {
				mod_res->shift_model(snum, bp); // shift in 2 contexts
				// decode bit
				const int bt = decoder_->decode(*mod_res);
				// update absv
				absv = absv << 1;
				if (bt)
					absv |= 1;
			}
			// decode sign
			const int sgn = decoder_->decode(*mod_sgn);
			// copy to colldata
			coeffs[dpos] = (sgn == 0) ? absv : -absv;
			// store absolute value/sign
			context.abs_coeffs_[dpos] = absv;
		}
	}
}

void PjgDecoder::ac_high(Component& component) {
	// decide segmentation setting
	const auto& segm_tab = pjg::segm_tables[component.segm_cnt - 1];

	// init models for bitlenghts and -patterns
	auto mod_len = std::make_unique<UniversalModel>(11, std::max(component.segm_cnt, 11), 2);
	auto mod_res = std::make_unique<BinaryModel>(std::max(component.segm_cnt, 16), 2);
	auto mod_sgn = std::make_unique<BinaryModel>(9, 1);

	// set width/height of each band
	const int bc = component.bc;
	const int w = component.bch;

	std::vector<std::uint8_t> signs(bc); // sign storage for context	
	auto zdstls = component.zdstdata; // copy of zero distribution list

	// set up quick access arrays for signs context
	//std::uint8_t* sgn_nbh = sgn_store.data() - 1; // Left signs neighbor.
	//std::uint8_t* sgn_nbv = sgn_store.data() - w; // Upper signs neighbor.

	auto& eob_x = component.eobxhigh;
	auto& eob_y = component.eobyhigh;

	// set up average context quick access arrays
	PjgContext context(component);

	// work through lower 7x7 bands in order of pjg::freqscan
	for (int i = 1; i < 64; i++) {
		// work through blocks in order of frequency scan
		const int bpos = static_cast<int>(component.freqscan[i]);
		const int b_x = pjg::unzigzag[bpos] % 8;
		const int b_y = pjg::unzigzag[bpos] / 8;

		if ((b_x == 0) || (b_y == 0))
			continue; // process remaining coefficients elsewhere

		// preset absolute values/sign storage
		context.reset_store();
		std::fill(std::begin(signs), std::end(signs), static_cast<std::uint8_t>(0));

		auto& coeffs = component.colldata[bpos]; // Current coefficient data.

		// get max bit length
		const int max_val = component.max_v(bpos); // Max value.
		const int max_len = pjg::bitlen1024p(max_val); // Max bitlength.

		// arithmetic compression loop
		for (int dpos = 0; dpos < bc; dpos++) {
			// skip if beyound eob
			if (zdstls[dpos] == 0)
				continue;

			//calculate x/y positions in band
			const int p_y = dpos / w;
			const int p_x = dpos % w;
			const int r_x = w - (p_x + 1);

			// get segment-number from zero distribution list and segmentation set
			const int snum = segm_tab[zdstls[dpos]];
			// calculate contexts (for bit length)
			const int ctx_avr = context.aavrg_context(dpos, p_y, p_x, r_x); // Average context.
			const int ctx_len = pjg::bitlen1024p(ctx_avr); // Bitlength context.
			// shift context / do context modelling (segmentation is done per context)
			mod_len->shift_model(ctx_len, snum);
			mod_len->exclude_symbols_above(max_len);

			// decode bit length of current coefficient
			const int clen = decoder_->decode(*mod_len);
			// simple treatment if coefficient is zero
			if (clen == 0) {
				// coeffs[ dpos ] = 0;
			} else {
				// decoding of residual
				int absv = 1;
				// first set bit must be 1, so we start at clen - 2
				for (int bp = clen - 2; bp >= 0; bp--) {
					mod_res->shift_model(snum, bp); // shift in 2 contexts
					// decode bit
					const int bt = decoder_->decode(*mod_res);
					// update absv
					absv = absv << 1;
					if (bt)
						absv |= 1;
				}
				// decode sign
				int ctx_sgn = p_x > 0 ? signs[dpos - 1] : 0; // Sign context.
				if (p_y > 0) {
					ctx_sgn += 3 * signs[dpos - w]; // IMPROVE! !!!!!!!!!!!
				}
				mod_sgn->shift_context(ctx_sgn);
				const int sgn = decoder_->decode(*mod_sgn);
				// copy to colldata
				coeffs[dpos] = (sgn == 0) ? absv : -absv;
				// store absolute value/sign, decrement zdst
				context.abs_coeffs_[dpos] = absv;
				signs[dpos] = sgn + 1;
				zdstls[dpos]--;
				// recalculate x/y eob
				if (b_x > eob_x[dpos])
					eob_x[dpos] = b_x;
				if (b_y > eob_y[dpos])
					eob_y[dpos] = b_y;
			}
		}
		// flush models
		mod_len->flush_model();
		mod_res->flush_model();
		mod_sgn->flush_model();
	}
}

void PjgDecoder::ac_low(Component& component) {
	std::array<int16_t*, 8> coeffs_x{nullptr}; // prediction coeffs - current block
	std::array<int16_t*, 8> coeffs_a{nullptr}; // prediction coeffs - neighboring block
	std::array<int, 8> pred_cf{}; // prediction multipliers

	// init models for bitlenghts and -patterns
	auto mod_len = std::make_unique<UniversalModel>(11, std::max(component.segm_cnt, 11), 2);
	auto mod_res = std::make_unique<BinaryModel>(1 << 4, 2);
	auto mod_top = std::make_unique<BinaryModel>(1 << std::max(4, component.nois_trs), 3);
	auto mod_sgn = std::make_unique<BinaryModel>(11, 1);

	// set width/height of each band
	const int bc = component.bc;
	const int w = component.bch;

	// work through each first row / first collumn band
	for (int i = 2; i < 16; i++) {
		// alternate between first row and first collumn
		int b_x = (i % 2 == 0) ? i / 2 : 0;
		int b_y = (i % 2 == 1) ? i / 2 : 0;
		const int bpos = static_cast<int>(pjg::zigzag[b_x + (8 * b_y)]);

		auto& coeffs = component.colldata[bpos]; // Current coefficient data.
		// store pointers to prediction coefficients
		int p_x, p_y;
		int* edge_c; // edge criteria
		auto& zdstls = b_x == 0 ? component.zdstylow : component.zdstxlow; // Pointer to row/col # of non-zeroes.
		if (b_x == 0) {
			for (; b_x < 8; b_x++) {
				coeffs_x[b_x] = component.colldata[pjg::zigzag[b_x + (8 * b_y)]].data();
				coeffs_a[b_x] = component.colldata[pjg::zigzag[b_x + (8 * b_y)]].data() - 1;
				pred_cf[b_x] = dct::icos_base_8x8[b_x * 8] * component.quant(pjg::zigzag[b_x + (8 * b_y)]);
			}
			edge_c = &p_x;
		} else { // if ( b_y == 0 )
			for (; b_y < 8; b_y++) {
				coeffs_x[b_y] = component.colldata[pjg::zigzag[b_x + (8 * b_y)]].data();
				coeffs_a[b_y] = component.colldata[pjg::zigzag[b_x + (8 * b_y)]].data() - w;
				pred_cf[b_y] = dct::icos_base_8x8[b_y * 8] * component.quant(pjg::zigzag[b_x + (8 * b_y)]);
			}
			edge_c = &p_y;
		}

		// get max bit length / other info
		const int max_valp = component.max_v(bpos); // Max value (positive).
		const int max_valn = -max_valp; // Max value (negative).
		const int max_len = pjg::bitlen1024p(max_valp); // Max bitlength.
		const int thrs_bp = std::max(0, max_len - component.nois_trs); // Residual threshold bitplane.

		// arithmetic compression loop
		for (int dpos = 0; dpos < bc; dpos++) {
			// skip if beyound eob
			if (zdstls[dpos] == 0)
				continue;

			//calculate x/y positions in band
			p_y = dpos / w;
			p_x = dpos % w;

			// edge treatment / calculate LAKHANI context
			int ctx_lak; // Lakhani context.
			if ((*edge_c) > 0)
				ctx_lak = PjgContext::lakh_context(coeffs_x, coeffs_a, pred_cf, dpos);
			else
				ctx_lak = 0;
			ctx_lak = bitops::clamp(ctx_lak, max_valn, max_valp);
			const int ctx_len = pjg::bitlen2048n(ctx_lak); // Bitlength context.				
			// shift context / do context modelling (segmentation is done per context)
			mod_len->shift_model(ctx_len, zdstls[dpos]);
			mod_len->exclude_symbols_above(max_len);

			// decode bit length of current coefficient
			const int clen = decoder_->decode(*mod_len);
			// simple treatment if coefficients == 0
			if (clen == 0) {
				// coeffs[ dpos ] = 0;
			} else {
				// decoding of residual
				int bp = clen - 2; // first set bit must be 1, so we start at clen - 2
				int ctx_res = (bp >= thrs_bp) ? 1 : 0; // Bit plane context for residual.
				const int ctx_abs = std::abs(ctx_lak); // Absolute context.
				const int ctx_sgn = (ctx_lak == 0) ? 0 : (ctx_lak > 0) ? 1 : 2; // Context for sign.
				for (; bp >= thrs_bp; bp--) {
					mod_top->shift_model(ctx_abs >> thrs_bp, ctx_res, clen - thrs_bp); // shift in 3 contexts
					// decode bit
					const int bt = decoder_->decode(*mod_top);
					// update context
					ctx_res = ctx_res << 1;
					if (bt)
						ctx_res |= 1;
				}
				int absv = (ctx_res == 0) ? 1 : ctx_res; // !!!!
				for (; bp >= 0; bp--) {
					mod_res->shift_model(zdstls[dpos], bp); // shift in 2 contexts
					// decode bit
					const int bt = decoder_->decode(*mod_res);
					// update absv
					absv = absv << 1;
					if (bt)
						absv |= 1;
				}
				// decode sign
				mod_sgn->shift_model(zdstls[dpos], ctx_sgn);
				const int sgn = decoder_->decode(*mod_sgn);
				// copy to colldata
				coeffs[dpos] = (sgn == 0) ? absv : -absv;
				// decrement # of non zeroes
				zdstls[dpos]--;
			}
		}
		// flush models
		mod_len->flush_model();
		mod_res->flush_model();
		mod_top->flush_model();
		mod_sgn->flush_model();
	}
}

std::vector<std::uint8_t> PjgDecoder::generic() {
	std::vector<std::uint8_t> generic_data;
	auto model = std::make_unique<UniversalModel>(256 + 1, 256, 1);
	for (int c = decoder_->decode(*model); c != 256; c = decoder_->decode(*model)) {
		generic_data.emplace_back(static_cast<std::uint8_t>(c));
		model->shift_context(c);
	}

	return generic_data;
}

std::uint8_t PjgDecoder::bit() {
	auto model = std::make_unique<BinaryModel>(1, -1);
	std::uint8_t bit = decoder_->decode(*model); // This conversion is okay since there are only 2 symbols in the model.
	return bit;
}
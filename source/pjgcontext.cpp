#include "pjgcontext.h"
#include "pjpgtbl.h"

// Context weighting factors:
namespace pjg {
constexpr std::array<int, 6> weights{
	pjg::abs_ctx_weights_lum[0][0][2], // top-top
	pjg::abs_ctx_weights_lum[0][1][1], // top-left
	pjg::abs_ctx_weights_lum[0][1][2], // top
	pjg::abs_ctx_weights_lum[0][1][3], // top-right
	pjg::abs_ctx_weights_lum[0][2][0], // left-left
	pjg::abs_ctx_weights_lum[0][2][1] // left
};
}

PjgContext::PjgContext(const Component& component) : abs_coeffs_(component.bc) {
	const auto w = component.bch;
	quick_abs_coeffs_[0] = abs_coeffs_.data() + (0 + ((-2) * w)); // top-top
	quick_abs_coeffs_[1] = abs_coeffs_.data() + (-1 + ((-1) * w)); // top-left
	quick_abs_coeffs_[2] = abs_coeffs_.data() + (0 + ((-1) * w)); // top
	quick_abs_coeffs_[3] = abs_coeffs_.data() + (1 + ((-1) * w)); // top-right
	quick_abs_coeffs_[4] = abs_coeffs_.data() + (-2 + ((0) * w)); // left-left
	quick_abs_coeffs_[5] = abs_coeffs_.data() + (-1 + ((0) * w)); // left
}

void PjgContext::reset_store() {
	std::fill(std::begin(abs_coeffs_), std::end(abs_coeffs_), static_cast<std::uint16_t>(0));
}

int PjgContext::aavrg_context(int pos, int band_width) {
	// Calculate x/y positions in band:
	const int p_y = pos / band_width;
	const int p_x = pos % band_width;
	const int r_x = band_width - (p_x + 1);

	int average_context = 0;
	int tot_context_weight = 0;
	int curr_context_weight;

	// different cases due to edge treatment
	if (p_y >= 2) {
		curr_context_weight = std::get<0>(pjg::weights);
		average_context += quick_abs_coeffs_[0][pos] * curr_context_weight;
		tot_context_weight += curr_context_weight;
		curr_context_weight = std::get<2>(pjg::weights);
		average_context += quick_abs_coeffs_[2][pos] * curr_context_weight;
		tot_context_weight += curr_context_weight;
		if (p_x >= 2) {
			curr_context_weight = std::get<1>(pjg::weights);
			average_context += quick_abs_coeffs_[1][pos] * curr_context_weight;
			tot_context_weight += curr_context_weight;
			curr_context_weight = std::get<4>(pjg::weights);
			average_context += quick_abs_coeffs_[4][pos] * curr_context_weight;
			tot_context_weight += curr_context_weight;
			curr_context_weight = std::get<5>(pjg::weights);
			average_context += quick_abs_coeffs_[5][pos] * curr_context_weight;
			tot_context_weight += curr_context_weight;
		} else if (p_x == 1) {
			curr_context_weight = std::get<1>(pjg::weights);
			average_context += quick_abs_coeffs_[1][pos] * curr_context_weight;
			tot_context_weight += curr_context_weight;
			curr_context_weight = std::get<5>(pjg::weights);
			average_context += quick_abs_coeffs_[5][pos] * curr_context_weight;
			tot_context_weight += curr_context_weight;
		}
		if (r_x >= 1) {
			curr_context_weight = std::get<3>(pjg::weights);
			average_context += quick_abs_coeffs_[3][pos] * curr_context_weight;
			tot_context_weight += curr_context_weight;
		}
	} else if (p_y == 1) {
		curr_context_weight = std::get<2>(pjg::weights);
		average_context += quick_abs_coeffs_[2][pos] * curr_context_weight;
		tot_context_weight += curr_context_weight;
		if (p_x >= 2) {
			curr_context_weight = std::get<1>(pjg::weights);
			average_context += quick_abs_coeffs_[1][pos] * curr_context_weight;
			tot_context_weight += curr_context_weight;
			curr_context_weight = std::get<4>(pjg::weights);
			average_context += quick_abs_coeffs_[4][pos] * curr_context_weight;
			tot_context_weight += curr_context_weight;
			curr_context_weight = std::get<5>(pjg::weights);
			average_context += quick_abs_coeffs_[5][pos] * curr_context_weight;
			tot_context_weight += curr_context_weight;
		} else if (p_x == 1) {
			curr_context_weight = std::get<1>(pjg::weights);
			average_context += quick_abs_coeffs_[1][pos] * curr_context_weight;
			tot_context_weight += curr_context_weight;
			curr_context_weight = std::get<5>(pjg::weights);
			average_context += quick_abs_coeffs_[5][pos] * curr_context_weight;
			tot_context_weight += curr_context_weight;
		}
		if (r_x >= 1) {
			curr_context_weight = std::get<3>(pjg::weights);
			average_context += quick_abs_coeffs_[3][pos] * curr_context_weight;
			tot_context_weight += curr_context_weight;
		}
	} else {
		if (p_x >= 2) {
			curr_context_weight = std::get<4>(pjg::weights);
			average_context += quick_abs_coeffs_[4][pos] * curr_context_weight;
			tot_context_weight += curr_context_weight;
			curr_context_weight = std::get<5>(pjg::weights);
			average_context += quick_abs_coeffs_[5][pos] * curr_context_weight;
			tot_context_weight += curr_context_weight;
		} else if (p_x == 1) {
			curr_context_weight = std::get<5>(pjg::weights);
			average_context += quick_abs_coeffs_[5][pos] * curr_context_weight;
			tot_context_weight += curr_context_weight;
		}
	}

	// return average context
	return (tot_context_weight != 0) ? (average_context + (tot_context_weight / 2)) / tot_context_weight : 0;
}

int PjgContext::lakh_context(const std::array<int16_t*, 8>& coeffs_x, const std::array<int16_t*, 8>& coeffs_a, const std::array<int, 8>& pred_cf, int pos) {
	int pred = 0;

	// calculate partial prediction
	pred -= (coeffs_x[1][pos] + coeffs_a[1][pos]) * pred_cf[1];
	pred -= (coeffs_x[2][pos] - coeffs_a[2][pos]) * pred_cf[2];
	pred -= (coeffs_x[3][pos] + coeffs_a[3][pos]) * pred_cf[3];
	pred -= (coeffs_x[4][pos] - coeffs_a[4][pos]) * pred_cf[4];
	pred -= (coeffs_x[5][pos] + coeffs_a[5][pos]) * pred_cf[5];
	pred -= (coeffs_x[6][pos] - coeffs_a[6][pos]) * pred_cf[6];
	pred -= (coeffs_x[7][pos] + coeffs_a[7][pos]) * pred_cf[7];
	// normalize / quantize partial prediction
	pred = ((pred > 0) ? (pred + (pred_cf[0] / 2)) : (pred - (pred_cf[0] / 2))) / pred_cf[0];
	// complete prediction
	pred += coeffs_a[0][pos];

	return pred;
}

std::pair<int, int> PjgContext::get_context_nnb(int pos, int w) {
	std::pair<int, int> coords;
	if (pos == 0) {
		coords = std::make_pair<int, int>(-1, -1);
	} else if ((pos % w) == 0) {
		if (pos >= (w << 1)) {
			coords = std::make_pair<int, int>(pos - (w << 1), pos - w);
		} else {
			coords = std::make_pair<int, int>(pos - w, pos - w);
		}
	} else if (pos < w) {
		if (pos >= 2) {
			coords = std::make_pair<int, int>(pos - 1, pos - 2);
		} else {
			coords = std::make_pair<int, int>(pos - 1, pos - 1);
		}
	} else {
		coords = std::make_pair<int, int>(pos - 1, pos - w);
	}
	return coords;
}
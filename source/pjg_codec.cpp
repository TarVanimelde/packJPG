#include "pjg_codec.h"

#include "jpeg_codec.h"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <limits>
#include <stdexcept>

#if defined(__GNUC__) || defined(__clang__)
#define PJG_RESTRICT __restrict__
#elif defined(_MSC_VER)
#define PJG_RESTRICT __restrict
#else
#define PJG_RESTRICT
#endif

namespace {

constexpr std::size_t kMaxGenericStreamBytes      = 64u << 20; // 64 MiB
constexpr std::size_t kMaxGenericExpansionPerByte = 2048u;     // > 1411:1 theoretical max of c_lim=255
constexpr std::size_t kMaxRstErrBytes             = 65536u;

inline int skip_zero_zdst( const unsigned char* zdstls, int px, int row_w ) noexcept
{
	if ( ( px < row_w ) && ( zdstls[ px ] != 0 ) ) {
		return px;
	}
	while ( px + 8 <= row_w ) {
		std::uint64_t z8 = 0;
		std::memcpy( &z8, zdstls + px, 8 );
		if ( z8 != 0 ) break;
		px += 8;
	}
	while ( ( px < row_w ) && ( zdstls[ px ] == 0 ) ) {
		++px;
	}
	return px;
}

void get_context_nnb( int pos, int w, int* a, int* b )
{
	if ( pos == 0 ) {
		*a = -1;
		*b = -1;
	}
	else if ( ( pos % w ) == 0 ) {
		*b = pos - w;
		*a = ( pos >= ( w << 1 ) ) ? ( pos - ( w << 1 ) ) : *b;
	}
	else if ( pos < w ) {
		*a = pos - 1;
		*b = ( pos >= 2 ) ? ( pos - 2 ) : *a;
	}
	else {
		*a = pos - 1;
		*b = pos - w;
	}
}

void pjg_aavrg_prepare( const CodecContext& ctx, int* offsets, int* weights, int cmp )
{
	const int w = ctx.cmpnfo[ cmp ].bch;

	offsets[ 0 ] =  0 + ( ( -2 ) * w ); // top-top
	offsets[ 1 ] = -1 + ( ( -1 ) * w ); // top-left
	offsets[ 2 ] =  0 + ( ( -1 ) * w ); // top
	offsets[ 3 ] =  1 + ( ( -1 ) * w ); // top-right
	offsets[ 4 ] = -2 + ( (  0 ) * w ); // left-left
	offsets[ 5 ] = -1 + ( (  0 ) * w ); // left

	weights[ 0 ] = abs_ctx_weights_lum[ 0 ][ 0 ][ 2 ];
	weights[ 1 ] = abs_ctx_weights_lum[ 0 ][ 1 ][ 1 ];
	weights[ 2 ] = abs_ctx_weights_lum[ 0 ][ 1 ][ 2 ];
	weights[ 3 ] = abs_ctx_weights_lum[ 0 ][ 1 ][ 3 ];
	weights[ 4 ] = abs_ctx_weights_lum[ 0 ][ 2 ][ 0 ];
	weights[ 5 ] = abs_ctx_weights_lum[ 0 ][ 2 ][ 1 ];
}

inline int pjg_aavrg_context( const unsigned short* abs_store, const int* offsets, const int* weights, int pos, int p_y, int p_x, int r_x )
{
	constexpr int kW0   = abs_ctx_weights_lum[ 0 ][ 0 ][ 2 ];
	constexpr int kW1   = abs_ctx_weights_lum[ 0 ][ 1 ][ 1 ];
	constexpr int kW2   = abs_ctx_weights_lum[ 0 ][ 1 ][ 2 ];
	constexpr int kW3   = abs_ctx_weights_lum[ 0 ][ 1 ][ 3 ];
	constexpr int kW4   = abs_ctx_weights_lum[ 0 ][ 2 ][ 0 ];
	constexpr int kW5   = abs_ctx_weights_lum[ 0 ][ 2 ][ 1 ];
	constexpr unsigned int kWAll = static_cast<unsigned int>( kW0 + kW1 + kW2 + kW3 + kW4 + kW5 );

	if ( ( p_y >= 2 ) && ( p_x >= 2 ) && ( r_x >= 1 ) ) {
		const unsigned int ctx_avr =
			static_cast<unsigned int>( abs_store[ pos + offsets[ 0 ] ] ) * kW0 +
			static_cast<unsigned int>( abs_store[ pos + offsets[ 1 ] ] ) * kW1 +
			static_cast<unsigned int>( abs_store[ pos + offsets[ 2 ] ] ) * kW2 +
			static_cast<unsigned int>( abs_store[ pos + offsets[ 3 ] ] ) * kW3 +
			static_cast<unsigned int>( abs_store[ pos + offsets[ 4 ] ] ) * kW4 +
			static_cast<unsigned int>( abs_store[ pos + offsets[ 5 ] ] ) * kW5;
		return static_cast<int>( ( ctx_avr + ( kWAll / 2u ) ) / kWAll );
	}

	int ctx_avr = 0;
	int w_ctx   = 0;
	int w_curr  = 0;

	if ( p_y >= 2 ) {
		w_curr = weights[ 0 ]; ctx_avr += abs_store[ pos + offsets[ 0 ] ] * w_curr; w_ctx += w_curr;
		w_curr = weights[ 2 ]; ctx_avr += abs_store[ pos + offsets[ 2 ] ] * w_curr; w_ctx += w_curr;
		if ( p_x >= 2 ) {
			w_curr = weights[ 1 ]; ctx_avr += abs_store[ pos + offsets[ 1 ] ] * w_curr; w_ctx += w_curr;
			w_curr = weights[ 4 ]; ctx_avr += abs_store[ pos + offsets[ 4 ] ] * w_curr; w_ctx += w_curr;
			w_curr = weights[ 5 ]; ctx_avr += abs_store[ pos + offsets[ 5 ] ] * w_curr; w_ctx += w_curr;
		}
		else if ( p_x == 1 ) {
			w_curr = weights[ 1 ]; ctx_avr += abs_store[ pos + offsets[ 1 ] ] * w_curr; w_ctx += w_curr;
			w_curr = weights[ 5 ]; ctx_avr += abs_store[ pos + offsets[ 5 ] ] * w_curr; w_ctx += w_curr;
		}
		if ( r_x >= 1 ) {
			w_curr = weights[ 3 ]; ctx_avr += abs_store[ pos + offsets[ 3 ] ] * w_curr; w_ctx += w_curr;
		}
	}
	else if ( p_y == 1 ) {
		w_curr = weights[ 2 ]; ctx_avr += abs_store[ pos + offsets[ 2 ] ] * w_curr; w_ctx += w_curr;
		if ( p_x >= 2 ) {
			w_curr = weights[ 1 ]; ctx_avr += abs_store[ pos + offsets[ 1 ] ] * w_curr; w_ctx += w_curr;
			w_curr = weights[ 4 ]; ctx_avr += abs_store[ pos + offsets[ 4 ] ] * w_curr; w_ctx += w_curr;
			w_curr = weights[ 5 ]; ctx_avr += abs_store[ pos + offsets[ 5 ] ] * w_curr; w_ctx += w_curr;
		}
		else if ( p_x == 1 ) {
			w_curr = weights[ 1 ]; ctx_avr += abs_store[ pos + offsets[ 1 ] ] * w_curr; w_ctx += w_curr;
			w_curr = weights[ 5 ]; ctx_avr += abs_store[ pos + offsets[ 5 ] ] * w_curr; w_ctx += w_curr;
		}
		if ( r_x >= 1 ) {
			w_curr = weights[ 3 ]; ctx_avr += abs_store[ pos + offsets[ 3 ] ] * w_curr; w_ctx += w_curr;
		}
	}
	else {
		if ( p_x >= 2 ) {
			w_curr = weights[ 4 ]; ctx_avr += abs_store[ pos + offsets[ 4 ] ] * w_curr; w_ctx += w_curr;
			w_curr = weights[ 5 ]; ctx_avr += abs_store[ pos + offsets[ 5 ] ] * w_curr; w_ctx += w_curr;
		}
		else if ( p_x == 1 ) {
			w_curr = weights[ 5 ]; ctx_avr += abs_store[ pos + offsets[ 5 ] ] * w_curr; w_ctx += w_curr;
		}
	}

	return ( w_ctx != 0 ) ? ( ctx_avr + ( w_ctx / 2 ) ) / w_ctx : 0;
}

struct LakhBandSpec {
	const short* coeffs_x[ 8 ]{};
	int pred_cf[ 8 ]{};
	unsigned char* zdstls = nullptr;
	int neighbor_off = 0;
	bool check_x = false;
	bool fits_i32 = false;
};

LakhBandSpec setup_lakh_band( CodecContext& ctx, int cmp, int w, int b_x, int b_y ) noexcept
{
	LakhBandSpec s{};
	s.check_x        = ( b_x == 0 );
	s.neighbor_off   = s.check_x ? 1 : w;
	s.zdstls         = s.check_x ? ctx.zdstylow[ cmp ].data() : ctx.zdstxlow[ cmp ].data();
	std::int64_t max_numer = 0;
	int max_c0       = 0;
	for ( int k = 0; k < 8; k++ ) {
		const int zpos    = s.check_x ? zigzag[ k + ( 8 * b_y ) ] : zigzag[ b_x + ( 8 * k ) ];
		s.coeffs_x[ k ]   = ctx.colldata[ cmp ][ zpos ].data();
		const int q       = ctx.quant( cmp, zpos );
		s.pred_cf[ k ]    = icos_base_8x8[ static_cast<std::size_t>( k ) * 8u ] * q;
		const int max_mag = ( 1 << bitlen_1024p( ctx.max_v( cmp, zpos ) ) ) - 1;
		if ( k == 0 ) {
			max_c0 = max_mag;
		}
		else {
			max_numer += 2LL * max_mag * s.pred_cf[ k ];
		}
	}
	max_numer += s.pred_cf[ 0 ] / 2;
	s.fits_i32 = ( max_numer + max_c0 <= std::numeric_limits<int>::max() );
	return s;
}

inline int pjg_lakh_context( const LakhBandSpec& s, int pos ) noexcept
{
	const int npos = pos - s.neighbor_off;
	if ( s.fits_i32 ) {
		int pred = 0;
		pred -= ( static_cast<int>( s.coeffs_x[ 1 ][ pos ] ) + static_cast<int>( s.coeffs_x[ 1 ][ npos ] ) ) * s.pred_cf[ 1 ];
		pred -= ( static_cast<int>( s.coeffs_x[ 2 ][ pos ] ) - static_cast<int>( s.coeffs_x[ 2 ][ npos ] ) ) * s.pred_cf[ 2 ];
		pred -= ( static_cast<int>( s.coeffs_x[ 3 ][ pos ] ) + static_cast<int>( s.coeffs_x[ 3 ][ npos ] ) ) * s.pred_cf[ 3 ];
		pred -= ( static_cast<int>( s.coeffs_x[ 4 ][ pos ] ) - static_cast<int>( s.coeffs_x[ 4 ][ npos ] ) ) * s.pred_cf[ 4 ];
		pred -= ( static_cast<int>( s.coeffs_x[ 5 ][ pos ] ) + static_cast<int>( s.coeffs_x[ 5 ][ npos ] ) ) * s.pred_cf[ 5 ];
		pred -= ( static_cast<int>( s.coeffs_x[ 6 ][ pos ] ) - static_cast<int>( s.coeffs_x[ 6 ][ npos ] ) ) * s.pred_cf[ 6 ];
		pred -= ( static_cast<int>( s.coeffs_x[ 7 ][ pos ] ) + static_cast<int>( s.coeffs_x[ 7 ][ npos ] ) ) * s.pred_cf[ 7 ];
		if ( s.pred_cf[ 0 ] != 0 ) {
			const int half = s.pred_cf[ 0 ] / 2;
			pred = ( ( pred > 0 ) ? ( pred + half ) : ( pred - half ) ) / s.pred_cf[ 0 ];
		}
		else {
			pred = 0;
		}
		return pred + static_cast<int>( s.coeffs_x[ 0 ][ npos ] );
	}

	std::int64_t pred = 0;
	pred -= static_cast<std::int64_t>( static_cast<int>( s.coeffs_x[ 1 ][ pos ] ) + static_cast<int>( s.coeffs_x[ 1 ][ npos ] ) ) * s.pred_cf[ 1 ];
	pred -= static_cast<std::int64_t>( static_cast<int>( s.coeffs_x[ 2 ][ pos ] ) - static_cast<int>( s.coeffs_x[ 2 ][ npos ] ) ) * s.pred_cf[ 2 ];
	pred -= static_cast<std::int64_t>( static_cast<int>( s.coeffs_x[ 3 ][ pos ] ) + static_cast<int>( s.coeffs_x[ 3 ][ npos ] ) ) * s.pred_cf[ 3 ];
	pred -= static_cast<std::int64_t>( static_cast<int>( s.coeffs_x[ 4 ][ pos ] ) - static_cast<int>( s.coeffs_x[ 4 ][ npos ] ) ) * s.pred_cf[ 4 ];
	pred -= static_cast<std::int64_t>( static_cast<int>( s.coeffs_x[ 5 ][ pos ] ) + static_cast<int>( s.coeffs_x[ 5 ][ npos ] ) ) * s.pred_cf[ 5 ];
	pred -= static_cast<std::int64_t>( static_cast<int>( s.coeffs_x[ 6 ][ pos ] ) - static_cast<int>( s.coeffs_x[ 6 ][ npos ] ) ) * s.pred_cf[ 6 ];
	pred -= static_cast<std::int64_t>( static_cast<int>( s.coeffs_x[ 7 ][ pos ] ) + static_cast<int>( s.coeffs_x[ 7 ][ npos ] ) ) * s.pred_cf[ 7 ];
	if ( s.pred_cf[ 0 ] != 0 ) {
		const std::int64_t half = s.pred_cf[ 0 ] / 2;
		pred = ( ( pred > 0 ) ? ( pred + half ) : ( pred - half ) ) / s.pred_cf[ 0 ];
	}
	else {
		pred = 0;
	}
	pred += s.coeffs_x[ 0 ][ npos ];

	if ( pred < std::numeric_limits<int>::min() ) return std::numeric_limits<int>::min();
	if ( pred > std::numeric_limits<int>::max() ) return std::numeric_limits<int>::max();
	return static_cast<int>( pred );
}


void pjg_get_zerosort_scan( const CodecContext& ctx, unsigned char* sv, int cmp )
{
	unsigned int zdist[ 64 ];
	const int bc = ctx.cmpnfo[ cmp ].bc;

	for ( int i = 0; i < 64; i++ ) {
		sv[ i ]    = static_cast<unsigned char>( i );
		zdist[ i ] = 0;
	}

	for ( int bpos = 0; bpos < 64; bpos++ ) {
		const short* PJG_RESTRICT band = ctx.colldata[ cmp ][ bpos ].data();
		unsigned int zc = 0;
		for ( int dpos = 0; dpos < bc; ) {
			const int end = std::min( bc, dpos + 256 );
			unsigned short sub = 0;
			for ( int i = dpos; i < end; i++ ) {
				sub = static_cast<unsigned short>( sub + ( band[ i ] == 0 ) );
			}
			zc += sub;
			dpos = end;
		}
		zdist[ bpos ] = zc;
	}

	bool done = false;
	while ( !done ) {
		done = true;
		for ( int i = 2; i < 64; i++ ) {
			if ( zdist[ i ] < zdist[ i - 1 ] ) {
				std::swap( zdist[ i ], zdist[ i - 1 ] );
				std::swap( sv[ i ], sv[ i - 1 ] );
				done = false;
			}
		}
	}
}

bool pjg_optimize_header( CodecContext& ctx )
{
	unsigned char type = 0x00;
	unsigned int  len  = 0;
	unsigned int  hpos = 0;

	while ( hpos < ctx.hdrdata.size() ) {
		if ( !ctx.read_header_segment( hpos, type, len ) ) return false;
		if ( type == 0xC4 ) {
			const unsigned int fpos = hpos + len;
			hpos += 4;
			while ( hpos < fpos ) {
				if ( hpos + 17 > fpos ) {
					ctx.errormessage = "size mismatch in dht marker";
					ctx.errorlevel = 2;
					return false;
				}
				hpos++;
				int sym_count = 0;
				if ( !jpg_parse_dht_counts( &ctx.hdrdata[ hpos ], sym_count ) ) {
					ctx.errormessage = "size mismatch in dht marker";
					ctx.errorlevel = 2;
					return false;
				}
				const unsigned int skip = 16u + static_cast<unsigned int>( sym_count );
				if ( hpos + skip > fpos ) {
					ctx.errormessage = "size mismatch in dht marker";
					ctx.errorlevel = 2;
					return false;
				}
				for ( int i = 0; i < 4; i++ ) {
					if ( hpos + std_huff_lengths[ i ] > fpos )
						continue;
					unsigned int spos = 0;
					for ( ; spos < std_huff_lengths[ i ]; spos++ ) {
						if ( ctx.hdrdata[ hpos + spos ] != std_huff_tables[ i ][ spos ] )
							break;
					}
					if ( spos != std_huff_lengths[ i ] )
						continue;

					ctx.hdrdata[ hpos + 0 ] = static_cast<unsigned char>( std_huff_lengths[ i ] - 16 - i );
					ctx.hdrdata[ hpos + 1 ] = static_cast<unsigned char>( i );
					for ( spos = 2; spos < std_huff_lengths[ i ]; spos++ )
						ctx.hdrdata[ hpos + spos ] = 0x00;
					break;
				}
				hpos += skip;
			}
		}
		else if ( type == 0xDB ) {
			const unsigned int fpos = hpos + len;
			hpos += 4;
			while ( hpos < fpos ) {
				const int prec = static_cast<int>( lbits_u8( ctx.hdrdata[ hpos ], 4 ) );
				hpos++;
				if ( prec == 1 ) {
					if ( hpos + 128 > fpos ) {
						ctx.errormessage = "size mismatch in dqt marker";
						ctx.errorlevel = 2;
						return false;
					}
					hpos += 128;
					continue;
				}
				if ( hpos + 64 > fpos ) {
					ctx.errormessage = "size mismatch in dqt marker";
					ctx.errorlevel = 2;
					return false;
				}
				for ( unsigned int spos = 63; spos > 0; spos-- )
					ctx.hdrdata[ hpos + spos ] -= ctx.hdrdata[ hpos + spos - 1 ];
				hpos += 64;
			}
		}
		else {
			hpos += len;
		}
	}

	return true;
}

bool pjg_unoptimize_header( CodecContext& ctx )
{
	unsigned char type = 0x00;
	unsigned int  len  = 0;
	unsigned int  hpos = 0;

	while ( hpos < ctx.hdrdata.size() ) {
		if ( !ctx.read_header_segment( hpos, type, len ) ) return false;

		if ( type == 0xC4 ) {
			const unsigned int fpos = hpos + len;
			hpos += 4;
			while ( hpos < fpos ) {
				if ( hpos + 17 > fpos ) {
					ctx.errormessage = "size mismatch in dht marker";
					ctx.errorlevel = 2;
					return false;
				}
				hpos++;
				if ( ctx.hdrdata[ hpos ] > 2 ) {
					const int i = ctx.hdrdata[ hpos + 1 ];
					if ( ( i < 0 ) || ( i >= 4 ) || ( hpos + std_huff_lengths[ i ] > fpos ) ) {
						ctx.errormessage = "size mismatch in dht marker";
						ctx.errorlevel = 2;
						return false;
					}
					for ( unsigned int spos = 0; spos < std_huff_lengths[ i ]; spos++ ) {
						ctx.hdrdata[ hpos + spos ] = std_huff_tables[ i ][ spos ];
					}
				}

				unsigned int skip = 16;
				for ( int i = 0; i < 16; i++ )
					skip += static_cast<unsigned int>( ctx.hdrdata[ hpos + i ] );
				if ( hpos + skip > fpos ) {
					ctx.errormessage = "size mismatch in dht marker";
					ctx.errorlevel = 2;
					return false;
				}
				hpos += skip;
			}
		}
		else if ( type == 0xDB ) {
			const unsigned int fpos = hpos + len;
			hpos += 4;
			while ( hpos < fpos ) {
				const int prec = static_cast<int>( lbits_u8( ctx.hdrdata[ hpos ], 4 ) );
				hpos++;
				if ( prec == 1 ) {
					if ( hpos + 128 > fpos ) {
						ctx.errormessage = "size mismatch in dqt marker";
						ctx.errorlevel = 2;
						return false;
					}
					hpos += 128;
					continue;
				}
				if ( hpos + 64 > fpos ) {
					ctx.errormessage = "size mismatch in dqt marker";
					ctx.errorlevel = 2;
					return false;
				}
				for ( unsigned int spos = 1; spos < 64; spos++ )
					ctx.hdrdata[ hpos + spos ] += ctx.hdrdata[ hpos + spos - 1 ];
				hpos += 64;
			}
		}
		else {
			hpos += len;
		}
	}

	return true;
}

bool pjg_encode_zstscan( CodecContext& ctx, ArithmeticEncoder* enc, int cmp )
{
	unsigned char freqlist[ 64 ];

	pjg_get_zerosort_scan( ctx, ctx.zsrtscan[ cmp ], cmp );

	for ( int i = 0; i < 64; i++ )
		freqlist[ i ] = stdscan[ i ];

	auto model = init_model_s( 64, 64, 1 );

	for ( int i = 1; i < 64; i++ ) {
		model->exclude_symbols( 64 - i );

		int tpos = 0;
		int c = i;
		for ( ; c < 64; c++ ) {
			for ( tpos++; freqlist[ tpos ] == 0; tpos++ );
			if ( freqlist[ tpos ] != ctx.zsrtscan[ cmp ][ c ] ) break;
		}
		if ( c == 64 ) {
			encode_ari( enc, model, 0 );
			break;
		}

		int cpos = 1;
		for ( tpos = 0; freqlist[ tpos ] != ctx.zsrtscan[ cmp ][ i ]; tpos++ )
			if ( freqlist[ tpos ] != 0 ) cpos++;
		freqlist[ tpos ] = 0;

		encode_ari( enc, model, cpos );
		model->shift_context( cpos );
	}

	ctx.freqscan[ cmp ] = ctx.zsrtscan[ cmp ];
	return true;
}

bool pjg_encode_zdst_high( const CodecContext& ctx, ArithmeticEncoder* enc, int cmp )
{
	auto model = init_model_s( 49 + 1, 25 + 1, 1 );
	const unsigned char* zdstls = ctx.zdstdata[ cmp ].data();
	const int w  = ctx.cmpnfo[ cmp ].bch;
	const int bc = ctx.cmpnfo[ cmp ].bc;

	for ( int dpos = 0; dpos < bc; dpos++ ) {
		int a = 0;
		int b = 0;
		get_context_nnb( dpos, w, &a, &b );
		a = ( a >= 0 ) ? zdstls[ a ] : 0;
		b = ( b >= 0 ) ? zdstls[ b ] : 0;
		model->shift_context( ( a + b + 2 ) / 4 );
		encode_ari( enc, model, zdstls[ dpos ] );
	}

	return true;
}

bool pjg_encode_zdst_low( const CodecContext& ctx, ArithmeticEncoder* enc, int cmp )
{
	auto model = init_model_s( 8, 8, 2 );
	const unsigned char* zdstls_x = ctx.zdstxlow[ cmp ].data();
	const unsigned char* zdstls_y = ctx.zdstylow[ cmp ].data();
	const unsigned char* ctx_eobx = ctx.eobxhigh[ cmp ].data();
	const unsigned char* ctx_eoby = ctx.eobyhigh[ cmp ].data();
	const unsigned char* ctx_zdst = ctx.zdstdata[ cmp ].data();
	const int bc = ctx.cmpnfo[ cmp ].bc;

	for ( int dpos = 0; dpos < bc; dpos++ ) {
		model->shift_context( ( ctx_zdst[ dpos ] + 3 ) / 7 );
		model->shift_context( ctx_eobx[ dpos ] );
		encode_ari( enc, model, zdstls_x[ dpos ] );
	}
	for ( int dpos = 0; dpos < bc; dpos++ ) {
		model->shift_context( ( ctx_zdst[ dpos ] + 3 ) / 7 );
		model->shift_context( ctx_eoby[ dpos ] );
		encode_ari( enc, model, zdstls_y[ dpos ] );
	}

	return true;
}

bool pjg_encode_dc( const CodecContext& ctx, ArithmeticEncoder* enc, int cmp )
{
	const unsigned char* segm_tab = segm_tables[ ctx.segm_cnt[ cmp ] - 1 ];

	const int max_val = ctx.max_v( cmp, 0 );
	const int max_len = bitlen_1024p( max_val );

	auto mod_len = init_model_s( max_len + 1, ( ctx.segm_cnt[ cmp ] > max_len ) ? ctx.segm_cnt[ cmp ] : max_len + 1, 2 );
	auto mod_res = init_model_b( ( ctx.segm_cnt[ cmp ] < 16 ) ? ( 1 << 4 ) : ctx.segm_cnt[ cmp ], 2 );
	auto mod_sgn = init_model_b( 1, 0 );

	const int bc = ctx.cmpnfo[ cmp ].bc;
	const int w  = ctx.cmpnfo[ cmp ].bch;

	std::vector<unsigned short> absv_store_vec( static_cast<std::size_t>( bc ), 0 );
	unsigned short* absv_store = absv_store_vec.data();

	int c_absc_off[ 6 ];
	int c_weight[ 6 ];
	pjg_aavrg_prepare( ctx, c_absc_off, c_weight, cmp );

	const short* coeffs = ctx.colldata[ cmp ][ 0 ].data();
	const unsigned char* zdstls = ctx.zdstdata[ cmp ].data();

	for ( int dpos = 0, p_y = 0; dpos < bc; ++p_y ) {
		const int row_w = std::min( w, bc - dpos );
		for ( int p_x = 0; p_x < row_w; ++p_x, ++dpos ) {
			const int r_x = w - ( p_x + 1 );

			const int snum    = segm_tab[ zdstls[ dpos ] ];
			const int ctx_avr = pjg_aavrg_context( absv_store, c_absc_off, c_weight, dpos, p_y, p_x, r_x );
			const int ctx_len = bitlen_1024p( ctx_avr );
			shift_model( mod_len, ctx_len, snum );

			if ( coeffs[ dpos ] == 0 ) {
				encode_ari( enc, mod_len, 0 );
			}
			else {
				const int absv = abs_val( coeffs[ dpos ] );
				const int clen = bitlen_1024p( absv );
				const int sgn  = ( coeffs[ dpos ] > 0 ) ? 0 : 1;
				encode_ari( enc, mod_len, clen );
				for ( int bp = clen - 2; bp >= 0; bp-- ) {
					shift_model( mod_res, snum, bp );
					const int bt = bit_n( absv, bp );
					encode_ari( enc, mod_res, bt );
				}
				encode_ari( enc, mod_sgn, sgn );
				absv_store[ dpos ] = static_cast<unsigned short>( absv );
			}
		}
	}

	return true;
}

bool pjg_encode_ac_high( CodecContext& ctx, ArithmeticEncoder* enc, int cmp )
{
	const unsigned char* segm_tab = segm_tables[ ctx.segm_cnt[ cmp ] - 1 ];

	auto mod_len = init_model_s( 11, ( ctx.segm_cnt[ cmp ] > 11 ) ? ctx.segm_cnt[ cmp ] : 11, 2 );
	auto mod_res = init_model_b( ( ctx.segm_cnt[ cmp ] < 16 ) ? ( 1 << 4 ) : ctx.segm_cnt[ cmp ], 2 );
	auto mod_sgn = init_model_b( 9, 1 );

	const int bc = ctx.cmpnfo[ cmp ].bc;
	const int w  = ctx.cmpnfo[ cmp ].bch;

	std::vector<unsigned short> absv_store_vec( static_cast<std::size_t>( bc ), 0 );
	std::vector<unsigned char>  sgn_store_vec( static_cast<std::size_t>( bc ), 0 );
	std::vector<unsigned char>  zdstls_vec = ctx.zdstdata[ cmp ];
	unsigned short* absv_store = absv_store_vec.data();
	unsigned char*  sgn_store  = sgn_store_vec.data();
	unsigned char*  zdstls     = zdstls_vec.data();

	unsigned char* eob_x = ctx.eobxhigh[ cmp ].data();
	unsigned char* eob_y = ctx.eobyhigh[ cmp ].data();
	std::memset( eob_x, 0x00, static_cast<std::size_t>( bc ) );
	std::memset( eob_y, 0x00, static_cast<std::size_t>( bc ) );

	int c_absc_off[ 6 ];
	int c_weight[ 6 ];

	for ( int i = 1; i < 64; i++ ) {
		const int bpos = static_cast<int>( ctx.freqscan[ cmp ][ i ] );
		const int b_x  = unzigzag[ bpos ] % 8;
		const int b_y  = unzigzag[ bpos ] / 8;

		if ( ( b_x == 0 ) || ( b_y == 0 ) )
			continue;

		std::memset( absv_store, 0x00, static_cast<std::size_t>( bc ) * sizeof( unsigned short ) );
		std::memset( sgn_store, 0x00, static_cast<std::size_t>( bc ) );

		pjg_aavrg_prepare( ctx, c_absc_off, c_weight, cmp );

		const short* coeffs = ctx.colldata[ cmp ][ bpos ].data();
		const int max_val = ctx.max_v( cmp, bpos );
		const int max_len = bitlen_1024p( max_val );

		for ( int row_base = 0, p_y = 0; row_base < bc; row_base += w, ++p_y ) {
			const int row_w = std::min( w, bc - row_base );
			const unsigned char* row_zdst = zdstls + row_base;
			for ( int p_x = skip_zero_zdst( row_zdst, 0, row_w ); p_x < row_w; p_x = skip_zero_zdst( row_zdst, p_x + 1, row_w ) ) {
				const int dpos = row_base + p_x;
				const int r_x  = w - ( p_x + 1 );

				const int snum    = segm_tab[ zdstls[ dpos ] ];
				const int ctx_avr = pjg_aavrg_context( absv_store, c_absc_off, c_weight, dpos, p_y, p_x, r_x );
				const int ctx_len = bitlen_1024p( ctx_avr );
				shift_model( mod_len, ctx_len, snum );
				mod_len->exclude_symbols( max_len );

				if ( coeffs[ dpos ] == 0 ) {
					encode_ari( enc, mod_len, 0 );
				}
				else {
					const int absv = abs_val( coeffs[ dpos ] );
					const int clen = bitlen_1024p( absv );
					const int sgn  = ( coeffs[ dpos ] > 0 ) ? 0 : 1;
					encode_ari( enc, mod_len, clen );
					for ( int bp = clen - 2; bp >= 0; bp-- ) {
						shift_model( mod_res, snum, bp );
						const int bt = bit_n( absv, bp );
						encode_ari( enc, mod_res, bt );
					}
					int ctx_sgn = ( p_x > 0 ) ? sgn_store[ dpos - 1 ] : 0;
					if ( p_y > 0 ) ctx_sgn += 3 * sgn_store[ dpos - w ];
					mod_sgn->shift_context( ctx_sgn );
					encode_ari( enc, mod_sgn, sgn );

					absv_store[ dpos ] = static_cast<unsigned short>( absv );
					sgn_store[ dpos ]  = static_cast<unsigned char>( sgn + 1 );
					zdstls[ dpos ]--;
					if ( b_x > eob_x[ dpos ] ) eob_x[ dpos ] = static_cast<unsigned char>( b_x );
					if ( b_y > eob_y[ dpos ] ) eob_y[ dpos ] = static_cast<unsigned char>( b_y );
				}
			}
		}

		mod_len->flush_model();
		mod_res->flush_model();
		mod_sgn->flush_model();
	}

	return true;
}

bool pjg_encode_ac_low( CodecContext& ctx, ArithmeticEncoder* enc, int cmp )
{
	auto mod_len = init_model_s( 11, ( ctx.segm_cnt[ cmp ] > 11 ) ? ctx.segm_cnt[ cmp ] : 11, 2 );
	auto mod_res = init_model_b( 1 << 4, 2 );
	auto mod_top = init_model_b( ( ctx.nois_trs[ cmp ] > 4 ) ? ( 1 << ctx.nois_trs[ cmp ] ) : ( 1 << 4 ), 3, 16 );
	auto mod_sgn = init_model_b( 11, 1 );

	const int bc = ctx.cmpnfo[ cmp ].bc;
	const int w  = ctx.cmpnfo[ cmp ].bch;

	for ( int i = 2; i < 16; i++ ) {
		const int b_x  = ( i % 2 == 0 ) ? ( i / 2 ) : 0;
		const int b_y  = ( i % 2 == 1 ) ? ( i / 2 ) : 0;
		const int bpos = static_cast<int>( zigzag[ b_x + ( 8 * b_y ) ] );

		const short* coeffs = ctx.colldata[ cmp ][ bpos ].data();
		const LakhBandSpec lspec = setup_lakh_band( ctx, cmp, w, b_x, b_y );
		unsigned char* zdstls = lspec.zdstls;

		const int max_valp = ctx.max_v( cmp, bpos );
		const int max_valn = -max_valp;
		const int max_len  = bitlen_1024p( max_valp );
		const int thrs_bp  = ( max_len > ctx.nois_trs[ cmp ] ) ? ( max_len - ctx.nois_trs[ cmp ] ) : 0;

		for ( int row_base = 0, p_y = 0; row_base < bc; row_base += w, ++p_y ) {
			const int row_w = std::min( w, bc - row_base );
			const unsigned char* row_zdst = zdstls + row_base;
			for ( int p_x = skip_zero_zdst( row_zdst, 0, row_w ); p_x < row_w; p_x = skip_zero_zdst( row_zdst, p_x + 1, row_w ) ) {
				const int dpos = row_base + p_x;
				int ctx_lak = ( ( lspec.check_x ? p_x : p_y ) > 0 )
					? pjg_lakh_context( lspec, dpos )
					: 0;
				ctx_lak = clamped( max_valn, max_valp, ctx_lak );
				const int ctx_len = bitlen_2048n( ctx_lak );

				shift_model( mod_len, ctx_len, zdstls[ dpos ] );
				mod_len->exclude_symbols( max_len );

				if ( coeffs[ dpos ] == 0 ) {
					encode_ari( enc, mod_len, 0 );
				}
				else {
					const int absv = abs_val( coeffs[ dpos ] );
					const int clen = bitlen_2048n( absv );
					const int sgn  = ( coeffs[ dpos ] > 0 ) ? 0 : 1;
					encode_ari( enc, mod_len, clen );
					int bp = clen - 2;
					int ctx_res = ( bp >= thrs_bp ) ? 1 : 0;
					const int ctx_abs = abs_val( ctx_lak );
					const int ctx_sgn = ( ctx_lak == 0 ) ? 0 : ( ( ctx_lak > 0 ) ? 1 : 2 );
					for ( ; bp >= thrs_bp; bp-- ) {
						mod_top->select_context_3( ctx_abs >> thrs_bp, ctx_res, clen - thrs_bp );
						const int bt = bit_n( absv, bp );
						encode_ari( enc, mod_top, bt );
						ctx_res = ( ctx_res << 1 ) | ( bt ? 1 : 0 );
					}
					for ( ; bp >= 0; bp-- ) {
						shift_model( mod_res, zdstls[ dpos ], bp );
						const int bt = bit_n( absv, bp );
						encode_ari( enc, mod_res, bt );
					}
					mod_sgn->shift_context( ctx_sgn );
					encode_ari( enc, mod_sgn, sgn );
					zdstls[ dpos ]--;
				}
			}
		}

		mod_len->flush_model();
		mod_res->flush_model();
		mod_top->flush_model();
		mod_sgn->flush_model();
	}

	return true;
}

bool pjg_encode_generic( ArithmeticEncoder* enc, const unsigned char* data, std::size_t len )
{
	auto model = init_model_s( 256 + 1, 256, 1 );
	for ( std::size_t i = 0; i < len; i++ ) {
		encode_ari( enc, model, data[ i ] );
		model->shift_context( data[ i ] );
	}
	encode_ari( enc, model, 256 );
	return true;
}

bool pjg_encode_bit( ArithmeticEncoder* enc, unsigned char bit )
{
	auto model = init_model_b( 1, -1 );
	encode_ari( enc, model, bit );
	return true;
}

bool pjg_decode_zstscan( CodecContext& ctx, ArithmeticDecoder* dec, int cmp )
{
	unsigned char freqlist[ 64 ];

	ctx.zsrtscan[ cmp ][ 0 ] = 0;
	for ( int i = 0; i < 64; i++ )
		freqlist[ i ] = stdscan[ i ];

	auto model = init_model_s( 64, 64, 1 );

	for ( int i = 1; i < 64; i++ ) {
		model->exclude_symbols( 64 - i );

		int cpos = decode_ari( dec, model );
		model->shift_context( cpos );

		if ( cpos == 0 ) {
			for ( int tpos = 0; i < 64; i++ ) {
				do {
					++tpos;
				} while ( ( tpos < 64 ) && ( freqlist[ tpos ] == 0 ) );
				if ( tpos >= 64 ) {
					throw std::runtime_error( "corrupted arithmetic stream" );
				}
				ctx.zsrtscan[ cmp ][ i ] = freqlist[ tpos ];
			}
			break;
		}

		int tpos = 0;
		for ( ; tpos < 64; tpos++ ) {
			if ( freqlist[ tpos ] != 0 ) cpos--;
			if ( cpos == 0 ) break;
		}
		if ( tpos >= 64 ) {
			throw std::runtime_error( "corrupted arithmetic stream" );
		}

		ctx.zsrtscan[ cmp ][ i ] = freqlist[ tpos ];
		freqlist[ tpos ] = 0;
	}

	ctx.freqscan[ cmp ] = ctx.zsrtscan[ cmp ];
	return true;
}

bool pjg_decode_zdst_high( CodecContext& ctx, ArithmeticDecoder* dec, int cmp )
{
	auto model = init_model_s( 49 + 1, 25 + 1, 1 );
	unsigned char* zdstls = ctx.zdstdata[ cmp ].data();
	const int w  = ctx.cmpnfo[ cmp ].bch;
	const int bc = ctx.cmpnfo[ cmp ].bc;

	for ( int dpos = 0; dpos < bc; dpos++ ) {
		int a = 0;
		int b = 0;
		get_context_nnb( dpos, w, &a, &b );
		a = ( a >= 0 ) ? zdstls[ a ] : 0;
		b = ( b >= 0 ) ? zdstls[ b ] : 0;
		model->shift_context( ( a + b + 2 ) / 4 );
		zdstls[ dpos ] = static_cast<unsigned char>( decode_ari( dec, model ) );
	}

	return true;
}

bool pjg_decode_zdst_low( CodecContext& ctx, ArithmeticDecoder* dec, int cmp )
{
	auto model = init_model_s( 8, 8, 2 );
	unsigned char* zdstls_x = ctx.zdstxlow[ cmp ].data();
	unsigned char* zdstls_y = ctx.zdstylow[ cmp ].data();
	const unsigned char* ctx_eobx = ctx.eobxhigh[ cmp ].data();
	const unsigned char* ctx_eoby = ctx.eobyhigh[ cmp ].data();
	const unsigned char* ctx_zdst = ctx.zdstdata[ cmp ].data();
	const int bc = ctx.cmpnfo[ cmp ].bc;

	for ( int dpos = 0; dpos < bc; dpos++ ) {
		model->shift_context( ( ctx_zdst[ dpos ] + 3 ) / 7 );
		model->shift_context( ctx_eobx[ dpos ] );
		zdstls_x[ dpos ] = static_cast<unsigned char>( decode_ari( dec, model ) );
	}
	for ( int dpos = 0; dpos < bc; dpos++ ) {
		model->shift_context( ( ctx_zdst[ dpos ] + 3 ) / 7 );
		model->shift_context( ctx_eoby[ dpos ] );
		zdstls_y[ dpos ] = static_cast<unsigned char>( decode_ari( dec, model ) );
	}

	return true;
}

bool pjg_decode_dc( CodecContext& ctx, ArithmeticDecoder* dec, int cmp )
{
	const unsigned char* segm_tab = segm_tables[ ctx.segm_cnt[ cmp ] - 1 ];

	const int max_val = ctx.max_v( cmp, 0 );
	const int max_len = bitlen_1024p( max_val );

	auto mod_len = init_model_s( max_len + 1, ( ctx.segm_cnt[ cmp ] > max_len ) ? ctx.segm_cnt[ cmp ] : max_len + 1, 2 );
	auto mod_res = init_model_b( ( ctx.segm_cnt[ cmp ] < 16 ) ? ( 1 << 4 ) : ctx.segm_cnt[ cmp ], 2 );
	auto mod_sgn = init_model_b( 1, 0 );

	const int bc = ctx.cmpnfo[ cmp ].bc;
	const int w  = ctx.cmpnfo[ cmp ].bch;

	std::vector<unsigned short> absv_store_vec( static_cast<std::size_t>( bc ), 0 );
	unsigned short* absv_store = absv_store_vec.data();

	int c_absc_off[ 6 ];
	int c_weight[ 6 ];
	pjg_aavrg_prepare( ctx, c_absc_off, c_weight, cmp );

	short* coeffs = ctx.colldata[ cmp ][ 0 ].data();
	const unsigned char* zdstls = ctx.zdstdata[ cmp ].data();

	for ( int dpos = 0, p_y = 0; dpos < bc; ++p_y ) {
		const int row_w = std::min( w, bc - dpos );
		for ( int p_x = 0; p_x < row_w; ++p_x, ++dpos ) {
			const int r_x = w - ( p_x + 1 );

			const int snum    = segm_tab[ zdstls[ dpos ] ];
			const int ctx_avr = pjg_aavrg_context( absv_store, c_absc_off, c_weight, dpos, p_y, p_x, r_x );
			const int ctx_len = bitlen_1024p( ctx_avr );
			shift_model( mod_len, ctx_len, snum );
			const int clen = decode_ari( dec, mod_len );

			if ( clen != 0 ) {
				int absv = 1;
				for ( int bp = clen - 2; bp >= 0; bp-- ) {
					shift_model( mod_res, snum, bp );
					const int bt = decode_ari( dec, mod_res );
					absv = ( absv << 1 ) | ( bt ? 1 : 0 );
				}
				const int sgn = decode_ari( dec, mod_sgn );
				coeffs[ dpos ]     = static_cast<short>( ( sgn == 0 ) ? absv : -absv );
				absv_store[ dpos ] = static_cast<unsigned short>( absv );
			}
		}
	}

	return true;
}

bool pjg_decode_ac_high( CodecContext& ctx, ArithmeticDecoder* dec, int cmp )
{
	const unsigned char* segm_tab = segm_tables[ ctx.segm_cnt[ cmp ] - 1 ];

	auto mod_len = init_model_s( 11, ( ctx.segm_cnt[ cmp ] > 11 ) ? ctx.segm_cnt[ cmp ] : 11, 2 );
	auto mod_res = init_model_b( ( ctx.segm_cnt[ cmp ] < 16 ) ? ( 1 << 4 ) : ctx.segm_cnt[ cmp ], 2 );
	auto mod_sgn = init_model_b( 9, 1 );

	const int bc = ctx.cmpnfo[ cmp ].bc;
	const int w  = ctx.cmpnfo[ cmp ].bch;

	std::vector<unsigned short> absv_store_vec( static_cast<std::size_t>( bc ), 0 );
	std::vector<unsigned char>  sgn_store_vec( static_cast<std::size_t>( bc ), 0 );
	std::vector<unsigned char>  zdstls_vec = ctx.zdstdata[ cmp ];
	unsigned short* absv_store = absv_store_vec.data();
	unsigned char*  sgn_store  = sgn_store_vec.data();
	unsigned char*  zdstls     = zdstls_vec.data();

	unsigned char* eob_x = ctx.eobxhigh[ cmp ].data();
	unsigned char* eob_y = ctx.eobyhigh[ cmp ].data();
	std::memset( eob_x, 0x00, static_cast<std::size_t>( bc ) );
	std::memset( eob_y, 0x00, static_cast<std::size_t>( bc ) );

	int c_absc_off[ 6 ];
	int c_weight[ 6 ];

	for ( int i = 1; i < 64; i++ ) {
		const int bpos = static_cast<int>( ctx.freqscan[ cmp ][ i ] );
		const int b_x  = unzigzag[ bpos ] % 8;
		const int b_y  = unzigzag[ bpos ] / 8;

		if ( ( b_x == 0 ) || ( b_y == 0 ) )
			continue;

		std::memset( absv_store, 0x00, static_cast<std::size_t>( bc ) * sizeof( unsigned short ) );
		std::memset( sgn_store, 0x00, static_cast<std::size_t>( bc ) );

		pjg_aavrg_prepare( ctx, c_absc_off, c_weight, cmp );

		short* coeffs = ctx.colldata[ cmp ][ bpos ].data();
		const int max_val = ctx.max_v( cmp, bpos );
		const int max_len = bitlen_1024p( max_val );

		for ( int row_base = 0, p_y = 0; row_base < bc; row_base += w, ++p_y ) {
			const int row_w = std::min( w, bc - row_base );
			const unsigned char* row_zdst = zdstls + row_base;
			for ( int p_x = skip_zero_zdst( row_zdst, 0, row_w ); p_x < row_w; p_x = skip_zero_zdst( row_zdst, p_x + 1, row_w ) ) {
				const int dpos = row_base + p_x;
				const int r_x  = w - ( p_x + 1 );

				const int snum    = segm_tab[ zdstls[ dpos ] ];
				const int ctx_avr = pjg_aavrg_context( absv_store, c_absc_off, c_weight, dpos, p_y, p_x, r_x );
				const int ctx_len = bitlen_1024p( ctx_avr );
				shift_model( mod_len, ctx_len, snum );
				mod_len->exclude_symbols( max_len );

				const int clen = decode_ari( dec, mod_len );
				if ( clen != 0 ) {
					int absv = 1;
					for ( int bp = clen - 2; bp >= 0; bp-- ) {
						shift_model( mod_res, snum, bp );
						const int bt = decode_ari( dec, mod_res );
						absv = ( absv << 1 ) | ( bt ? 1 : 0 );
					}
					int ctx_sgn = ( p_x > 0 ) ? sgn_store[ dpos - 1 ] : 0;
					if ( p_y > 0 ) ctx_sgn += 3 * sgn_store[ dpos - w ];
					mod_sgn->shift_context( ctx_sgn );
					const int sgn = decode_ari( dec, mod_sgn );

					coeffs[ dpos ]     = static_cast<short>( ( sgn == 0 ) ? absv : -absv );
					absv_store[ dpos ] = static_cast<unsigned short>( absv );
					sgn_store[ dpos ]  = static_cast<unsigned char>( sgn + 1 );
					zdstls[ dpos ]--;
					if ( b_x > eob_x[ dpos ] ) eob_x[ dpos ] = static_cast<unsigned char>( b_x );
					if ( b_y > eob_y[ dpos ] ) eob_y[ dpos ] = static_cast<unsigned char>( b_y );
				}
			}
		}

		mod_len->flush_model();
		mod_res->flush_model();
		mod_sgn->flush_model();
	}

	return true;
}

bool pjg_decode_ac_low( CodecContext& ctx, ArithmeticDecoder* dec, int cmp )
{
	auto mod_len = init_model_s( 11, ( ctx.segm_cnt[ cmp ] > 11 ) ? ctx.segm_cnt[ cmp ] : 11, 2 );
	auto mod_res = init_model_b( 1 << 4, 2 );
	auto mod_top = init_model_b( ( ctx.nois_trs[ cmp ] > 4 ) ? ( 1 << ctx.nois_trs[ cmp ] ) : ( 1 << 4 ), 3, 16 );
	auto mod_sgn = init_model_b( 11, 1 );

	const int bc = ctx.cmpnfo[ cmp ].bc;
	const int w  = ctx.cmpnfo[ cmp ].bch;

	for ( int i = 2; i < 16; i++ ) {
		const int b_x  = ( i % 2 == 0 ) ? ( i / 2 ) : 0;
		const int b_y  = ( i % 2 == 1 ) ? ( i / 2 ) : 0;
		const int bpos = static_cast<int>( zigzag[ b_x + ( 8 * b_y ) ] );

		short* coeffs = ctx.colldata[ cmp ][ bpos ].data();
		const LakhBandSpec lspec = setup_lakh_band( ctx, cmp, w, b_x, b_y );
		unsigned char* zdstls = lspec.zdstls;

		const int max_valp = ctx.max_v( cmp, bpos );
		const int max_valn = -max_valp;
		const int max_len  = bitlen_1024p( max_valp );
		const int thrs_bp  = ( max_len > ctx.nois_trs[ cmp ] ) ? ( max_len - ctx.nois_trs[ cmp ] ) : 0;

		for ( int row_base = 0, p_y = 0; row_base < bc; row_base += w, ++p_y ) {
			const int row_w = std::min( w, bc - row_base );
			const unsigned char* row_zdst = zdstls + row_base;
			for ( int p_x = skip_zero_zdst( row_zdst, 0, row_w ); p_x < row_w; p_x = skip_zero_zdst( row_zdst, p_x + 1, row_w ) ) {
				const int dpos = row_base + p_x;
				int ctx_lak = ( ( lspec.check_x ? p_x : p_y ) > 0 )
					? pjg_lakh_context( lspec, dpos )
					: 0;
				ctx_lak = clamped( max_valn, max_valp, ctx_lak );
				const int ctx_len = bitlen_2048n( ctx_lak );

				shift_model( mod_len, ctx_len, zdstls[ dpos ] );
				mod_len->exclude_symbols( max_len );

				const int clen = decode_ari( dec, mod_len );
				if ( clen != 0 ) {
					int bp = clen - 2;
					int ctx_res = ( bp >= thrs_bp ) ? 1 : 0;
					const int ctx_abs = abs_val( ctx_lak );
					const int ctx_sgn = ( ctx_lak == 0 ) ? 0 : ( ( ctx_lak > 0 ) ? 1 : 2 );
					for ( ; bp >= thrs_bp; bp-- ) {
						mod_top->select_context_3( ctx_abs >> thrs_bp, ctx_res, clen - thrs_bp );
						const int bt = decode_ari( dec, mod_top );
						ctx_res = ( ctx_res << 1 ) | ( bt ? 1 : 0 );
					}
					int absv = ( ctx_res == 0 ) ? 1 : ctx_res;
					for ( ; bp >= 0; bp-- ) {
						shift_model( mod_res, zdstls[ dpos ], bp );
						const int bt = decode_ari( dec, mod_res );
						absv = ( absv << 1 ) | ( bt ? 1 : 0 );
					}
					mod_sgn->shift_context( ctx_sgn );
					const int sgn = decode_ari( dec, mod_sgn );
					coeffs[ dpos ] = static_cast<short>( ( sgn == 0 ) ? absv : -absv );
					zdstls[ dpos ]--;
				}
			}
		}

		mod_len->flush_model();
		mod_res->flush_model();
		mod_top->flush_model();
		mod_sgn->flush_model();
	}

	return true;
}

bool pjg_decode_generic( ArithmeticDecoder* dec, std::vector<unsigned char>& out, std::size_t max_bytes )
{
	MemoryWriter bwrt;

	auto model = init_model_s( 256 + 1, 256, 1 );
	while ( true ) {
		const int c = decode_ari( dec, model );
		if ( c == 256 ) break;
		if ( bwrt.num_bytes_written() >= max_bytes ) {
			return false;
		}
		bwrt.write_byte( static_cast<unsigned char>( c ) );
		model->shift_context( c );
	}

	out = bwrt.take_data();
	return true;
}

bool pjg_decode_bit( ArithmeticDecoder* dec, unsigned char* bit )
{
	auto model = init_model_b( 1, -1 );
	*bit = static_cast<unsigned char>( decode_ari( dec, model ) );
	return true;
}

#if defined(USE_PLOCOI)
inline int plocoi( int a, int b, int c )
{
	const int min_v = ( a < b ) ? a : b;
	const int max_v = ( a > b ) ? a : b;
	if ( c >= max_v ) return min_v;
	if ( c <= min_v ) return max_v;
	return a + b - c;
}

int dc_coll_predictor( const CodecContext& ctx, int cmp, int dpos )
{
	const short* coefs = ctx.colldata[ cmp ][ 0 ].data();
	const int w = ctx.cmpnfo[ cmp ].bch;
	int a = 0;
	int b = 0;
	int c = 0;

	if ( dpos < w ) {
		a = coefs[ dpos - 1 ];
	}
	else if ( ( dpos % w ) == 0 ) {
		b = coefs[ dpos - w ];
	}
	else {
		a = coefs[ dpos - 1 ];
		b = coefs[ dpos - w ];
		c = coefs[ dpos - 1 - w ];
	}

	return plocoi( a, b, c );
}
#else
struct Dc1dTables {
	const short* b8x1[ 8 ]{};
	const short* b1x8[ 8 ]{};
	const int* w8x1_0 = nullptr;
	const int* w8x1_7 = nullptr;
	const int* w1x8_0 = nullptr;
	const int* w1x8_7 = nullptr;
};

Dc1dTables make_dc1d_tables( const CodecContext& ctx, int cmp ) noexcept
{
	static constexpr int kB8x1[ 8 ] = { 0, 1, 5, 6, 14, 15, 27, 28 };
	static constexpr int kB1x8[ 8 ] = { 0, 2, 3, 9, 10, 20, 21, 35 };
	Dc1dTables t{};
	for ( int k = 0; k < 8; k++ ) {
		t.b8x1[ k ] = ctx.colldata[ cmp ][ kB8x1[ k ] ].data();
		t.b1x8[ k ] = ctx.colldata[ cmp ][ kB1x8[ k ] ].data();
	}
	t.w8x1_0 = &ctx.adpt_idct_8x1[ cmp ][ 0 ];
	t.w8x1_7 = &ctx.adpt_idct_8x1[ cmp ][ 56 ];
	t.w1x8_0 = &ctx.adpt_idct_1x8[ cmp ][ 0 ];
	t.w1x8_7 = &ctx.adpt_idct_1x8[ cmp ][ 56 ];
	return t;
}

int eval_dc1d_top_only( const Dc1dTables& t, int rs, int w ) noexcept
{
	int pb = t.b1x8[ 0 ][ rs - w ] * t.w1x8_7[ 0 ];
	int xb = 0;
	for ( int k = 1; k < 8; k++ ) {
		pb += t.b1x8[ k ][ rs - w ] * t.w1x8_7[ k ];
		xb += t.b1x8[ k ][ rs ]     * t.w1x8_0[ k ];
	}
	return ( pb - xb ) * 8;
}

void accum_dc1d_row_base(
	const Dc1dTables& t,
	int rs,
	int w,
	int n,
	bool has_top,
	int* PJG_RESTRICT rbuf ) noexcept
{
	if ( has_top ) {
		const short* PJG_RESTRICT p_dc_t = t.b1x8[ 0 ] + rs - w + 1;
		const int wt0 = t.w1x8_7[ 0 ];
		for ( int i = 0; i < n; i++ ) {
			rbuf[ i ] = p_dc_t[ i ] * wt0;
		}
		for ( int k = 1; k < 8; k++ ) {
			const short* PJG_RESTRICT bh_l = t.b8x1[ k ] + rs;
			const short* PJG_RESTRICT bh_c = t.b8x1[ k ] + rs + 1;
			const short* PJG_RESTRICT bv_t = t.b1x8[ k ] + rs - w + 1;
			const short* PJG_RESTRICT bv_c = t.b1x8[ k ] + rs + 1;
			const int wh7 = t.w8x1_7[ k ], wh0 = t.w8x1_0[ k ];
			const int wv7 = t.w1x8_7[ k ], wv0 = t.w1x8_0[ k ];
			for ( int i = 0; i < n; i++ ) {
				rbuf[ i ] += ( bh_l[ i ] * wh7 - bh_c[ i ] * wh0 ) + ( bv_t[ i ] * wv7 - bv_c[ i ] * wv0 );
			}
		}
	}
	else {
		std::fill( rbuf, rbuf + n, 0 );
		for ( int k = 1; k < 8; k++ ) {
			const short* PJG_RESTRICT bh_l = t.b8x1[ k ] + rs;
			const short* PJG_RESTRICT bh_c = t.b8x1[ k ] + rs + 1;
			const int wh7 = t.w8x1_7[ k ], wh0 = t.w8x1_0[ k ];
			for ( int i = 0; i < n; i++ ) {
				rbuf[ i ] += bh_l[ i ] * wh7 - bh_c[ i ] * wh0;
			}
		}
	}
}

inline short finish_dc_diff( int raw_pred, short dc_val, int q0, int absmaxp, int corr_f, int sign ) noexcept
{
	int pred = clamped( -( 1024 * DCT_RSC_FACTOR ), ( 1016 * DCT_RSC_FACTOR ), raw_pred );
	pred = ( q0 > 0 ) ? dct_rescale( pred / q0 ) : 0;
	int val = static_cast<int>( dc_val ) + ( sign * pred );
	if ( val > absmaxp ) val -= corr_f;
	else if ( val < -absmaxp ) val += corr_f;
	return static_cast<short>( val );
}
#endif

} // namespace


/* -----------------------------------------------
	adapt ICOS tables for quantizer tables
	----------------------------------------------- */
bool adapt_icos( CodecContext& ctx )
{
	unsigned short quant[ 64 ];
	if ( ( ctx.action == A_PGM_DUMP ) && !ctx.adpt_idct_8x8 ) {
		ctx.adpt_idct_8x8 = std::make_unique<std::array<std::array<int, 4096>, 4>>();
	}

	for ( int cmp = 0; cmp < ctx.cmpc; cmp++ ) {
		for ( int ipos = 0; ipos < 64; ipos++ ) {
			quant[ ipos ] = ctx.quant( cmp, zigzag[ ipos ] );
			if ( quant[ ipos ] >= 2048 )
				quant[ ipos ] = 0;
		}
		if ( ctx.action == A_PGM_DUMP ) {
			for ( std::size_t ipos = 0; ipos < 4096u; ipos++ )
				( *ctx.adpt_idct_8x8 )[ cmp ][ ipos ] = icos_idct_8x8[ ipos ] * quant[ ipos % 64u ];
		}
		for ( std::size_t ipos = 0; ipos < 64u; ipos++ )
			ctx.adpt_idct_1x8[ cmp ][ ipos ] = icos_idct_1x8[ ipos ] * quant[ ( ipos % 8u ) * 8u ];
		for ( std::size_t ipos = 0; ipos < 64u; ipos++ )
			ctx.adpt_idct_8x1[ cmp ][ ipos ] = icos_idct_1x8[ ipos ] * quant[ ipos % 8u ];
	}

	return true;
}


/* -----------------------------------------------
	filter DC coefficients
	----------------------------------------------- */
bool predict_dc( CodecContext& ctx )
{
#if defined(USE_PLOCOI)
	for ( int cmp = 0; cmp < ctx.cmpc; cmp++ ) {
		const int absmaxp = ctx.max_v( cmp, 0 );
		const int absmaxn = -absmaxp;
		const int corr_f  = ( 2 * absmaxp ) + 1;

		for ( int dpos = ctx.cmpnfo[ cmp ].bc - 1; dpos > 0; dpos-- ) {
			short* coef = &ctx.colldata[ cmp ][ 0 ][ dpos ];
			int val = static_cast<int>( *coef ) - dc_coll_predictor( ctx, cmp, dpos );
			if ( val > absmaxp ) val -= corr_f;
			else if ( val < absmaxn ) val += corr_f;
			*coef = static_cast<short>( val );
		}
	}
#else
	std::vector<int> raw_row;
	for ( int cmp = 0; cmp < ctx.cmpc; cmp++ ) {
		const int w  = ctx.cmpnfo[ cmp ].bch;
		const int bc = ctx.cmpnfo[ cmp ].bc;
		if ( ( bc <= 1 ) || ( w <= 0 ) ) continue;

		const int absmaxp   = ctx.max_v( cmp, 0 );
		const int corr_f    = ( 2 * absmaxp ) + 1;
		const int q0        = ctx.quant( cmp, 0 );
		const int rows      = ( bc + w - 1 ) / w;
		const Dc1dTables t  = make_dc1d_tables( ctx, cmp );
		const int wl0       = t.w8x1_7[ 0 ];
		short* dc           = ctx.colldata[ cmp ][ 0 ].data();

		if ( w > 1 ) {
			raw_row.resize( static_cast<std::size_t>( w - 1 ) );
		}
		for ( int py = rows - 1; py >= 0; --py ) {
			const int rs = py * w;
			const int rw = std::min( w, bc - rs );
			if ( rw > 1 ) {
				const int n     = rw - 1;
				const int scale = ( py > 0 ) ? ( 8 / 2 ) : 8;
				accum_dc1d_row_base( t, rs, w, n, py > 0, raw_row.data() );
				for ( int i = n - 1; i >= 0; --i ) {
					const int raw = ( raw_row[ static_cast<std::size_t>( i ) ] + static_cast<int>( dc[ rs + i ] ) * wl0 ) * scale;
					dc[ rs + 1 + i ] = finish_dc_diff( raw, dc[ rs + 1 + i ], q0, absmaxp, corr_f, -1 );
				}
			}
			if ( ( py > 0 ) && ( rw > 0 ) ) {
				dc[ rs ] = finish_dc_diff( eval_dc1d_top_only( t, rs, w ), dc[ rs ], q0, absmaxp, corr_f, -1 );
			}
		}
	}
#endif

	return true;
}


/* -----------------------------------------------
	unpredict DC coefficients
	----------------------------------------------- */
bool unpredict_dc( CodecContext& ctx )
{
#if defined(USE_PLOCOI)
	for ( int cmp = 0; cmp < ctx.cmpc; cmp++ ) {
		const int absmaxp = ctx.max_v( cmp, 0 );
		const int absmaxn = -absmaxp;
		const int corr_f  = ( 2 * absmaxp ) + 1;

		for ( int dpos = 1; dpos < ctx.cmpnfo[ cmp ].bc; dpos++ ) {
			short* coef = &ctx.colldata[ cmp ][ 0 ][ dpos ];
			int val = static_cast<int>( *coef ) + dc_coll_predictor( ctx, cmp, dpos );
			if ( val > absmaxp ) val -= corr_f;
			else if ( val < absmaxn ) val += corr_f;
			*coef = static_cast<short>( val );
		}
	}
#else
	std::vector<int> raw_row;
	for ( int cmp = 0; cmp < ctx.cmpc; cmp++ ) {
		const int w  = ctx.cmpnfo[ cmp ].bch;
		const int bc = ctx.cmpnfo[ cmp ].bc;
		if ( ( bc <= 1 ) || ( w <= 0 ) ) continue;

		const int absmaxp   = ctx.max_v( cmp, 0 );
		const int corr_f    = ( 2 * absmaxp ) + 1;
		const int q0        = ctx.quant( cmp, 0 );
		const int rows      = ( bc + w - 1 ) / w;
		const Dc1dTables t  = make_dc1d_tables( ctx, cmp );
		const int wl0       = t.w8x1_7[ 0 ];
		short* dc           = ctx.colldata[ cmp ][ 0 ].data();

		if ( w > 1 ) {
			raw_row.resize( static_cast<std::size_t>( w - 1 ) );
		}
		for ( int py = 0; py < rows; ++py ) {
			const int rs = py * w;
			const int rw = std::min( w, bc - rs );
			if ( ( py > 0 ) && ( rw > 0 ) ) {
				dc[ rs ] = finish_dc_diff( eval_dc1d_top_only( t, rs, w ), dc[ rs ], q0, absmaxp, corr_f, +1 );
			}
			if ( rw > 1 ) {
				const int n     = rw - 1;
				const int scale = ( py > 0 ) ? ( 8 / 2 ) : 8;
				accum_dc1d_row_base( t, rs, w, n, py > 0, raw_row.data() );
				for ( int i = 0; i < n; ++i ) {
					const int raw = ( raw_row[ static_cast<std::size_t>( i ) ] + static_cast<int>( dc[ rs + i ] ) * wl0 ) * scale;
					dc[ rs + 1 + i ] = finish_dc_diff( raw, dc[ rs + 1 + i ], q0, absmaxp, corr_f, +1 );
				}
			}
		}
	}
#endif

	return true;
}


/* -----------------------------------------------
	checks range of values, error if out of bounds
	----------------------------------------------- */
bool check_value_range( CodecContext& ctx )
{
	for ( int cmp = 0; cmp < ctx.cmpc; cmp++ ) {
		const int bc = ctx.cmpnfo[ cmp ].bc;
		for ( int bpos = 0; bpos < 64; bpos++ ) {
			const int absmax = ctx.max_v( cmp, bpos );
			const short hi = static_cast<short>( std::min( absmax, 32767 ) );
			const short lo = static_cast<short>( -hi );
			const short* PJG_RESTRICT band = ctx.colldata[ cmp ][ bpos ].data();
			unsigned char out_of_range = 0;
			for ( int dpos = 0; dpos < bc; dpos++ ) {
				const short v = band[ dpos ];
				out_of_range |= static_cast<unsigned char>( ( v > hi ) | ( v < lo ) );
			}
			if ( out_of_range != 0 ) {
				for ( int dpos = 0; dpos < bc; dpos++ ) {
					if ( ( band[ dpos ] > absmax ) || ( band[ dpos ] < -absmax ) ) {
						ctx.set_errormessage( "value out of range error: cmp%i, frq%i, val %i, max %i",
							cmp, bpos, static_cast<int>( band[ dpos ] ), absmax );
						ctx.errorlevel = 2;
						return false;
					}
				}
			}
		}
	}

	return true;
}


/* -----------------------------------------------
	calculate zero distribution lists
	----------------------------------------------- */
bool calc_zdst_lists( CodecContext& ctx )
{
	for ( int cmp = 0; cmp < ctx.cmpc; cmp++ ) {
		const int bc = ctx.cmpnfo[ cmp ].bc;
		std::fill( ctx.zdstdata[ cmp ].begin(), ctx.zdstdata[ cmp ].end(), 0 );
		std::fill( ctx.zdstxlow[ cmp ].begin(), ctx.zdstxlow[ cmp ].end(), 0 );
		std::fill( ctx.zdstylow[ cmp ].begin(), ctx.zdstylow[ cmp ].end(), 0 );

		for ( int bpos = 1; bpos < 64; bpos++ ) {
			const int b_x = unzigzag[ bpos ] % 8;
			const int b_y = unzigzag[ bpos ] / 8;
			const short* PJG_RESTRICT band = ctx.colldata[ cmp ][ bpos ].data();
			unsigned char* PJG_RESTRICT dst =
				( b_x == 0 ) ? ctx.zdstylow[ cmp ].data() :
				( b_y == 0 ) ? ctx.zdstxlow[ cmp ].data() :
				               ctx.zdstdata[ cmp ].data();
			for ( int dpos = 0; dpos < bc; dpos++ ) {
				dst[ dpos ] = static_cast<unsigned char>( dst[ dpos ] + ( band[ dpos ] != 0 ) );
			}
		}
	}

	return true;
}


/* -----------------------------------------------
	packs all parts to compressed pjg
	----------------------------------------------- */
bool pack_pjg( CodecContext& ctx )
{
	if ( ctx.disc_meta ) {
		if ( !jpg_rebuild_header( ctx ) ) return false;
	}

	if ( ( ctx.hdrdata.size() > kMaxGenericStreamBytes ) ||
	     ( ctx.grbgdata.size() > kMaxGenericStreamBytes ) ||
	     ( static_cast<std::size_t>( std::max( 0, ctx.scnc ) ) > kMaxRstErrBytes ) ) {
		ctx.errormessage = "header or trailing data exceeds maximum supported size";
		ctx.errorlevel   = 2;
		return false;
	}

	std::vector<unsigned char> saved_hdr = ctx.hdrdata;
	if ( !pjg_optimize_header( ctx ) ) {
		ctx.hdrdata = std::move( saved_hdr );
		return false;
	}
	std::vector<unsigned char> opt_hdr = std::move( ctx.hdrdata );
	ctx.hdrdata = std::move( saved_hdr );

	if ( ctx.padbit == -1 ) ctx.padbit = 1;

	ctx.str_out->write( kPjgMagic, 2 );

	if ( !ctx.auto_set ) {
		ctx.str_out->write_byte( 0x00 );
		ctx.str_out->write( ctx.nois_trs.data(), 4 );
		ctx.str_out->write( ctx.segm_cnt.data(), 4 );
	}

	ctx.str_out->write_byte( kAppVersion );

	ArithmeticEncoder encoder_obj( *ctx.str_out );
	ArithmeticEncoder* encoder = &encoder_obj;
	auto fail = [&]() {
		encoder_obj.abandon();
		return false;
	};

	if ( !pjg_encode_generic( encoder, opt_hdr.data(), opt_hdr.size() ) ) return fail();
	if ( !pjg_encode_bit( encoder, static_cast<unsigned char>( ctx.padbit ) ) ) return fail();
	if ( !pjg_encode_bit( encoder, ctx.has_rst_err ? 1 : 0 ) ) return fail();
	if ( ctx.has_rst_err ) {
		if ( ctx.rst_err.size() < static_cast<std::size_t>( ctx.scnc ) )
			ctx.rst_err.resize( static_cast<std::size_t>( ctx.scnc ), 0 );
		if ( !pjg_encode_generic( encoder, ctx.rst_err.data(), static_cast<std::size_t>( ctx.scnc ) ) )
			return fail();
	}

	for ( int cmp = 0; cmp < ctx.cmpc; cmp++ ) {
		if ( !pjg_encode_zstscan( ctx, encoder, cmp ) ) return fail();
		if ( !pjg_encode_zdst_high( ctx, encoder, cmp ) ) return fail();
		if ( !pjg_encode_ac_high( ctx, encoder, cmp ) ) return fail();
		if ( !pjg_encode_zdst_low( ctx, encoder, cmp ) ) return fail();
		if ( !pjg_encode_ac_low( ctx, encoder, cmp ) ) return fail();
		if ( !pjg_encode_dc( ctx, encoder, cmp ) ) return fail();
	}

	if ( !pjg_encode_bit( encoder, ctx.grbgdata.empty() ? 0 : 1 ) ) return fail();
	if ( !ctx.grbgdata.empty() ) {
		if ( !pjg_encode_generic( encoder, ctx.grbgdata.data(), ctx.grbgdata.size() ) )
			return fail();
	}

	encoder_obj.finalize();

	if ( ctx.str_out->error() ) {
		ctx.errormessage = "write error, possibly drive is full";
		ctx.errorlevel = 2;
		return false;
	}

	ctx.pjgfilesize = static_cast<std::int64_t>( ctx.str_out->num_bytes_written() );
	return true;
}


/* -----------------------------------------------
	unpacks compressed pjg to colldata
	----------------------------------------------- */
bool unpack_pjg( CodecContext& ctx )
{
	unsigned char hcode = 0;
	unsigned char cb    = 0;

	ctx.auto_set = true;

	while ( true ) {
		if ( !ctx.str_in->read_byte( &hcode ) ) {
			ctx.errormessage = "unexpected end of data encountered";
			ctx.errorlevel = 2;
			return false;
		}
		if ( hcode == 0x00 ) {
			if ( ( ctx.str_in->read( ctx.nois_trs.data(), 4 ) != 4 ) ||
				 ( ctx.str_in->read( ctx.segm_cnt.data(), 4 ) != 4 ) ) {
				ctx.errormessage = "unexpected end of data encountered";
				ctx.errorlevel = 2;
				return false;
			}
			for ( int cmp = 0; cmp < 4; cmp++ ) {
				if ( ( ctx.nois_trs[ cmp ] > 10 ) || ( ctx.segm_cnt[ cmp ] < 1 ) || ( ctx.segm_cnt[ cmp ] > 49 ) ) {
					ctx.errormessage = "invalid compression settings in header";
					ctx.errorlevel = 2;
					return false;
				}
			}
			ctx.auto_set = false;
		}
		else if ( hcode >= 0x14 ) {
			if ( hcode != kAppVersion ) {
				ctx.set_errormessage( "incompatible file, use %s v%i.%i",
					kAppName, hcode / 10, hcode % 10 );
				ctx.errorlevel = 2;
				return false;
			}
			break;
		}
		else {
			ctx.set_errormessage( "unknown header code, use newer version of %s", kAppName );
			ctx.errorlevel = 2;
			return false;
		}
	}

	const std::size_t in_bytes = ( ctx.str_in != nullptr ) ? ctx.str_in->get_size() : 0;
	const std::size_t max_generic = ( ctx.str_in == nullptr )
		? kMaxGenericStreamBytes
		: std::min<std::size_t>(
			kMaxGenericStreamBytes,
			std::max<std::size_t>(
				4096u,
				( in_bytes <= ( kMaxGenericStreamBytes / kMaxGenericExpansionPerByte ) )
					? ( in_bytes * kMaxGenericExpansionPerByte )
					: kMaxGenericStreamBytes ) );

	ArithmeticDecoder decoder_obj( *ctx.str_in );
	ArithmeticDecoder* decoder = &decoder_obj;

	if ( !pjg_decode_generic( decoder, ctx.hdrdata, max_generic ) ) {
		ctx.errormessage = "corrupted or oversized header data";
		ctx.errorlevel = 2;
		return false;
	}
	if ( !pjg_decode_bit( decoder, &cb ) ) return false;
	ctx.padbit = cb;

	if ( !pjg_decode_bit( decoder, &cb ) ) return false;
	ctx.has_rst_err = ( cb == 1 );
	if ( ctx.has_rst_err ) {
		if ( !pjg_decode_generic( decoder, ctx.rst_err, kMaxRstErrBytes ) ) {
			ctx.errormessage = "corrupted restart error data";
			ctx.errorlevel = 2;
			return false;
		}
	}
	else {
		ctx.rst_err.clear();
	}

	if ( !pjg_unoptimize_header( ctx ) ) return false;
	if ( ctx.disc_meta ) {
		if ( !jpg_rebuild_header( ctx ) ) return false;
	}
	if ( !jpg_setup_imginfo( ctx, true ) ) return false;

	for ( int cmp = 0; cmp < ctx.cmpc; cmp++ ) {
		if ( !pjg_decode_zstscan( ctx, decoder, cmp ) ) return false;
		if ( !pjg_decode_zdst_high( ctx, decoder, cmp ) ) return false;
		if ( !pjg_decode_ac_high( ctx, decoder, cmp ) ) return false;
		if ( !pjg_decode_zdst_low( ctx, decoder, cmp ) ) return false;
		if ( !pjg_decode_ac_low( ctx, decoder, cmp ) ) return false;
		if ( !pjg_decode_dc( ctx, decoder, cmp ) ) return false;
	}

	if ( !pjg_decode_bit( decoder, &cb ) ) return false;
	if ( cb == 0 ) {
		ctx.grbgdata.clear();
	}
	else if ( !pjg_decode_generic( decoder, ctx.grbgdata, max_generic ) ) {
		ctx.errormessage = "corrupted or oversized trailing data";
		ctx.errorlevel = 2;
		return false;
	}

	ctx.pjgfilesize = static_cast<std::int64_t>( ( ctx.str_in != nullptr ) ? ctx.str_in->get_size() : 0 );
	return true;
}


/* -----------------------------------------------
	inverse DCT transform using precalc tables (fast)
	----------------------------------------------- */
int idct_2d_fst_8x8( const CodecContext& ctx, int cmp, int dpos, int ix, int iy )
{
	if ( !ctx.adpt_idct_8x8 ) {
		return 0;
	}
	std::int64_t idct = 0;
	const int ixy = ( ( iy << 3 ) + ix ) << 6;
	const auto& tbl = ( *ctx.adpt_idct_8x8 )[ cmp ];

	for ( int i = 0; i < 64; i++ ) {
		idct += static_cast<std::int64_t>( ctx.colldata[ cmp ][ zigzag[ i ] ][ dpos ] ) * tbl[ ixy + i ];
	}

	if ( idct < std::numeric_limits<int>::min() ) return std::numeric_limits<int>::min();
	if ( idct > std::numeric_limits<int>::max() ) return std::numeric_limits<int>::max();
	return static_cast<int>( idct );
}

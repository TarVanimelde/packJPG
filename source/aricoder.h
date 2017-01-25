#ifndef aricoder_h
#define aricoder_h

#include <cstdint>

#include "bitops.h"
#include <vector>

// defines for coder
constexpr uint32_t CODER_USE_BITS = 31; // Must never be above 31.
constexpr uint32_t CODER_LIMIT100 = uint32_t(1 << CODER_USE_BITS);
constexpr uint32_t CODER_LIMIT025 = CODER_LIMIT100 / 4;
constexpr uint32_t CODER_LIMIT050 = (CODER_LIMIT100 / 4) * 2;
constexpr uint32_t CODER_LIMIT075 = (CODER_LIMIT100 / 4) * 3;
constexpr uint32_t CODER_MAXSCALE = CODER_LIMIT025 - 1;
constexpr uint32_t ESCAPE_SYMBOL = CODER_LIMIT025;

// symbol struct, used in arithmetic coding
struct symbol {
    unsigned int low_count;
	unsigned int high_count;
	unsigned int scale;
};

// table struct, used in in statistical models,
// holding all info needed for one context
struct table {
	// counts for each symbol contained in the table
	std::vector<uint16_t> counts;
	// links to higher order contexts
	std::vector<table*> links;
	// accumulated counts
	unsigned int scale = unsigned int(0);
};

// special table struct, used in in model_s,
// holding additional info for a speedier 'totalize_table'
struct table_s {
	// counts for each symbol contained in the table
	std::vector<uint16_t> counts;
	// links to higher order contexts
	std::vector<table_s*> links;
	// speedup info
	unsigned short max_count = unsigned short(0);
	unsigned short max_symbol = unsigned short(0);
};


/* -----------------------------------------------
	class for arithmetic coding of data to/from iostream
	----------------------------------------------- */
	
class aricoder
{
	public:
	aricoder( iostream* stream, int iomode );
	~aricoder();
	void encode( symbol* s );
	unsigned int decode_count( symbol* s );
	void decode( symbol* s );
	
	private:

	template<uint8_t bit>
	void write_bit() {
		// add bit at last position
		bbyte = (bbyte << 1) | bit;
		// increment bit position
		cbit++;

		// write bit if done
		if (cbit == 8) {
			sptr->write(&bbyte, 1, 1);
			cbit = 0;
		}
	}
	
	void writeNrbitsAsZero();
	void writeNrbitsAsOne();
	unsigned char read_bit();
	
	// i/o variables
	iostream* sptr;
	int mode;
	unsigned char bbyte;
	unsigned char cbit;
	
	// arithmetic coding variables
	unsigned int ccode;
	unsigned int clow;
	unsigned int chigh;
	unsigned int cstep;
	unsigned int nrbits;
};


/* -----------------------------------------------
	universal statistical model for arithmetic coding
	----------------------------------------------- */
	
class model_s
{	
	public:
	
	model_s( int max_s, int max_c, int max_o, int c_lim );
	~model_s();
	
	void update_model( int symbol );
	void shift_context( int c );
	void flush_model( int scale_factor );
	void exclude_symbols( char rule, int c );
	
	int  convert_int_to_symbol( int c, symbol *s );
	void get_symbol_scale( symbol *s );
	int  convert_symbol_to_int( int count, symbol *s );
	
	private:
	
	std::vector<uint32_t> totals;
	bool* scoreboard;
	int sb0_count;
	std::vector<table_s*> contexts;
	
	int max_symbol;
	int max_context;
	int current_order;
	int max_order;
	int max_count;
	
	inline void totalize_table(table_s* context );
	inline void rescale_table(table_s* context, int scale_factor );
	inline void recursive_flush(table_s* context, int scale_factor );
	inline void recursive_cleanup(table_s* context );
};


/* -----------------------------------------------
	binary statistical model for arithmetic coding
	----------------------------------------------- */
	
class model_b
{	
	public:
	
	model_b( int max_c, int max_o, int c_lim );
	~model_b();
	
	void update_model( int symbol );
	void shift_context( int c );
	void flush_model( int scale_factor );
	
	int  convert_int_to_symbol( int c, symbol *s );
	void get_symbol_scale( symbol *s );
	int  convert_symbol_to_int( int count, symbol *s );	
	
	private:
	
	std::vector<table*> contexts;
	
	int max_context;
	int max_order;
	int max_count;
	
	inline void check_counts( table *context );
	inline void rescale_table( table* context, int scale_factor );
	inline void recursive_flush( table* context, int scale_factor );
	inline void recursive_cleanup( table *context );
};

// Base case for shifting an arbitrary number of contexts into the model.
template <typename M>
static void shift_model(M model) {}

// Shift an arbitrary number of contexts into the model (at most max_c contexts).
template <typename M, typename C, typename... Cargs>
static void shift_model(M model, C context, Cargs ... contextList) {
	model->shift_context(context);
	shift_model(model, contextList...);
}

/* -----------------------------------------------
	generic model_s encoder function
	----------------------------------------------- */
static inline void encode_ari( aricoder* encoder, model_s* model, int c )
{
	symbol s;
	int esc;
	
	do {		
		esc = model->convert_int_to_symbol( c, &s );
		encoder->encode( &s );
	} while ( esc );
	model->update_model( c );
}

/* -----------------------------------------------
	generic model_s decoder function
	----------------------------------------------- */	
static inline int decode_ari( aricoder* decoder, model_s* model )
{
	symbol s;
	uint32_t count;
	int c;
	
	do{
		model->get_symbol_scale( &s );
		count = decoder->decode_count( &s );
		c = model->convert_symbol_to_int( count, &s );
		decoder->decode( &s );	
	} while ( c == ESCAPE_SYMBOL );
	model->update_model( c );
	
	return c;
}

/* -----------------------------------------------
	generic model_b encoder function
	----------------------------------------------- */	
static inline void encode_ari( aricoder* encoder, model_b* model, int c )
{
	symbol s;
	
	model->convert_int_to_symbol( c, &s );
	encoder->encode( &s );
	model->update_model( c );
}

/* -----------------------------------------------
	generic model_b decoder function
	----------------------------------------------- */	
static inline int decode_ari( aricoder* decoder, model_b* model )
{
	symbol s;
	
	model->get_symbol_scale( &s );
	uint32_t count = decoder->decode_count( &s );
	int c = model->convert_symbol_to_int( count, &s );
	decoder->decode( &s );	
	model->update_model( c );
	
	return c;
}

#endif
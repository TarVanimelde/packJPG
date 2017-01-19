// defines for coder
#define CODER_USE_BITS		31 // must never be above 31
#define CODER_LIMIT100		( (unsigned int) ( 1 << CODER_USE_BITS ) )
#define CODER_LIMIT025		( ( CODER_LIMIT100 / 4 ) * 1 )
#define CODER_LIMIT050		( ( CODER_LIMIT100 / 4 ) * 2 )
#define CODER_LIMIT075		( ( CODER_LIMIT100 / 4 ) * 3 )
#define CODER_MAXSCALE		CODER_LIMIT025 - 1
#define ESCAPE_SYMBOL		CODER_LIMIT025


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
	unsigned short* counts = nullptr;
	// links to higher order contexts
	struct table** links = nullptr;
	// link to lower order context
	struct table* lesser = nullptr;
	// accumulated counts
	unsigned int scale = 0;
};

// special table struct, used in in model_s,
// holding additional info for a speedier 'totalize_table'
struct table_s {
	// counts for each symbol contained in the table
	unsigned short* counts = nullptr;
	// links to higher order contexts
	struct table_s** links = nullptr;
	// link to lower order context
	struct table_s* lesser = nullptr;
	// speedup info
	unsigned short max_count = unsigned short(0);
	unsigned short max_symbol = unsigned short(0);
	// unsigned short esc_prob;
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
	
	bool error;
	
	
	private:
	
	// unsigned short* totals;
	unsigned int* totals;
	bool* scoreboard;
	int sb0_count;
	table_s **contexts;
	table_s **storage;
	
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
	
	bool error;
	
	
	private:
	
	table **contexts;
	table **storage;
	
	int max_context;
	int max_order;
	int max_count;
	
	inline void check_counts( table *context );
	inline void rescale_table( table* context, int scale_factor );
	inline void recursive_flush( table* context, int scale_factor );
	inline void recursive_cleanup( table *context );
};


/* -----------------------------------------------
	shift context x2 model_s function
	----------------------------------------------- */
static inline void shift_model( model_s* model, int ctx1, int ctx2 )
{
	model->shift_context( ctx1 );
	model->shift_context( ctx2 );
}


/* -----------------------------------------------
	shift context x3 model_s function
	----------------------------------------------- */
static inline void shift_model( model_s* model, int ctx1, int ctx2, int ctx3 )
{
	model->shift_context( ctx1 );
	model->shift_context( ctx2 );
	model->shift_context( ctx3 );
}


/* -----------------------------------------------
	shift context x2 model_b function
	----------------------------------------------- */
static inline void shift_model( model_b* model, int ctx1, int ctx2 )
{
	model->shift_context( ctx1 );
	model->shift_context( ctx2 );
}


/* -----------------------------------------------
	shift context x3 model_b function
	----------------------------------------------- */
static inline void shift_model( model_b* model, int ctx1, int ctx2, int ctx3 )
{
	model->shift_context( ctx1 );
	model->shift_context( ctx2 );
	model->shift_context( ctx3 );
}


/* -----------------------------------------------
	generic model_s encoder function
	----------------------------------------------- */
static inline void encode_ari( aricoder* encoder, model_s* model, int c )
{
	static symbol s;
	static int esc;
	
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
	static symbol s;
	static unsigned int count;
	static int c;
	
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
	static symbol s;
	
	model->convert_int_to_symbol( c, &s );
	encoder->encode( &s );
	model->update_model( c );
}

/* -----------------------------------------------
	generic model_b decoder function
	----------------------------------------------- */	
static inline int decode_ari( aricoder* decoder, model_b* model )
{
	static symbol s;
	static unsigned int count;
	static int c;
	
	model->get_symbol_scale( &s );
	count = decoder->decode_count( &s );
	c = model->convert_symbol_to_int( count, &s );
	decoder->decode( &s );	
	model->update_model( c );
	
	return c;
}

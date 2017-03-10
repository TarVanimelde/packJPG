#ifndef BITOPS_H
#define BITOPS_H

#define RBITS( c, n )		( c & ( 0xFF >> (8 - n) ) )
#define LBITS( c, n )		( c >> (8 - n) )
#define MBITS( c, l, r )	( RBITS( c,l ) >> r )
#define RBITS32( c, n )		( c & ( 0xFFFFFFFF >> (32 - n) ) )
#define MBITS32( c, l, r )	( RBITS32( c,l ) >> r )
#define BITN( c, n )		( (c >> n) & 0x1 )

#include <cstdint>
#include <memory>
#include <vector>

enum StreamType {
	kFile = 0,
	kMemory = 1,
	kStream = 2
};

enum StreamMode {
	kRead = 0,
	kWrite = 1
};

	
/* -----------------------------------------------
	class to read arrays bitwise
	----------------------------------------------- */

class abitreader {
public:
	abitreader(const std::vector<std::uint8_t>& bits);
	~abitreader();
	unsigned int read(int nbits);
	unsigned char read_bit();
	unsigned char unpad(unsigned char fillbit);
	int getpos();
	int getbitp();
	void setpos(int pbyte, int pbit);
	void rewind_bits(int nbits);
	bool eof();
	int peof();

private:
	std::vector<std::uint8_t> data;
	int cbyte = 0;
	int cbit = 8;
	int peof_ = 0;
	bool eof_ = false;
};


/* -----------------------------------------------
	class to write arrays bitwise
	----------------------------------------------- */

class abitwriter {
public:
	abitwriter(int size);
	~abitwriter();
	void write(unsigned int val, int nbits);
	void write_bit(unsigned char bit);
	void set_fillbit(unsigned char fillbit);
	void pad();
	std::vector<std::uint8_t> get_data();
	int getpos();
	int getbitp();

private:
	unsigned char fillbit_ = 1;
	std::vector<std::uint8_t> data;
	int cbyte = 0;
	int cbit = 8;
};


/* -----------------------------------------------
	class to read arrays bytewise
	----------------------------------------------- */

class abytereader {
public:
	abytereader(const std::vector<std::uint8_t>& bytes);
	~abytereader();
	int read(unsigned char* byte);
	int read_n(unsigned char* byte, int n);
	void seek(int pos);
	int getsize();
	int getpos();
	bool eof();

private:
	const std::vector<std::uint8_t> data;
	int cbyte = 0;
	bool _eof = false;
};


/* -----------------------------------------------
	class to write arrays bytewise
	----------------------------------------------- */

class abytewriter
{
public:
	abytewriter( int size );
	~abytewriter();	
	void write( unsigned char byte );
	void write_n(const unsigned char* byte, int n );
	std::vector<std::uint8_t> get_data();
	int getpos();
	void reset();
	
private:
	std::vector<std::uint8_t> data;
	int cbyte = 0;
};


/* -----------------------------------------------
	class for input and output from file or memory
	----------------------------------------------- */

class iostream
{
public:
	// Memory iostream
	iostream(const std::vector<std::uint8_t>& bytes, StreamMode iomode);
	// File iostream
	iostream(const std::string& file_path, StreamMode iomode);
	// stdout/in iostream.
	iostream(StreamMode iomode);
	~iostream();
	void switch_mode();
	int read(unsigned char* to, int dtsize);
	bool read_byte(unsigned char* to);
	int write(const unsigned char* from, int dtsize);
	int write_byte(unsigned char byte);
	int flush();
	int rewind();
	int getpos();
	int getsize();
	std::vector<std::uint8_t> get_data();
	bool chkerr();
	bool chkeof();
	
private:
	void open_file();
	void open_mem();
	void open_stream();
	
	int write_file(const unsigned char* from, int dtsize);
	int write_file_byte(unsigned char byte);
	int read_file(unsigned char* to, int dtsize);
	bool read_file_byte(unsigned char* to);
	int write_mem(const unsigned char* from, int dtsize );
	int write_mem_byte(unsigned char byte);
	int read_mem(unsigned char* to, int dtsize);
	bool read_mem_byte(unsigned char* to);
	
	FILE* fptr;
	std::vector<char> file_buffer; // Used to replace the default file buffer for reads/writes to improve performance.
  
	std::unique_ptr<abytewriter> mwrt;
	std::unique_ptr<abytereader> mrdr;
	
	const std::string file_path;
	std::vector<std::uint8_t> data;
	StreamMode mode;
	StreamType srct;
};

#endif
package uk.ac.ebi.ena.sra.cram.io;

import java.io.IOException;

/**
 * Discard data and count number of bits flying through. 
 * 
 * @author vadim
 * 
 */
public class CountingNullBitOutputStream implements BitOutputStream {
	private long size = 0L;
	
	public void reset () {
		size = 0L ;
	}
	
	public long getSize () {
		return size ;
	}

	@Override
	public final void write(int b, int nbits) throws IOException {
		size += nbits ;
	}

	@Override
	public final void write(long b, int nbits) throws IOException {
		size += nbits ;
	}

	@Override
	public final void write(byte b, int nbits) throws IOException {
		size += nbits ;
	}

	@Override
	public final void write(boolean bit) throws IOException {
		size ++ ;
	}

	@Override
	public final void write(boolean bit, long repeat) throws IOException {
		size += repeat ;
	}

	@Override
	public final void flush() throws IOException {
	}

}

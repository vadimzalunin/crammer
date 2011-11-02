package uk.ac.ebi.ena.sra.cram.io;

import java.io.IOException;

/** Does nothing. 
 * @author vadim
 *
 */
public class NullBitOutputStream implements BitOutputStream {
	
	public final static NullBitOutputStream INSTANCE = new NullBitOutputStream();

	@Override
	public final void write(int b, int nbits) throws IOException {
	}

	@Override
	public final void write(long b, int nbits) throws IOException {
	}

	@Override
	public final void write(byte b, int nbits) throws IOException {
	}

	@Override
	public final void write(boolean bit) throws IOException {
	}

	@Override
	public final void write(boolean bit, long repeat) throws IOException {
	}

	@Override
	public final void flush() throws IOException {
	}

	@Override
	public void close() throws IOException {
	}

}

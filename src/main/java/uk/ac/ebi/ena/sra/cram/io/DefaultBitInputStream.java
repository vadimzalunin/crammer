package uk.ac.ebi.ena.sra.cram.io;

import java.io.DataInputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;

/**
 * Must not read from delegate unless no bits left in the buffer!!!
 * 
 * @author vadim
 * 
 */
public class DefaultBitInputStream extends DataInputStream implements
		BitInputStream {

	private int nofBufferedBits = 0;
	private int byteBuffer = 0;
	private boolean endOfStream = false;

	public DefaultBitInputStream(InputStream in) {
		super(in);
	}

	public final boolean readBit() throws IOException {
		if (--nofBufferedBits >= 0)
			return ((byteBuffer >>> nofBufferedBits) & 1) == 1;

		nofBufferedBits = 7;
		byteBuffer = in.read();
		if (byteBuffer == -1) {
			endOfStream = true;
			throw new EOFException("End of stream.");
		}

		return ((byteBuffer >>> 7) & 1) == 1;
	}

	public final int readBits(int n) throws IOException {
		if (n == 0) return 0 ;
		int x = 0;
		while (n > nofBufferedBits) {
			n -= nofBufferedBits;
			x |= rightBits(nofBufferedBits, byteBuffer) << n;
			byteBuffer = in.read();
			if (byteBuffer == -1) {
				endOfStream = true;
				throw new EOFException("End of stream.");
			}

			nofBufferedBits = 8;
		}
		nofBufferedBits -= n;
		return x | rightBits(n, byteBuffer >>> nofBufferedBits);
	}

	private static final int rightBits(int n, int x) {
		return x & ((1 << n) - 1);
	}

	public final long readLongBits(int len) throws IOException {
		if (len > 64)
			throw new RuntimeException(
					"More then 64 bits are requested in one read from bit stream.");

		long result = 0;
		final long last = len - 1;
		for (long bi = 0; bi <= last; bi++) {
			final boolean frag = readBit();
			if (frag)
				result |= 1L << (last - bi);
		}
		return result;
	}

	public void reset() {
		nofBufferedBits = 0;
		byteBuffer = 0;
	}

	@Override
	public boolean endOfStream() throws IOException {
		return endOfStream;
	}

	public int getNofBufferedBits() {
		return nofBufferedBits;
	}
}
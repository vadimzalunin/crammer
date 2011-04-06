package uk.ac.ebi.ena.sra.cram.io;

import java.io.DataInputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;

public class DefaultBitInputStream extends DataInputStream implements
		BitInputStream {

	private int leftBits = 0;
	private int byteBuffer = 0;

	public DefaultBitInputStream(InputStream in) {
		super(in);
	}

	public final boolean readBit() throws IOException {
		if (--leftBits >= 0)
			return ((byteBuffer >>> leftBits) & 1) == 1;

		leftBits = 7;
		byteBuffer = in.read();
		if (byteBuffer == -1)
			throw new EOFException("End of stream.");

		return ((byteBuffer >>> 7) & 1) == 1;
	}

	public final int readBits(int n) throws IOException {
		int x = 0;
		while (n > leftBits) {
			n -= leftBits;
			x |= rightBits(leftBits, byteBuffer) << n;
			byteBuffer = in.read();
			if (byteBuffer == -1)
				throw new EOFException("End of stream.");

			leftBits = 8;
		}
		leftBits -= n;
		return x | rightBits(n, byteBuffer >>> leftBits);
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
}
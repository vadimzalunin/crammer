package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.IOException;

import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;

public class GolombCodec implements BitCodec<Long> {
	private long m;
	private int log_m;
	private boolean quotientBit = false ;

	public GolombCodec(int log_m) {
		this.m = 1 << log_m;
		this.log_m = log_m;
	}

	public final Long read(final BitInputStream bis) throws IOException {

		long unary = 0;
		while (bis.readBit() == quotientBit)
			unary++;

		long remainder = readBits(bis, log_m);

		long result = unary * m + remainder;
		return result;
	}

	private static final long readBits(final BitInputStream bis, final long len)
			throws IOException {
		if (len > 64)
			throw new RuntimeException(
					"More then 64 bits are requested in one read from bit stream.");

		long result = 0;
		final long last = len - 1;
		for (long bi = 0; bi <= last; bi++) {
			final boolean frag = bis.readBit();
			if (frag)
				result |= 1L << (last - bi);
		}
		return result;
	}

	@Override
	public long write(final BitOutputStream bos, final Long value)
			throws IOException {
		long quotient = value / m;
		if (quotient > 0x7fffffffL)
			for (long i = 0; i < quotient; i++)
//				bos.write(false);
		bos.write(quotientBit);

		else if (quotient > 0) {
			final int qi = (int) quotient;
			for (int i = 0; i < qi; i++)
//				bos.write(false);
			bos.write(quotientBit);
		}
//		bos.write((byte) 1, 1);
		bos.write(!quotientBit);
		long remainder = value % m;
		long mask = 1 << (log_m - 1);
		for (int i = log_m - 1; i >= 0; i--) {
			final long b = remainder & mask;
			bos.write(b != 0L);
			mask >>>= 1;
		}
		long bits = quotient + 1 + log_m;
		return bits;
	}

	@Override
	public long numberOfBits(Long value) {
		return value / m + 1 + log_m;
	}

}
package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.IOException;

import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;

public class GolombRiceCodec implements BitCodec<Long> {
	private long m;
	private int log_m;
	private boolean quotientBit = false ;

	public GolombRiceCodec(int log_m) {
		this.m = 1 << log_m;
		this.log_m = log_m;
	}

	public final Long read(final BitInputStream bis) throws IOException {

		long unary = 0;
		while (bis.readBit() == quotientBit)
			unary++;

		long remainder = bis.readBits(log_m);

		long result = unary * m + remainder;
		return result;
	}

	@Override
	public long write(final BitOutputStream bos, final Long value)
			throws IOException {
		long quotient = value / m;
		if (quotient > 0x7fffffffL)
			for (long i = 0; i < quotient; i++)
		bos.write(quotientBit);

		else if (quotient > 0) {
			final int qi = (int) quotient;
			for (int i = 0; i < qi; i++)
			bos.write(quotientBit);
		}
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

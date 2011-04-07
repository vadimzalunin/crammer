package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.IOException;

import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;

public class GolombRiceCodec implements BitCodec<Long> {
	private long m;
	private int log2m;
	private boolean quotientBit = false;
	private long offset = 0L;

	public GolombRiceCodec(int log2m, boolean quotientBit, long offset) {
		super();
		this.log2m = log2m;
		m = 1<< log2m ;
		this.quotientBit = quotientBit;
		this.offset = offset;
	}

	public GolombRiceCodec(int log2m) {
		this(log2m, false, 0L) ;
	}

	public final Long read(final BitInputStream bis) throws IOException {

		long unary = 0;
		while (bis.readBit() == quotientBit)
			unary++;

		long remainder = bis.readBits(log2m);

		long result = unary * m + remainder;
		return result - offset;
	}

	@Override
	public final long write(final BitOutputStream bos, final Long value)
			throws IOException {
		long newValue = value + offset;
		long quotient = newValue / m;
		if (quotient > 0x7fffffffL)
			for (long i = 0; i < quotient; i++)
				bos.write(quotientBit);

		else if (quotient > 0) {
			final int qi = (int) quotient;
			for (int i = 0; i < qi; i++)
				bos.write(quotientBit);
		}
		bos.write(!quotientBit);
		long remainder = newValue % m;
		long mask = 1 << (log2m - 1);
		for (int i = log2m - 1; i >= 0; i--) {
			final long b = remainder & mask;
			bos.write(b != 0L);
			mask >>>= 1;
		}
		long bits = quotient + 1 + log2m;
		return bits;
	}

	@Override
	public final long numberOfBits(Long value) {
		return (value+offset) / m + 1 + log2m;
	}

	public int getLog2m() {
		return log2m;
	}

	public void setLog2m(int log2m) {
		this.log2m = log2m;
		m = 1<< log2m ;
	}

	public boolean isQuotientBit() {
		return quotientBit;
	}

	public void setQuotientBit(boolean quotientBit) {
		this.quotientBit = quotientBit;
	}

	public long getOffset() {
		return offset;
	}

	public void setOffset(long offset) {
		this.offset = offset;
	}

}

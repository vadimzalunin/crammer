package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.IOException;

import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;

public class GolombCodec implements BitCodec<Long> {
	private long m;
	private boolean quotientBit = true;
	private Long offset = 0L;

	public GolombCodec(int m) {
		this(m, true, 0L);
	}

	public GolombCodec(long m, boolean quotientBit, Long offset) {
		this.m = m;
		this.quotientBit = quotientBit;
		this.offset = offset;
	}

	public final Long read(final BitInputStream bis) throws IOException {
		long quotient = 0L;
		while (bis.readBit() == quotientBit)
			quotient++;

		long numbits = quotient + 1;
		long ceiling = (long) Math.ceil(Math.log(m) / Math.log(2));
		long reminder = bis.readBits((int) (ceiling - 1));
		numbits += ceiling - 1;
		if (reminder >= Math.pow(2, ceiling) - m) {
			reminder <<= 1;
			reminder |= bis.readBits(1);
			reminder -= Math.pow(2, ceiling) - m;
		}

		return (quotient * m + reminder) - offset;
	}

	@Override
	public final long write(final BitOutputStream bos, final Long value)
			throws IOException {
		long newValue = value + offset;
		long quotient = (long) Math.floor(newValue / m);
		long reminder = newValue % m;
		long ceiling = (long) Math.ceil(Math.log(m) / Math.log(2));

		long len = quotient + 1;
		bos.write(quotientBit, quotient);
		bos.write(!quotientBit);

		if (reminder < Math.pow(2, ceiling) - m) {
			bos.write(reminder, (int) ceiling - 1);
			len += ceiling - 1;
		} else {
			bos.write((int) (reminder + Math.pow(2, ceiling) - m),
					(int) ceiling);
			len += ceiling;
		}
		return len;
	}

	@Override
	public final long numberOfBits(Long value) {
		long newValue = value + offset;
		long quotient = (long) Math.floor(newValue / m);
		long reminder = newValue % m;
		long ceiling = (long) Math.ceil(Math.log(m) / Math.log(2));
		long l = quotient + 1;

		if (reminder < Math.pow(2, ceiling) - m)
			l += ceiling - 1;
		else
			l += ceiling;

		return l;
	}

	public long getM() {
		return m;
	}

	public boolean isQuotientBit() {
		return quotientBit;
	}

	public Long getOffset() {
		return offset;
	}

	public void setM(long m) {
		this.m = m;
	}

	public void setQuotientBit(boolean quotientBit) {
		this.quotientBit = quotientBit;
	}

	public void setOffset(Long offset) {
		this.offset = offset;
	}

}

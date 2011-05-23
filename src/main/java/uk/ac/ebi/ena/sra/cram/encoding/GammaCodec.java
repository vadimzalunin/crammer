package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.IOException;

import org.apache.commons.math.util.MathUtils;

import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;

public class GammaCodec implements BitCodec<Long> {
	private long offset = 0L;
	private boolean lenCodingBit = false;

	public GammaCodec(long offset, boolean lenCodingBit) {
		super();
		this.offset = offset;
		this.lenCodingBit = lenCodingBit;
	}

	public GammaCodec() {
		this(0L, false);
	}

	@Override
	public final Long read(BitInputStream bis) throws IOException {
		int len = 1;
		while (bis.readBit() == lenCodingBit)
			len++;
		int readBits = bis.readBits(len - 1);
		long value = readBits | 1 << (len - 1);
		return value - offset;
	}

	@Override
	public final long write(BitOutputStream bos, Long value) throws IOException {
		if (value + offset < 1)
			throw new IllegalArgumentException(
					"Gamma codec handles only positive values: " + value);

		long newValue = value + offset;
		int betaCodeLength = 1 + (int) MathUtils.log(2, newValue);
		if (betaCodeLength > 1)
			bos.write(0L, betaCodeLength - 1);

		bos.write(newValue, betaCodeLength);
		return betaCodeLength * 2 - 1;
	}

	@Override
	public final long numberOfBits(Long value) {
		long newValue = value + offset;
		int betaCodeLength = 1 + (int) MathUtils.log(2, newValue);
		return betaCodeLength * 2 - 1;
	}

	public long getOffset() {
		return offset;
	}

	public boolean isLenCodingBit() {
		return lenCodingBit;
	}

	public void setOffset(long offset) {
		this.offset = offset;
	}

	public void setLenCodingBit(boolean lenCodingBit) {
		this.lenCodingBit = lenCodingBit;
	}

}

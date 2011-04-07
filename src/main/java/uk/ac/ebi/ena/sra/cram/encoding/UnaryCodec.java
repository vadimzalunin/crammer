package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.IOException;

import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;

public class UnaryCodec implements BitCodec<Long> {
	private boolean stopBit = false;
	private long offset = 0L;

	public UnaryCodec() {
		this(false, 0L);
	}

	public UnaryCodec(boolean stopBit, long offset) {
		this.stopBit = stopBit;
		this.offset = offset;
	}

	@Override
	public final Long read(BitInputStream bis) throws IOException {
		long bits = 0;
		while (bis.readBit() != stopBit)
			bits++;

		return bits - offset - 1;
	}

	@Override
	public final long write(BitOutputStream bos, Long value) throws IOException {
		long bits = value + 1;

		while (bits-- > 0)
			bos.write(stopBit ? 0 : 1, 1);

		bos.write(stopBit ? 1 : 0, 1);

		return value + 1;
	}

	@Override
	public final long numberOfBits(Long value) {
		return value + 1;
	}

	public boolean isStopBit() {
		return stopBit;
	}

	public long getOffset() {
		return offset;
	}

	public void setStopBit(boolean stopBit) {
		this.stopBit = stopBit;
	}

	public void setOffset(long offset) {
		this.offset = offset;
	}

}

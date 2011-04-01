package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.IOException;

import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;

public class UnaryCodec implements BitCodec<Long> {
	private boolean stopBit = false;

	public UnaryCodec(boolean stopBit, long offset) {
		super();
		this.stopBit = stopBit;
		this.offset = offset;
	}

	private long offset = 0;

	@Override
	public Long read(BitInputStream bis) throws IOException {
		long bits = 0;
		while (bis.readBit() != stopBit)
			bits++;

		return bits - offset - 1;
	}

	@Override
	public long write(BitOutputStream bos, Long value) throws IOException {
		long bits = value + 1;

		while (bits-- > 0)
			bos.writeBits(stopBit ? 0 : 1, 1);

		bos.writeBits(stopBit ? 1 : 0, 1);

		return value + 1;
	}

}

package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.IOException;

import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;

public class BetaCodec implements BitCodec<Long> {
	private long offset = 0L;
	private int readNofBits;

	@Override
	public final Long read(BitInputStream bis) throws IOException {
		return bis.readLongBits(readNofBits) - offset;
	}

	@Override
	public final long write(BitOutputStream bos, Long value) throws IOException {
		if (value + offset < 0)
			throw new IllegalArgumentException("Value is less then offset: "
					+ value);

		int nofBits = (int) numberOfBits(value);
		long newValue = value + offset;
		bos.write(newValue, nofBits);
		return nofBits;
	}

	@Override
	public final long numberOfBits(Long value) {
		long newValue = value + offset;
		return (int) Math.floor(Math.log(newValue) / Math.log(2));
	}

	public long getOffset() {
		return offset;
	}

	public void setOffset(long offset) {
		this.offset = offset;
	}

	public int getReadNofBits() {
		return readNofBits;
	}

	public void setReadNofBits(int readNofBits) {
		this.readNofBits = readNofBits;
	}

}

package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.IOException;

import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;

public class LongDeltaCodec implements BitCodec<Long> {
	private BitCodec<Long> delegate;
	private long previousValue = 0L;

	public LongDeltaCodec(BitCodec<Long> delegate, long offset) {
		this.delegate = delegate;
		this.previousValue = offset;
	}

	@Override
	public Long read(BitInputStream bis) throws IOException {
		long value = previousValue + delegate.read(bis);
		previousValue = value;
		return value;
	}

	@Override
	public long write(BitOutputStream bos, Long object) throws IOException {
		long value = delegate.write(bos, object - previousValue);
		previousValue = value;
		return value;
	}

	@Override
	public long numberOfBits(Long object) {
		return delegate.numberOfBits(object - previousValue);
	}

	public void reset() {
		previousValue = 0L ;
	}

}

package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.IOException;

import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;

public class MeasuringCodec<T> implements BitCodec<T> {
	private BitCodec<T> delegate;
	private long readObjects = 0L;
	private long writtenObjects = 0L;
	private long writtenBits = 0L;
	private String name;

	public MeasuringCodec(BitCodec<T> delegate, String name) {
		this.delegate = delegate;
		this.name = name;
	}
	
	public BitCodec<T> getDelegate () {
		return delegate ;
	}

	@Override
	public T read(BitInputStream bis) throws IOException {
		T object = delegate.read(bis);
		readObjects++;
		return object;
	}

	@Override
	public long write(BitOutputStream bos, T object) throws IOException {
		long len = delegate.write(bos, object);
		writtenObjects++;
		writtenBits += len;
		return len;
	}

	@Override
	public long numberOfBits(T object) {
		return delegate.numberOfBits(object);
	}

	public void reset() {
		readObjects = 0L;
		writtenObjects = 0L;
		writtenBits = 0L;
	}

	public long getReadObjects() {
		return readObjects;
	}

	public long getWrittenObjects() {
		return writtenObjects;
	}

	public long getWrittenBits() {
		return writtenBits;
	}

	@Override
	public String toString() {
		return String
				.format("%s: written objects=%d; written bits=%d; written bits per object=%.2f",
						name, getWrittenObjects(),
						getWrittenBits(), (double) getWrittenBits()
								/ getWrittenObjects());
	}

	public String getName() {
		return name;
	}
}

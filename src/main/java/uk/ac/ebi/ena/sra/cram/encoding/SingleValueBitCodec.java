package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.IOException;

import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;

public class SingleValueBitCodec<T> implements BitCodec<T> {
	private T value;

	@Override
	public T read(BitInputStream bis) throws IOException {
		return value;
	}

	@Override
	public long write(BitOutputStream bos, T object) throws IOException {
		return 0;
	}

	@Override
	public long numberOfBits(T object) {
		return 0;
	}

	public T getValue() {
		return value;
	}

	public void setValue(T value) {
		this.value = value;
	}

}

package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.IOException;

import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;

public class NullBitCodec<T> implements BitCodec<T> {

	@Override
	public T read(BitInputStream bis) throws IOException {
		return null;
	}

	@Override
	public long write(BitOutputStream bos, T object) throws IOException {
		return 0;
	}

	@Override
	public long numberOfBits(T object) {
		return 0;
	}

}

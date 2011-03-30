package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.IOException;

import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;

public interface BitCodec<T> {

	public T read(BitInputStream bis) throws IOException;

	public long write(BitOutputStream bis, T object) throws IOException;
}

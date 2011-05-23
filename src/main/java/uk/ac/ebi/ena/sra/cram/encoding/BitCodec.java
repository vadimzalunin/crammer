package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.IOException;

import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;

public interface BitCodec<T> {

	public T read(BitInputStream bis) throws IOException;

	public long write(BitOutputStream bos, T object) throws IOException;

	public long numberOfBits(T object);

//	/**
//	 * Resets any state information, for example information about previously
//	 * seen values.
//	 * 
//	 */
//	public void reset();

}

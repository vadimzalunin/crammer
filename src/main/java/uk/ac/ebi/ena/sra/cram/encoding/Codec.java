package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;

import uk.ac.ebi.ena.sra.cram.format.CramFormatException;

public interface Codec<T> {

	public T read(DataInputStream dis) throws IOException, CramFormatException;

	public void write(DataOutputStream dos, T object) throws IOException;
}

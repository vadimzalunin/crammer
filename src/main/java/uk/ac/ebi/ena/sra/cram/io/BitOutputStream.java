package uk.ac.ebi.ena.sra.cram.io;

import java.io.IOException;

public interface BitOutputStream {

	public void writeByte(int b) throws IOException;

	public void write(int b, int nbits) throws IOException;

	public void write(long b, int nbits) throws IOException;

	public void write(byte b, int nbits) throws IOException;

	public void flush() throws IOException;

}

package uk.ac.ebi.ena.sra.cram.io;

import java.io.IOException;

public interface BitInputStream {

	public boolean readBit() throws IOException;

	public int readBits(int n) throws IOException;
}

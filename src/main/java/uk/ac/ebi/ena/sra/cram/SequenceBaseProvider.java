package uk.ac.ebi.ena.sra.cram;

import java.io.IOException;

public interface SequenceBaseProvider {

	public byte getBaseAt(String sequenceName, long position) throws IOException;

	public void copyBases(String sequenceName, long from, int len, byte[] dest)
			throws IOException;
}

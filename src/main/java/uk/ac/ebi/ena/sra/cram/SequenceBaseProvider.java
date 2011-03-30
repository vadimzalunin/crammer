package uk.ac.ebi.ena.sra.cram;

import java.io.IOException;

public interface SequenceBaseProvider {

	public byte getBaseAt (String sequenceName, long l) throws IOException ;
}

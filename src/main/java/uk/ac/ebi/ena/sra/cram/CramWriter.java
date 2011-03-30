package uk.ac.ebi.ena.sra.cram;

import java.io.IOException;

import uk.ac.ebi.ena.sra.cram.format.CramHeader;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;

public interface CramWriter {
	public void setCram(CramHeader cram) throws IOException;

	public void addRecord(CramRecord record) throws IOException;
	
	public void close () throws IOException ;
}

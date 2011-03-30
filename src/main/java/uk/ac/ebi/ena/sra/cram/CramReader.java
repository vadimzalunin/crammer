package uk.ac.ebi.ena.sra.cram;

import java.util.Iterator;

import uk.ac.ebi.ena.sra.cram.format.CramHeader;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;

public interface CramReader extends Iterator<CramRecord> {
	

	public CramHeader getCram () ;
}

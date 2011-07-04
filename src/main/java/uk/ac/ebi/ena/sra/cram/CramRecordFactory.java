package uk.ac.ebi.ena.sra.cram;

import uk.ac.ebi.ena.sra.cram.format.CramRecord;

public interface CramRecordFactory<T> {

	public CramRecord createCramRecord(T object);
}

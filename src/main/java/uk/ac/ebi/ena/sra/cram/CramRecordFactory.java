package uk.ac.ebi.ena.sra.cram;

import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.mask.PositionMask;

public interface CramRecordFactory<T> {

	public CramRecord createCramRecord(T object);

}

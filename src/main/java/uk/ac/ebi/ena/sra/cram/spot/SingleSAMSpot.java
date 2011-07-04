package uk.ac.ebi.ena.sra.cram.spot;

import java.util.Collections;
import java.util.Iterator;

class SingleSAMSpot implements SAMSpot {
	private SAMRecordHolder record;

	public SingleSAMSpot(SAMRecordHolder record) {
		super();
		this.record = record;
	}

	@Override
	public Iterator<SAMRecordHolder> iterator() {
		return Collections.singletonList(record).iterator();
	}

	@Override
	public boolean isComplete() {
		return record != null;
	}

	@Override
	public String getName() {
		return record.getSamRecord().getReadName();
	}

	@Override
	public int getAlignmentStart() {
		return record.getSamRecord().getAlignmentStart();
	}

	@Override
	public void addRecord(SAMRecordHolder record) {
		if (this.record != null)
			throw new IllegalAccessError(
					"Record already set in single-layout spot.");

		this.record = record;
	}

}

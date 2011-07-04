package uk.ac.ebi.ena.sra.cram.spot;

import net.sf.samtools.SAMRecord;

class SAMRecordHolder {

	private SAMRecord samRecord;
	private long index;

	public SAMRecordHolder(SAMRecord samRecord, long index) {
		super();
		this.samRecord = samRecord;
		this.index = index;
	}

	public SAMRecord getSamRecord() {
		return samRecord;
	}

	public long getIndex() {
		return index;
	}
}

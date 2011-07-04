package uk.ac.ebi.ena.sra.cram.spot;

interface SAMSpot extends Iterable<SAMRecordHolder> {

	public boolean isComplete();

	public String getName();

	public int getAlignmentStart();

	public void addRecord(SAMRecordHolder holder);

}

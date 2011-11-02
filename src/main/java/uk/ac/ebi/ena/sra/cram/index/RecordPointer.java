package uk.ac.ebi.ena.sra.cram.index;

public class RecordPointer extends BitPointer implements
		Comparable<RecordPointer> {
	private long blockStart;
	private long alignmentStart;
	private long recordNumber;

	@Override
	public int compareTo(RecordPointer o) {
		return (int) (alignmentStart - o.alignmentStart);
	}

	public long getBlockStart() {
		return blockStart;
	}

	public void setBlockStart(long blockStart) {
		this.blockStart = blockStart;
	}

	public long getAlignmentStart() {
		return alignmentStart;
	}

	public void setAlignmentStart(long alignmentStart) {
		this.alignmentStart = alignmentStart;
	}

	public long getRecordNumber() {
		return recordNumber;
	}

	public void setRecordNumber(long recordNumber) {
		this.recordNumber = recordNumber;
	}

	@Override
	public String toString() {
		return blockStart + ">" +recordNumber + ":"+ super.toString() + "@" + alignmentStart;
	}
}

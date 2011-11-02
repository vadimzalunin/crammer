package uk.ac.ebi.ena.sra.cram.format;

import java.io.Serializable;
import java.util.Collection;

public class CramRecordBlock implements Serializable {
	private transient Collection<CramRecord> records;
	private String sequenceName;
	private int sequenceLength;
	private long firstRecordPosition;
	private long recordCount;
	private int readLength;
	private CramCompression compression;
	private boolean positiveStrandBasePositionReversed = false;
	private boolean negativeStrandBasePositionReversed = !positiveStrandBasePositionReversed;

	private boolean unmappedReadQualityScoresIncluded = false;
	private boolean substitutionQualityScoresIncluded = true;
	private boolean maskedQualityScoresIncluded = true;

	public String getSequenceName() {
		return sequenceName;
	}

	public void setSequenceName(String sequenceName) {
		this.sequenceName = sequenceName;
	}

	public Collection<CramRecord> getRecords() {
		return records;
	}

	public void setRecords(Collection<CramRecord> records) {
		this.records = records;
	}

	public long getFirstRecordPosition() {
		return firstRecordPosition;
	}

	public void setFirstRecordPosition(long firstRecordPosition) {
		this.firstRecordPosition = firstRecordPosition;
	}

	public long getRecordCount() {
		return recordCount;
	}

	public void setRecordCount(long recordCount) {
		this.recordCount = recordCount;
	}

	public int getReadLength() {
		return readLength;
	}

	public void setReadLength(int readLength) {
		this.readLength = readLength;
	}

	public CramCompression getCompression() {
		return compression;
	}

	public void setCompression(CramCompression compression) {
		this.compression = compression;
	}

	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof CramRecordBlock))
			return false;

		CramRecordBlock foe = (CramRecordBlock) obj;

		if (sequenceName == null)
			if (foe.sequenceName != null)
				return false;
		if (!sequenceName.equals(foe.sequenceName))
			return false;
		if (sequenceLength != foe.sequenceLength)
			return false;
		if (firstRecordPosition != foe.firstRecordPosition)
			return false;
		if (recordCount != foe.recordCount)
			return false;
		if (recordCount != foe.recordCount)
			return false;
		if (readLength != foe.readLength)
			return false;
		if (positiveStrandBasePositionReversed != foe.positiveStrandBasePositionReversed)
			return false;
		if (negativeStrandBasePositionReversed != foe.negativeStrandBasePositionReversed)
			return false;
		if (unmappedReadQualityScoresIncluded != foe.unmappedReadQualityScoresIncluded)
			return false;
		if (substitutionQualityScoresIncluded != foe.substitutionQualityScoresIncluded)
			return false;
		if (maskedQualityScoresIncluded != foe.maskedQualityScoresIncluded)
			return false;

		if (!compression.equals(foe.compression))
			return false;

		return true;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append(getClass().getSimpleName()).append(" [");
		sb.append("sequenceName=").append(sequenceName).append(", ");
		sb.append("sequenceLength=").append(sequenceLength).append(", ");
		sb.append("firstRecordPosition=").append(firstRecordPosition)
				.append(", ");
		sb.append("recordCount=").append(recordCount).append(", ");
		sb.append("readLength=").append(readLength).append(", ");
		sb.append("in-read pos reversed=(")
				.append(positiveStrandBasePositionReversed).append(",")
				.append(negativeStrandBasePositionReversed).append("), ");
		sb.append("unmappedReadQualityScoresIncluded=")
				.append(unmappedReadQualityScoresIncluded).append(", ");
		sb.append("substitutionQualityScoresIncluded=")
				.append(substitutionQualityScoresIncluded).append(", ");
		sb.append("maskedQualityScoresIncluded=")
				.append(maskedQualityScoresIncluded).append(", ");
		sb.append("compression=").append(compression).append("]");
		return sb.toString();
	}

	public boolean isPositiveStrandBasePositionReversed() {
		return positiveStrandBasePositionReversed;
	}

	public void setPositiveStrandBasePositionReversed(
			boolean positiveStrandBasePositionReversed) {
		this.positiveStrandBasePositionReversed = positiveStrandBasePositionReversed;
	}

	public boolean isNegativeStrandBasePositionReversed() {
		return negativeStrandBasePositionReversed;
	}

	public void setNegativeStrandBasePositionReversed(
			boolean negativeStrandBasePositionReversed) {
		this.negativeStrandBasePositionReversed = negativeStrandBasePositionReversed;
	}

	public int getSequenceLength() {
		return sequenceLength;
	}

	public void setSequenceLength(int sequenceLength) {
		this.sequenceLength = sequenceLength;
	}

	public boolean isUnmappedReadQualityScoresIncluded() {
		return unmappedReadQualityScoresIncluded;
	}

	public void setUnmappedReadQualityScoresIncluded(
			boolean unmappedReadQualityScoresIncluded) {
		this.unmappedReadQualityScoresIncluded = unmappedReadQualityScoresIncluded;
	}

	public boolean isSubstitutionQualityScoresIncluded() {
		return substitutionQualityScoresIncluded;
	}

	public void setSubstitutionQualityScoresIncluded(
			boolean substitutionQualityScoresIncluded) {
		this.substitutionQualityScoresIncluded = substitutionQualityScoresIncluded;
	}

	public boolean isMaskedQualityScoresIncluded() {
		return maskedQualityScoresIncluded;
	}

	public void setMaskedQualityScoresIncluded(
			boolean maskedQualityScoresIncluded) {
		this.maskedQualityScoresIncluded = maskedQualityScoresIncluded;
	}

}

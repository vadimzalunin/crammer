package uk.ac.ebi.ena.sra.cram.format;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;

public class CramRecord {

	private long alignmentStart;
	private boolean perfectMatch;
	private boolean negativeStrand;
	private boolean readMapped;

	private long readLength;

	private boolean lastFragment;
	private long recordsToNextFragment;

	private byte[] readBases;
	private byte[] qualityScores;

	private List<ReadFeature> readFeatures;

	public long getAlignmentStart() {
		return alignmentStart;
	}

	public void setAlignmentStart(long alignmentStart) {
		this.alignmentStart = alignmentStart;
	}

	public boolean isPerfectMatch() {
		return perfectMatch;
	}

	public void setPerfectMatch(boolean perfectMatch) {
		this.perfectMatch = perfectMatch;
	}

	public boolean isNegativeStrand() {
		return negativeStrand;
	}

	public void setNegativeStrand(boolean negativeStrand) {
		this.negativeStrand = negativeStrand;
	}

	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof CramRecord))
			return false;

		CramRecord r = (CramRecord) obj;

		if (alignmentStart != r.alignmentStart)
			return false;
		if (negativeStrand != r.negativeStrand)
			return false;
		if (perfectMatch != r.perfectMatch)
			return false;
		if (readMapped != r.readMapped)
			return false;
		if (readLength != r.readLength)
			return false;
		if (lastFragment != r.lastFragment)
			return false;
		if (recordsToNextFragment != r.recordsToNextFragment)
			return false;

		if (!deepEquals(readFeatures, r.readFeatures))
			return false;

		if (!Arrays.equals(readBases, r.readBases))
			return false;
		if (!Arrays.equals(qualityScores, r.qualityScores))
			return false;

		return true;
	}

	private boolean deepEquals(Collection<?> c1, Collection<?> c2) {
		if ((c1 == null || c1.isEmpty()) && (c2 == null || c2.isEmpty()))
			return true;
		if (c1 != null)
			return c1.equals(c2);
		return false;
	}

	@Override
	public String toString() {
		StringBuffer sb = new StringBuffer("[");
		sb.append("readMapped=").append(readMapped);
		sb.append("; alignmentStart=").append(alignmentStart);
		sb.append("; negativeStrand=").append(negativeStrand);
		sb.append("; perfectMatch=").append(perfectMatch);
		
		if (readFeatures != null)
			for (ReadFeature feature : readFeatures)
				sb.append("; ").append(feature.toString());

		if (readBases != null)
			sb.append("; ").append("bases: ").append(new String(readBases));
		if (qualityScores != null)
			sb.append("; ").append("qscores: ")
					.append(new String(qualityScores));

		sb.append("]");
		return sb.toString();
	}

	public long getReadLength() {
		return readLength;
	}

	public void setReadLength(long readLength) {
		this.readLength = readLength;
	}

	public boolean isLastFragment() {
		return lastFragment;
	}

	public void setLastFragment(boolean lastFragment) {
		this.lastFragment = lastFragment;
	}

	public long getRecordsToNextFragment() {
		return recordsToNextFragment;
	}

	public void setRecordsToNextFragment(long recordsToNextFragment) {
		this.recordsToNextFragment = recordsToNextFragment;
	}

	public boolean isReadMapped() {
		return readMapped;
	}

	public void setReadMapped(boolean readMapped) {
		this.readMapped = readMapped;
	}

	public List<ReadFeature> getReadFeatures() {
		return readFeatures;
	}

	public void setReadFeatures(List<ReadFeature> readFeatures) {
		this.readFeatures = readFeatures;
	}

	public byte[] getReadBases() {
		return readBases;
	}

	public void setReadBases(byte[] readBases) {
		this.readBases = readBases;
	}

	public byte[] getQualityScores() {
		return qualityScores;
	}

	public void setQualityScores(byte[] qualityScores) {
		this.qualityScores = qualityScores;
	}
}

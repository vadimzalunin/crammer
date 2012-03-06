package uk.ac.ebi.ena.sra.cram.format;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;

public class CramRecord {

	private Collection<ReadAnnotation> annotations;
	public Collection<ReadTag> tags;

	private long alignmentStart;
	public boolean vendorFiltered = false;
	private boolean negativeStrand;
	private boolean readMapped;

	private long readLength;

	private boolean lastFragment;
	private long recordsToNextFragment = -1;

	private byte[] readBases;
	private byte[] qualityScores;

	private List<ReadFeature> readFeatures;

	private int readGroupID = 0;

	private boolean firstInPair = false;
	private boolean properPair = false;
	private boolean duplicate = false;

	private byte mappingQuality;

	private String sequenceName;
	private String readName;
	public CramRecord next, previous;
	public boolean detached = false;
	public int insertSize;
	public long counter =1 ;

	private Byte flags = null;

	public byte getFlags() {
		if (flags == null) {
			byte b = 0;
			b |= detached ? 1 : 0;
			b <<= 1;
			b |= duplicate ? 1 : 0;
			b <<= 1;
			b |= vendorFiltered ? 1 : 0;
			b <<= 1;
			b |= properPair ? 1 : 0;
			b <<= 1;
			b |= firstInPair ? 1 : 0;
			b <<= 1;
			b |= lastFragment ? 1 : 0;
			b <<= 1;
			b |= negativeStrand ? 1 : 0;
			b <<= 1;
			b |= readMapped ? 1 : 0;
			flags = new Byte(b);
		}
		return flags;
	}

	public void setFlags(byte value) {
		int b = 0xFF & value;
		readMapped = ((b & 1) == 0) ? false : true;
		negativeStrand = ((b & 2) == 0) ? false : true;
		lastFragment = ((b & 4) == 0) ? false : true;
		firstInPair = ((b & 8) == 0) ? false : true;

		properPair = ((b & 16) == 0) ? false : true;
		vendorFiltered = ((b & 32) == 0) ? false : true;
		duplicate = ((b & 64) == 0) ? false : true;
		detached = ((b & 128) == 0) ? false : true;
	}

	public void resetFlags() {
		flags = null;
	}

	public long getAlignmentStart() {
		return alignmentStart;
	}

	public void setAlignmentStart(long alignmentStart) {
		this.alignmentStart = alignmentStart;
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
		if (vendorFiltered != r.vendorFiltered)
			return false;
		if (readMapped != r.readMapped)
			return false;
		if (readLength != r.readLength)
			return false;
		if (lastFragment != r.lastFragment)
			return false;
		if (recordsToNextFragment != r.recordsToNextFragment)
			return false;
		if (firstInPair != r.firstInPair)
			return false;
		if (mappingQuality != r.mappingQuality)
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
		sb.append("; vendorFiltered=").append(vendorFiltered);
		sb.append("; firstInPair=").append(firstInPair);
		sb.append("; mappingQuality=").append(mappingQuality);

		if (readFeatures != null)
			for (ReadFeature feature : readFeatures)
				sb.append("; ").append(feature.toString());

		if (readBases != null)
			sb.append("; ").append("bases: ").append(new String(readBases));
		if (qualityScores != null)
			sb.append("; ").append("qscores: ").append(new String(qualityScores));
		// .append(SAMUtils.phredToFastq(qualityScores));

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

	public Collection<ReadAnnotation> getAnnotations() {
		return annotations;
	}

	public void setAnnotations(Collection<ReadAnnotation> annotations) {
		this.annotations = annotations;
	}

	public int getReadGroupID() {
		return readGroupID;
	}

	public void setReadGroupID(int readGroupID) {
		this.readGroupID = readGroupID;
	}

	public boolean isFirstInPair() {
		return firstInPair;
	}

	public void setFirstInPair(boolean firstInPair) {
		this.firstInPair = firstInPair;
	}

	public byte getMappingQuality() {
		return mappingQuality;
	}

	public void setMappingQuality(byte mappingQuality) {
		this.mappingQuality = mappingQuality;
	}

	public boolean isProperPair() {
		return properPair;
	}

	public void setProperPair(boolean properPair) {
		this.properPair = properPair;
	}

	public boolean isDuplicate() {
		return duplicate;
	}

	public void setDuplicate(boolean duplicate) {
		this.duplicate = duplicate;
	}

	public String getSequenceName() {
		return sequenceName;
	}

	public void setSequenceName(String sequenceName) {
		this.sequenceName = sequenceName;
	}

	public String getReadName() {
		return readName;
	}

	public void setReadName(String readName) {
		this.readName = readName;
	}
}

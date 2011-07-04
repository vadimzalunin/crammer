package uk.ac.ebi.ena.sra.cram.format;

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

	// @Deprecated
	// private Collection<InsertionVariation> insertionVariations;
	// @Deprecated
	// private Collection<SubstitutionVariation> substitutionVariations;
	// @Deprecated
	// private Collection<DeletionVariation> deletionVariations;

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

	// @Deprecated
	// public Collection<InsertionVariation> getInsertionVariations() {
	// return insertionVariations;
	// }
	//
	// @Deprecated
	// public void setInsertionVariations(
	// Collection<InsertionVariation> insertionVariations) {
	// this.insertionVariations = insertionVariations;
	// }
	//
	// @Deprecated
	// public Collection<SubstitutionVariation> getSubstitutionVariations() {
	// return substitutionVariations;
	// }
	//
	// @Deprecated
	// public void setSubstitutionVariations(
	// Collection<SubstitutionVariation> substitutionVariations) {
	// this.substitutionVariations = substitutionVariations;
	// }
	//
	// @Deprecated
	// public Collection<DeletionVariation> getDeletionVariations() {
	// return deletionVariations;
	// }
	//
	// @Deprecated
	// public void setDeletionVariations(
	// Collection<DeletionVariation> deletionVariations) {
	// this.deletionVariations = deletionVariations;
	// }

	// @Deprecated
	// public Collection<ReadBase> getReadBases() {
	// return readBases;
	// }
	//
	// @Deprecated
	// public void setReadBases(Collection<ReadBase> readBases) {
	// this.readBases = readBases;
	// }

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

		// if (!deepEquals(substitutionVariations, r.substitutionVariations))
		// return false;
		// if (!deepEquals(insertionVariations, r.insertionVariations))
		// return false;
		// if (!deepEquals(deletionVariations, r.deletionVariations))
		// return false;
		// if (!deepEquals(readBases, r.readBases))
		// return false;
		if (!deepEquals(readFeatures, r.readFeatures))
			return false;
		

		return true;
	}
	
	private boolean isEqual (byte[] array1, byte[] array2) {
		if (array1 == null && array2 == null) return true ;
		if (array1 != null && array2 != null) {
			
		} 
		return false ;
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
		sb.append("alignmentStart=").append(alignmentStart);
		sb.append("; negativeStrand=").append(negativeStrand);
		sb.append("; perfectMatch=").append(perfectMatch);

		// if (substitutionVariations != null)
		// for (SubstitutionVariation sub : substitutionVariations)
		// sb.append("; ").append(sub.toString());
		// if (insertionVariations != null)
		// for (InsertionVariation ins : insertionVariations)
		// sb.append("; ").append(ins.toString());
		// if (deletionVariations != null)
		// for (DeletionVariation del : deletionVariations)
		// sb.append("; ").append(del.toString());
		// if (readBases != null)
		// for (ReadBase base : readBases)
		// sb.append("; ").append(base.toString());
		if (readFeatures != null)
			for (ReadFeature feature : readFeatures)
				sb.append("; ").append(feature.toString());

		if (readBases != null)
			sb.append("bases: ").append(new String(readBases));
		if (qualityScores != null)
			sb.append("qscores: ").append(new String(qualityScores));

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

package uk.ac.ebi.ena.sra.cram.format;

import java.util.Collection;

public class CramRecord {

	private long alignmentStart;
	private boolean perfectMatch;
	private boolean negativeStrand;

	private int readLength;

	private boolean lastFragment;
	private int recordsToNextFragment;

	private Collection<InsertionVariation> insertionVariations;
	private Collection<SubstitutionVariation> substitutionVariations;
	private Collection<DeletionVariation> deletionVariations;

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

	public Collection<InsertionVariation> getInsertionVariations() {
		return insertionVariations;
	}

	public void setInsertionVariations(
			Collection<InsertionVariation> insertionVariations) {
		this.insertionVariations = insertionVariations;
	}

	public Collection<SubstitutionVariation> getSubstitutionVariations() {
		return substitutionVariations;
	}

	public void setSubstitutionVariations(
			Collection<SubstitutionVariation> substitutionVariations) {
		this.substitutionVariations = substitutionVariations;
	}

	public Collection<DeletionVariation> getDeletionVariations() {
		return deletionVariations;
	}

	public void setDeletionVariations(
			Collection<DeletionVariation> deletionVariations) {
		this.deletionVariations = deletionVariations;
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

		if (!deepEquals(substitutionVariations, r.substitutionVariations))
			return false;
		if (!deepEquals(insertionVariations, r.insertionVariations))
			return false;
		if (!deepEquals(deletionVariations, r.deletionVariations))
			return false;

		return true;
	}

	private boolean deepEquals(Collection<?> c1, Collection<?> c2) {
		if (c1 == null && c2 == null)
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

		if (substitutionVariations != null)
			for (SubstitutionVariation sub : substitutionVariations)
				sb.append("; ").append(sub.toString());
		if (insertionVariations != null)
			for (InsertionVariation ins : insertionVariations)
				sb.append("; ").append(ins.toString());
		if (deletionVariations != null)
			for (DeletionVariation del : deletionVariations)
				sb.append("; ").append(del.toString());

		sb.append("]");
		return sb.toString();
	}

	public int getReadLength() {
		return readLength;
	}

	public void setReadLength(int readLength) {
		this.readLength = readLength;
	}

	public boolean isLastFragment() {
		return lastFragment;
	}

	public void setLastFragment(boolean lastFragment) {
		this.lastFragment = lastFragment;
	}

	public int getRecordsToNextFragment() {
		return recordsToNextFragment;
	}

	public void setRecordsToNextFragment(int recordsToNextFragment) {
		this.recordsToNextFragment = recordsToNextFragment;
	}
}

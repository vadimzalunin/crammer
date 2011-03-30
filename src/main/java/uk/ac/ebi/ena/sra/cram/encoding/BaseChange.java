package uk.ac.ebi.ena.sra.cram.encoding;

import uk.ac.ebi.ena.sra.cram.impl.BaseTransitions;

public class BaseChange {
	private int change;

	public BaseChange(int change) {
		this.change = change;
	}

	public BaseChange(byte from, byte to) {
		change = toInt(from, to);
	}

	public byte getBaseForReference(byte refBase) {
		return BaseTransitions.getBaseForTransition(refBase, change);
	}

	public static int toInt(byte from, byte to) {
		return BaseTransitions.getBaseTransition(from, to);
	}

	public int getChange() {
		return change;
	}

	@Override
	public boolean equals(Object obj) {
		if (obj instanceof BaseChange && ((BaseChange) obj).change == change)
			return true;
		return false;
	}
}

package uk.ac.ebi.ena.sra.cram.format;

import java.io.Serializable;
import java.util.Arrays;

public class InsertionVariation implements Serializable, Variation{

	private int position;
	private byte[] sequence;

	
	@Override
	public byte getOperator() {
		return 'I';
	}
	public int getPosition() {
		return position;
	}

	public void setPosition(int position) {
		this.position = position;
	}

	public byte[] getSequence() {
		return sequence;
	}

	public void setSequence(byte[] sequence) {
		this.sequence = sequence;
	}

	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof InsertionVariation))
			return false;

		InsertionVariation v = (InsertionVariation) obj;

		if (position != v.position)
			return false;
		return Arrays.equals(sequence, v.sequence);
	}

	@Override
	public String toString() {
		StringBuffer sb = new StringBuffer(getClass().getSimpleName() + "[");
		sb.append("position=").append(position);
		sb.append("; sequence=").append(new String(sequence));
		sb.append("] ");
		return sb.toString();
	}
}

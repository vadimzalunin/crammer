package net.sf.cram.encoding.read_features;

import java.io.Serializable;
import java.util.Arrays;

public class Insertion implements Serializable, ReadFeature {

	private int position;
	private byte[] sequence;
	public static final byte operator = 'I';

	public Insertion() {
	}

	public Insertion(int position, byte[] sequence) {
		this.position = position;
		this.sequence = sequence;
	}

	@Override
	public byte getOperator() {
		return operator;
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
		if (!(obj instanceof Insertion))
			return false;

		Insertion v = (Insertion) obj;

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

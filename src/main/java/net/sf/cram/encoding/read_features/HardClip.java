package net.sf.cram.encoding.read_features;

import java.io.Serializable;
import java.util.Arrays;

public class HardClip implements Serializable, ReadFeature {

	private int position;
	private byte[] sequence;

	public byte[] getSequence() {
		return sequence;
	}

	public void setSequence(byte[] sequence) {
		this.sequence = sequence;
	}

	public HardClip() {
	}

	public HardClip(int position, byte[] sequence) {
		this.position = position;
		this.sequence = sequence;
	}

	public static final byte operator = 'H';

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

	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof HardClip))
			return false;

		HardClip v = (HardClip) obj;

		if (position != v.position)
			return false;
		if (Arrays.equals(sequence, v.sequence))
			return false;

		return true;
	}

	@Override
	public String toString() {
		StringBuffer sb = new StringBuffer(getClass().getSimpleName() + "[");
		sb.append("position=").append(position);
		sb.append("; bases=").append(new String(sequence));
		sb.append("] ");
		return sb.toString();
	}
}

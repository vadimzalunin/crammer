package net.sf.cram.encoding.read_features;

import java.io.Serializable;

import net.sf.samtools.CigarOperator;

public class Deletion implements Serializable, ReadFeature{

	private int position;
	private int length;
	public static final byte operator = 'D';

	public Deletion() {
	}

	public Deletion(int position, int length) {
		this.position = position;
		this.length = length;
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

	public int getLength() {
		return length;
	}

	public void setLength(int length) {
		this.length = length;
	}

	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof Deletion))
			return false;

		Deletion v = (Deletion) obj;

		if (position != v.position)
			return false;
		if (length != v.length)
			return false;

		return true;
	}

	@Override
	public String toString() {
		StringBuffer sb = new StringBuffer().append((char)operator).append('@');
		sb.append(position);
		sb.append('+').append(length);
		return sb.toString();
	}
}

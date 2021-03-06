package net.sf.cram.encoding.read_features;

import java.io.Serializable;

public class BaseQualityScore implements Serializable, ReadFeature {

	private int position;
	private byte qualityScore;

	public static final byte operator = 'Q';

	public BaseQualityScore(int position, byte qualityScore) {
		this.position = position;
		this.qualityScore = qualityScore;
	}

	@Override
	public byte getOperator() {
		return operator;
	}

	@Override
	public int getPosition() {
		return position;
	}

	public void setPosition(int position) {
		this.position = position;
	}

	public byte getQualityScore() {
		return qualityScore;
	}

	public void setQualityScore(byte qualityScore) {
		this.qualityScore = qualityScore;
	}

	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof BaseQualityScore))
			return false;

		BaseQualityScore v = (BaseQualityScore) obj;

		if (position != v.position)
			return false;
		if (qualityScore != v.qualityScore)
			return false;

		return true;
	}

	@Override
	public String toString() {
		StringBuffer sb = new StringBuffer().append((char)operator).append('@');
		sb.append(position);
		sb.append('#').appendCodePoint(qualityScore);
		return sb.toString();
	}

}

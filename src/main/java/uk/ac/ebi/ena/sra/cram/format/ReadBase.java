package uk.ac.ebi.ena.sra.cram.format;

import java.io.Serializable;

public class ReadBase implements Serializable, ReadFeature {

	private int position;
	private byte qualityScore;

	public static final byte operator = 'N';

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
		if (!(obj instanceof ReadBase))
			return false;

		ReadBase v = (ReadBase) obj;

		if (position != v.position)
			return false;
		if (qualityScore != v.qualityScore)
			return false;

		return true;
	}

	@Override
	public String toString() {
		StringBuffer sb = new StringBuffer(getClass().getSimpleName() + "[");
		sb.append("position=").append(position);
		sb.append("; score=").append((char) qualityScore);
		sb.append("] ");
		return sb.toString();
	}

}

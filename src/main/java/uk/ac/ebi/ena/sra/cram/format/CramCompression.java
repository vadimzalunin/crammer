package uk.ac.ebi.ena.sra.cram.format;

import java.io.Serializable;
import java.util.Arrays;

import uk.ac.ebi.ena.sra.cram.format.compression.EncodingAlgorithm;

public class CramCompression implements Serializable {
	private Encoding inSeqPosEncoding;
	private Encoding inReadPosEncoding;
	private Encoding readLengthEncoding;
	private Encoding delLengthEncoding;
	private Encoding recordsToNextFragmentEncoding;

	private byte[] baseAlphabet;
	private int[] baseFrequencies;

	private byte[] scoreAlphabet;
	private int[] scoreFrequencies;

	private int[] readLengthAlphabet;
	private int[] readLengthFrequencies;

	private byte[] readFeatureAlphabet;
	private int[] readFeatureFrequencies;

	public CramCompression() {
		setBaseAlphabet("ACGTN".getBytes());
		setBaseFrequencies(new int[] { 100, 90, 80, 75, 10 });

		byte[] scoreAlphabet = new byte[73];
		int[] scoreFreqs = new int[scoreAlphabet.length];
		for (int i = 0; i < scoreAlphabet.length; i++) {
			scoreAlphabet[i] = (byte) (i);
			scoreFreqs[i] = 100 + i;
		}

		setReadLengthAlphabet(new int[] { 36 });
		setReadLengthFrequencies(new int[] { 100 });

		readFeatureAlphabet = "$SIDN".getBytes();
		readFeatureFrequencies = new int[] { 100, 90, 10, 10, 1 };

		setScoreFrequencies(scoreFreqs);

		setDelLengthEncoding(new Encoding(EncodingAlgorithm.GOLOMB_RICE,
				"1,0,1"));
		setInReadPosEncoding(new Encoding(EncodingAlgorithm.GOLOMB_RICE,
				"1,0,1"));
		setInSeqPosEncoding(new Encoding(EncodingAlgorithm.GOLOMB_RICE, "1,0,1"));
		setReadLengthEncoding(new Encoding(EncodingAlgorithm.GOLOMB_RICE,
				"1,0,1"));
		setRecordsToNextFragmentEncoding(new Encoding(
				EncodingAlgorithm.GOLOMB_RICE, "1,0,1"));
	}

	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof CramCompression))
			return false;
		CramCompression foe = (CramCompression) obj;
		if (!getInSeqPosEncoding().equals(foe.getInSeqPosEncoding()))
			return false;
		if (!getInReadPosEncoding().equals(foe.getInReadPosEncoding()))
			return false;
		if (!getReadLengthEncoding().equals(foe.getReadLengthEncoding()))
			return false;
		if (!getDelLengthEncoding().equals(foe.getDelLengthEncoding()))
			return false;
		if (!getRecordsToNextFragmentEncoding().equals(
				foe.getRecordsToNextFragmentEncoding()))
			return false;

		if (!Arrays.equals(getBaseFrequencies(), foe.getBaseFrequencies()))
			return false;
		if (!Arrays.equals(getBaseAlphabet(), foe.getBaseAlphabet()))
			return false;
		if (!Arrays.equals(getScoreFrequencies(), foe.getScoreFrequencies()))
			return false;
		if (!Arrays.equals(getScoreAlphabet(), foe.getScoreAlphabet()))
			return false;
		if (!Arrays.equals(getReadLengthFrequencies(),
				foe.getReadLengthFrequencies()))
			return false;
		if (!Arrays
				.equals(getReadLengthAlphabet(), foe.getReadLengthAlphabet()))
			return false;
		if (!Arrays.equals(getReadFeatureFrequencies(),
				foe.getReadFeatureFrequencies()))
			return false;
		if (!Arrays.equals(getReadFeatureAlphabet(),
				foe.getReadFeatureAlphabet()))
			return false;

		return true;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();

		sb.append("inSeqPosEncoding=").append(getInSeqPosEncoding())
				.append(", ");
		sb.append("inReadPosEncoding=").append(getInReadPosEncoding())
				.append(", ");
		sb.append("readLengthEncoding=").append(getReadLengthEncoding())
				.append(", ");
		sb.append("delLengthEncoding=").append(getDelLengthEncoding())
				.append(", ");
		sb.append("recordsToNextFragmentEncoding=")
				.append(getRecordsToNextFragmentEncoding()).append(", ");

		sb.append("baseAlphabet=").append(Arrays.toString(getBaseAlphabet()))
				.append(", ");
		sb.append("baseFrequencies=")
				.append(Arrays.toString(getBaseFrequencies())).append(", ");

		sb.append("scoreAlphabet=").append(Arrays.toString(getScoreAlphabet()))
				.append(", ");
		sb.append("scoreFrequencies=")
				.append(Arrays.toString(getScoreFrequencies())).append("], ");

		sb.append("readLengthAlphabet=")
				.append(Arrays.toString(getReadLengthAlphabet())).append(", ");
		sb.append("readLengthFrequencies=")
				.append(Arrays.toString(getReadLengthFrequencies()))
				.append("], ");

		sb.append("readFeatureAlphabet=")
				.append(Arrays.toString(getReadFeatureAlphabet())).append(", ");
		sb.append("readFeatureFrequencies=")
				.append(Arrays.toString(getReadFeatureFrequencies()))
				.append("]");

		return sb.toString();
	}

	public void setInSeqPosEncoding(Encoding inSeqPosEncoding) {
		this.inSeqPosEncoding = inSeqPosEncoding;
	}

	public Encoding getInSeqPosEncoding() {
		return inSeqPosEncoding;
	}

	public void setInReadPosEncoding(Encoding inReadPosEncoding) {
		this.inReadPosEncoding = inReadPosEncoding;
	}

	public Encoding getInReadPosEncoding() {
		return inReadPosEncoding;
	}

	public void setReadLengthEncoding(Encoding readLengthEncoding) {
		this.readLengthEncoding = readLengthEncoding;
	}

	public Encoding getReadLengthEncoding() {
		return readLengthEncoding;
	}

	public void setDelLengthEncoding(Encoding delLengthEncoding) {
		this.delLengthEncoding = delLengthEncoding;
	}

	public Encoding getDelLengthEncoding() {
		return delLengthEncoding;
	}

	public void setRecordsToNextFragmentEncoding(
			Encoding recordsToNextFragmentEncoding) {
		this.recordsToNextFragmentEncoding = recordsToNextFragmentEncoding;
	}

	public Encoding getRecordsToNextFragmentEncoding() {
		return recordsToNextFragmentEncoding;
	}

	public void setBaseAlphabet(byte[] baseAlphabet) {
		this.baseAlphabet = baseAlphabet;
	}

	public byte[] getBaseAlphabet() {
		return baseAlphabet;
	}

	public void setBaseFrequencies(int[] baseFrequencies) {
		this.baseFrequencies = baseFrequencies;
	}

	public int[] getBaseFrequencies() {
		return baseFrequencies;
	}

	public void setScoreAlphabet(byte[] scoreAlphabet) {
		this.scoreAlphabet = scoreAlphabet;
	}

	public byte[] getScoreAlphabet() {
		return scoreAlphabet;
	}

	public void setScoreFrequencies(int[] scoreFrequencies) {
		this.scoreFrequencies = scoreFrequencies;
	}

	public int[] getScoreFrequencies() {
		return scoreFrequencies;
	}

	public int[] getReadLengthAlphabet() {
		return readLengthAlphabet;
	}

	public void setReadLengthAlphabet(int[] readLengthAlphabet) {
		this.readLengthAlphabet = readLengthAlphabet;
	}

	public int[] getReadLengthFrequencies() {
		return readLengthFrequencies;
	}

	public void setReadLengthFrequencies(int[] readLengthFrequencies) {
		this.readLengthFrequencies = readLengthFrequencies;
	}

	public byte[] getReadFeatureAlphabet() {
		return readFeatureAlphabet;
	}

	public void setReadFeatureAlphabet(byte[] readFeatureAlphabet) {
		this.readFeatureAlphabet = readFeatureAlphabet;
	}

	public int[] getReadFeatureFrequencies() {
		return readFeatureFrequencies;
	}

	public void setReadFeatureFrequencies(int[] readFeatureFrequencies) {
		this.readFeatureFrequencies = readFeatureFrequencies;
	}
}

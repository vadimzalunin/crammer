package uk.ac.ebi.ena.sra.cram.format;

import java.io.Serializable;
import java.util.Arrays;
import java.util.Map;
import java.util.TreeMap;

import uk.ac.ebi.ena.sra.cram.format.compression.EncodingAlgorithm;

public class CramCompression implements Serializable {
	private Encoding inSeqPosEncoding;
	private Encoding inReadPosEncoding;
	private Encoding readLengthEncoding;
	private Encoding delLengthEncoding;
	private Encoding recordsToNextFragmentEncoding;

	private byte[] baseAlphabet;
	private int[] baseFrequencies;

	private byte[] stopBaseAlphabet;
	private int[] stopBaseFrequencies;

	private byte[] scoreAlphabet;
	private int[] scoreFrequencies;

	private byte[] stopScoreAlphabet;
	private int[] stopScoreFrequencies;

	private int[] readLengthAlphabet;
	private int[] readLengthFrequencies;

	private byte[] readFeatureAlphabet;
	private int[] readFeatureFrequencies;

	private int[] readAnnotationIndexes;
	private int[] readAnnotationFrequencies;

	private int[] readGroupIndexes;
	private int[] readGroupFrequencies;

	private byte[] mappingQualityAlphabet;
	private int[] mappingQualityFrequencies;

	private byte[] heapByteAlphabet;
	private int[] heapByteFrequencies;

	public Map<String, ByteFrequencies> tagByteFrequencyMap = new TreeMap<String, ByteFrequencies>();
	public Map<String, IntFrequencies> tagByteLengthMap = new TreeMap<String, IntFrequencies>();
	public String[] tagKeyAlphabet;
	public int[] tagKeyFrequency;
	
	public ByteFrequencies flagStats ;

	public CramCompression() {
	}

	public static CramCompression createDefaultCramCompression() {
		CramCompression cc = new CramCompression();
		cc.setBaseAlphabet("ACGTN".getBytes());
		cc.setBaseFrequencies(new int[] { 100, 90, 80, 75, 10 });

		byte[] scoreAlphabet = new byte[73];
		int[] scoreFreqs = new int[scoreAlphabet.length];
		for (int i = 0; i < scoreAlphabet.length; i++) {
			scoreAlphabet[i] = (byte) (i);
			scoreFreqs[i] = 100 + i;
		}

		cc.setReadLengthAlphabet(new int[] { 36 });
		cc.setReadLengthFrequencies(new int[] { 100 });

		cc.readFeatureAlphabet = "$SIDN".getBytes();
		cc.readFeatureFrequencies = new int[] { 100, 90, 10, 10, 1 };

		cc.setScoreFrequencies(scoreFreqs);

		cc.setDelLengthEncoding(new Encoding(EncodingAlgorithm.GOLOMB_RICE, "1,0,1"));
		cc.setInReadPosEncoding(new Encoding(EncodingAlgorithm.GOLOMB_RICE, "1,0,1"));
		cc.setInSeqPosEncoding(new Encoding(EncodingAlgorithm.GOLOMB_RICE, "1,0,1"));
		cc.setReadLengthEncoding(new Encoding(EncodingAlgorithm.GOLOMB_RICE, "1,0,1"));
		cc.setRecordsToNextFragmentEncoding(new Encoding(EncodingAlgorithm.GOLOMB_RICE, "1,0,1"));

		return cc;
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
		if (!getRecordsToNextFragmentEncoding().equals(foe.getRecordsToNextFragmentEncoding()))
			return false;

		if (!Arrays.equals(getBaseFrequencies(), foe.getBaseFrequencies()))
			return false;
		if (!Arrays.equals(getBaseAlphabet(), foe.getBaseAlphabet()))
			return false;
		if (!Arrays.equals(getScoreFrequencies(), foe.getScoreFrequencies()))
			return false;
		if (!Arrays.equals(getScoreAlphabet(), foe.getScoreAlphabet()))
			return false;

		if (!Arrays.equals(getStopBaseFrequencies(), foe.getStopBaseFrequencies()))
			return false;
		if (!Arrays.equals(getStopBaseAlphabet(), foe.getStopBaseAlphabet()))
			return false;
		if (!Arrays.equals(getStopScoreFrequencies(), foe.getStopScoreFrequencies()))
			return false;
		if (!Arrays.equals(getStopScoreAlphabet(), foe.getStopScoreAlphabet()))
			return false;

		if (!Arrays.equals(getReadLengthFrequencies(), foe.getReadLengthFrequencies()))
			return false;
		if (!Arrays.equals(getReadLengthAlphabet(), foe.getReadLengthAlphabet()))
			return false;
		if (!Arrays.equals(getReadFeatureFrequencies(), foe.getReadFeatureFrequencies()))
			return false;
		if (!Arrays.equals(getReadFeatureAlphabet(), foe.getReadFeatureAlphabet()))
			return false;

		if (!Arrays.equals(getReadAnnotationFrequencies(), foe.getReadAnnotationFrequencies()))
			return false;
		if (!Arrays.equals(getReadAnnotationIndexes(), foe.getReadAnnotationIndexes()))
			return false;

		if (!Arrays.equals(getReadGroupFrequencies(), foe.getReadGroupFrequencies()))
			return false;
		if (!Arrays.equals(getReadGroupIndexes(), foe.getReadGroupIndexes()))
			return false;

		if (!Arrays.equals(getMappingQualityFrequencies(), foe.getMappingQualityFrequencies()))
			return false;
		if (!Arrays.equals(getMappingQualityAlphabet(), foe.getMappingQualityAlphabet()))
			return false;

		if (!Arrays.equals(getHeapByteFrequencies(), foe.getHeapByteFrequencies()))
			return false;
		if (!Arrays.equals(getHeapByteAlphabet(), foe.getHeapByteAlphabet()))
			return false;

		return true;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();

		sb.append("inSeqPosEncoding=").append(getInSeqPosEncoding()).append("\n");
		sb.append("inReadPosEncoding=").append(getInReadPosEncoding()).append("\n");
		sb.append("readLengthEncoding=").append(getReadLengthEncoding()).append("\n");
		sb.append("delLengthEncoding=").append(getDelLengthEncoding()).append("\n");
		sb.append("recordsToNextFragmentEncoding=").append(getRecordsToNextFragmentEncoding()).append("\n");

		sb.append("baseAlphabet=").append(Arrays.toString(getBaseAlphabet())).append("\n");
		sb.append("baseFrequencies=").append(Arrays.toString(getBaseFrequencies())).append("\n");

		sb.append("scoreAlphabet=").append(Arrays.toString(getScoreAlphabet())).append("\n");
		sb.append("scoreFrequencies=").append(Arrays.toString(getScoreFrequencies())).append("\n");

		sb.append("readLengthAlphabet=").append(Arrays.toString(getReadLengthAlphabet())).append("\n");
		sb.append("readLengthFrequencies=").append(Arrays.toString(getReadLengthFrequencies())).append("\n");

		sb.append("readFeatureAlphabet=").append(Arrays.toString(getReadFeatureAlphabet())).append("\n");
		sb.append("readFeatureFrequencies=").append(Arrays.toString(getReadFeatureFrequencies())).append("\n");

		sb.append("readAnnotationIndexes=").append(Arrays.toString(getReadAnnotationIndexes())).append("\n");
		sb.append("readAnnotationFrequencies=").append(Arrays.toString(getReadAnnotationFrequencies())).append("\n");

		sb.append("mappingQualityAlphabet=").append(Arrays.toString(getMappingQualityAlphabet())).append("\n");
		sb.append("mappingQualityFrequencies=").append(Arrays.toString(getMappingQualityFrequencies())).append("\n");

		sb.append("heapByteAlphabet=").append(Arrays.toString(getHeapByteAlphabet())).append("\n");
		sb.append("heapByteFrequencies=").append(Arrays.toString(getHeapByteFrequencies())).append("\n");

		sb.append("readGroupIndexes=").append(Arrays.toString(getReadGroupIndexes())).append("\n");
		sb.append("readGroupFrequencies=").append(Arrays.toString(getReadGroupFrequencies())).append("\n");

		if (tagKeyAlphabet != null && tagKeyAlphabet.length > 0) {
			sb.append("tagKeyAlphabet=").append(Arrays.toString(tagKeyAlphabet)).append("\n");
			sb.append("tagKeyFrequency=").append(Arrays.toString(tagKeyFrequency)).append("\n");

			for (int i = 0; i < tagKeyAlphabet.length; i++) {
				String tagKeyAndType = tagKeyAlphabet[i];

				sb.append(tagKeyAndType).append("\n");

				ByteFrequencies bf = tagByteFrequencyMap.get(tagKeyAndType);
				sb.append("Bytes:").append(Arrays.toString(bf.getValues())).append("\n");
				sb.append("Byte freqs:").append(Arrays.toString(bf.getFrequencies())).append("\n");

				IntFrequencies lf = tagByteLengthMap.get(tagKeyAndType);
				sb.append("Length:").append(Arrays.toString(lf.getValues())).append("\n");
				sb.append("Length freqs:").append(Arrays.toString(lf.getFrequencies())).append("\n");
			}

		}
		
		sb.append("Flags:").append(Arrays.toString(flagStats.getValues())).append("\n");
		sb.append("Flags freqs:").append(Arrays.toString(flagStats.getFrequencies())).append("\n");
		
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

	public void setRecordsToNextFragmentEncoding(Encoding recordsToNextFragmentEncoding) {
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

	public int[] getReadAnnotationIndexes() {
		return readAnnotationIndexes;
	}

	public void setReadAnnotationIndexes(int[] readAnnotationIndexes) {
		this.readAnnotationIndexes = readAnnotationIndexes;
	}

	public int[] getReadAnnotationFrequencies() {
		return readAnnotationFrequencies;
	}

	public void setReadAnnotationFrequencies(int[] readAnnotationFrequencies) {
		this.readAnnotationFrequencies = readAnnotationFrequencies;
	}

	public int[] getReadGroupIndexes() {
		return readGroupIndexes;
	}

	public void setReadGroupIndexes(int[] readGroupIndexes) {
		this.readGroupIndexes = readGroupIndexes;
	}

	public int[] getReadGroupFrequencies() {
		return readGroupFrequencies;
	}

	public void setReadGroupFrequencies(int[] readGroupFrequencies) {
		this.readGroupFrequencies = readGroupFrequencies;
	}

	public byte[] getMappingQualityAlphabet() {
		return mappingQualityAlphabet;
	}

	public void setMappingQualityAlphabet(byte[] mappingQualityAlphabet) {
		this.mappingQualityAlphabet = mappingQualityAlphabet;
	}

	public int[] getMappingQualityFrequencies() {
		return mappingQualityFrequencies;
	}

	public void setMappingQualityFrequencies(int[] mappingQualityFrequencies) {
		this.mappingQualityFrequencies = mappingQualityFrequencies;
	}

	public byte[] getStopBaseAlphabet() {
		return stopBaseAlphabet;
	}

	public void setStopBaseAlphabet(byte[] stopBaseAlphabet) {
		this.stopBaseAlphabet = stopBaseAlphabet;
	}

	public int[] getStopBaseFrequencies() {
		return stopBaseFrequencies;
	}

	public void setStopBaseFrequencies(int[] stopBaseFrequencies) {
		this.stopBaseFrequencies = stopBaseFrequencies;
	}

	public byte[] getStopScoreAlphabet() {
		return stopScoreAlphabet;
	}

	public void setStopScoreAlphabet(byte[] stopScoreAlphabet) {
		this.stopScoreAlphabet = stopScoreAlphabet;
	}

	public int[] getStopScoreFrequencies() {
		return stopScoreFrequencies;
	}

	public void setStopScoreFrequencies(int[] stopScoreFrequencies) {
		this.stopScoreFrequencies = stopScoreFrequencies;
	}

	public byte[] getHeapByteAlphabet() {
		return heapByteAlphabet;
	}

	public void setHeapByteAlphabet(byte[] heapByteAlphabet) {
		this.heapByteAlphabet = heapByteAlphabet;
	}

	public int[] getHeapByteFrequencies() {
		return heapByteFrequencies;
	}

	public void setHeapByteFrequencies(int[] heapByteFrequencies) {
		this.heapByteFrequencies = heapByteFrequencies;
	}

}

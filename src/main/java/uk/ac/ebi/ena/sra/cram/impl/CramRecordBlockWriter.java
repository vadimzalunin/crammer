package uk.ac.ebi.ena.sra.cram.impl;

import java.io.DataOutputStream;
import java.io.IOException;
import java.io.OutputStream;

import uk.ac.ebi.ena.sra.cram.format.CramCompression;
import uk.ac.ebi.ena.sra.cram.format.CramRecordBlock;
import uk.ac.ebi.ena.sra.cram.format.Encoding;

class CramRecordBlockWriter {
	private OutputStream delegate;

	public CramRecordBlockWriter(OutputStream delegate) {
		super();
		this.delegate = delegate;
	}

	public long write(CramRecordBlock block) throws IOException {
		DataOutputStream dos = new DataOutputStream(delegate);

		dos.write("BLOCKBEGIN".getBytes());

		dos.writeUTF(block.getSequenceName());
		dos.writeInt(block.getSequenceLength());

		dos.writeLong(block.getFirstRecordPosition());
		dos.writeLong(block.getRecordCount());
		dos.writeInt(block.getReadLength());
		dos.writeBoolean(block.isPositiveStrandBasePositionReversed());
		dos.writeBoolean(block.isNegativeStrandBasePositionReversed());
		dos.writeBoolean(block.isUnmappedReadQualityScoresIncluded());
		dos.writeBoolean(block.isSubstitutionQualityScoresIncluded());
		dos.writeBoolean(block.isMaskedQualityScoresIncluded());

		dos.write("COMPRESSIONBEGIN".getBytes());
		writeCramCompression(dos, block.getCompression());
		dos.write("BLOCKEND".getBytes());
		dos.flush();
		return dos.size() * 8;
	}

	private static void writeArray(DataOutputStream dos, int[] array)
			throws IOException {
		dos.writeInt(array.length);
		for (int i = 0; i < array.length; i++)
			dos.writeInt(array[i]);
	}

	private static void writeArray(DataOutputStream dos, byte[] array)
			throws IOException {
		dos.writeInt(array.length);
		for (int i = 0; i < array.length; i++)
			dos.writeByte(array[i]);
	}

	private static final void writeCramCompression(DataOutputStream os,
			CramCompression compression) throws IOException {

		writeEncoding(os, compression.getInSeqPosEncoding());
		writeEncoding(os, compression.getInReadPosEncoding());
		writeEncoding(os, compression.getReadLengthEncoding());
		writeEncoding(os, compression.getDelLengthEncoding());
		writeEncoding(os, compression.getRecordsToNextFragmentEncoding());

		writeArray(os, compression.getBaseAlphabet());
		writeArray(os, compression.getBaseFrequencies());

		writeArray(os, compression.getScoreAlphabet());
		writeArray(os, compression.getScoreFrequencies());

		writeArray(os, compression.getReadFeatureAlphabet());
		writeArray(os, compression.getReadFeatureFrequencies());

		writeArray(os, compression.getReadLengthAlphabet());
		writeArray(os, compression.getReadLengthFrequencies());

		writeArray(os, compression.getReadAnnotationIndexes());
		writeArray(os, compression.getReadAnnotationFrequencies());
		
		writeArray(os, compression.getReadGroupIndexes());
		writeArray(os, compression.getReadGroupFrequencies());
		
		writeArray(os, compression.getMappingQualityAlphabet());
		writeArray(os, compression.getMappingQualityFrequencies());
	}

	private static final void writeEncoding(DataOutputStream os,
			Encoding encoding) throws IOException {
		os.writeByte(encoding.getAlgorithm().ordinal());
		os.writeUTF(encoding.getParameters());
	}
}

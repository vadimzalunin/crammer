package uk.ac.ebi.ena.sra.cram.impl;

import java.io.DataOutputStream;
import java.io.IOException;
import java.io.OutputStream;

import uk.ac.ebi.ena.sra.cram.format.CramCompression;
import uk.ac.ebi.ena.sra.cram.format.CramRecordBlock;
import uk.ac.ebi.ena.sra.cram.format.Encoding;

public class CramRecordBlockWriter {
	private OutputStream delegate;

	public CramRecordBlockWriter(OutputStream delegate) {
		super();
		this.delegate = delegate;
	}

	public long write(CramRecordBlock block) throws IOException {
		DataOutputStream dos = new DataOutputStream(delegate);

		dos.writeUTF(block.getSequenceName());
		dos.writeInt(block.getSequenceLength());

		dos.writeLong(block.getFirstRecordPosition());
		dos.writeLong(block.getRecordCount());
		dos.writeInt(block.getReadLength());
		dos.writeBoolean(block.isPositiveStrandBasePositionReversed());
		dos.writeBoolean(block.isNegativeStrandBasePositionReversed());

		writeCramCompression(dos, block.getCompression());
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
	}

	private static final void writeEncoding(DataOutputStream os,
			Encoding encoding) throws IOException {
		os.writeByte(encoding.getAlgorithm().ordinal());
		os.writeUTF(encoding.getParameters());
	}
}

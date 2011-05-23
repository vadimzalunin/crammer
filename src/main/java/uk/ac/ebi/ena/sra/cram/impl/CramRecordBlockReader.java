package uk.ac.ebi.ena.sra.cram.impl;

import java.io.DataInputStream;
import java.io.EOFException;
import java.io.IOException;

import uk.ac.ebi.ena.sra.cram.format.CramCompression;
import uk.ac.ebi.ena.sra.cram.format.CramRecordBlock;
import uk.ac.ebi.ena.sra.cram.format.Encoding;
import uk.ac.ebi.ena.sra.cram.format.compression.EncodingAlgorithm;

public class CramRecordBlockReader {
	private DataInputStream dis;

	public CramRecordBlockReader(DataInputStream os) {
		super();
		this.dis = os;
	}

	/** Read next block header from the stream 
	 * @return next block header or null if EOF.
	 * @throws IOException
	 */
	public CramRecordBlock read() throws IOException {
		CramRecordBlock block = new CramRecordBlock();
		try {
			block.setSequenceName(dis.readUTF());
		} catch (EOFException e) {
			return null;
		}
		block.setSequenceLength(dis.readInt()) ;
		block.setFirstRecordPosition(dis.readLong());
		block.setRecordCount(dis.readLong());
		block.setReadLength(dis.readInt());
		block.setPositiveStrandBasePositionReversed(dis.readBoolean());
		block.setNegativeStrandBasePositionReversed(dis.readBoolean());

		block.setCompression(readCompression());
		return block;
	}

	private int[] readIntArray() throws IOException {
		int size = dis.readInt();
		int[] array = new int[size];
		for (int i = 0; i < size; i++)
			array[i] = dis.readInt();
		return array;
	}

	private byte[] readByteArray() throws IOException {
		int size = dis.readInt();
		byte[] array = new byte[size];
		for (int i = 0; i < size; i++)
			array[i] = dis.readByte();
		return array;
	}

	private CramCompression readCompression() throws IOException {
		CramCompression compression = new CramCompression();

		compression.setInSeqPosEncoding(readEncoding());
		compression.setInReadPosEncoding(readEncoding());
		compression.setReadLengthEncoding(readEncoding());
		compression.setDelLengthEncoding(readEncoding());
		compression.setRecordsToNextFragmentEncoding(readEncoding());

		compression.setBaseAlphabet(readByteArray());
		compression.setBaseFrequencies(readIntArray());

		compression.setScoreAlphabet(readByteArray());
		compression.setScoreFrequencies(readIntArray());

		compression.setReadFeatureAlphabet(readByteArray());
		compression.setReadFeatureFrequencies(readIntArray());

		compression.setReadLengthAlphabet(readIntArray());
		compression.setReadLengthFrequencies(readIntArray());

		return compression;
	}

	private Encoding readEncoding() throws IOException {
		byte ordinal = dis.readByte();
		Encoding encoding = new Encoding();
		encoding.setAlgorithm(EncodingAlgorithm.values()[ordinal]);
		encoding.setParameters(dis.readUTF());
		return encoding;
	}
}

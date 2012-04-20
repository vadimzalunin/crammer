package uk.ac.ebi.ena.sra.cram.impl;

import java.io.DataInputStream;
import java.io.EOFException;
import java.io.IOException;
import java.util.Arrays;

import uk.ac.ebi.ena.sra.cram.format.ByteFrequencies;
import uk.ac.ebi.ena.sra.cram.format.CramCompression;
import uk.ac.ebi.ena.sra.cram.format.CramFormatException;
import uk.ac.ebi.ena.sra.cram.format.CramRecordBlock;
import uk.ac.ebi.ena.sra.cram.format.DiByteFrequencies;
import uk.ac.ebi.ena.sra.cram.format.Encoding;
import uk.ac.ebi.ena.sra.cram.format.IntFrequencies;
import uk.ac.ebi.ena.sra.cram.format.compression.EncodingAlgorithm;

class CramRecordBlockReader {
	private DataInputStream dis;

	public CramRecordBlockReader(DataInputStream os) {
		super();
		this.dis = os;
	}

	private void ensureExpectedBytes(byte[] expectedBytes) throws CramFormatException, IOException {
		byte[] actuallBeginBytes = new byte[expectedBytes.length];
		dis.readFully(actuallBeginBytes);
		if (!Arrays.equals(expectedBytes, actuallBeginBytes))
			throw new CramFormatException("Expecting bytes " + Arrays.toString(expectedBytes) + " ("
					+ new String(expectedBytes) + ")" + " but got " + Arrays.toString(actuallBeginBytes) + " ("
					+ new String(actuallBeginBytes) + ")");
	}

	/**
	 * Read next block header from the stream
	 * 
	 * @return next block header or null if EOF.
	 * @throws IOException
	 * @throws CramFormatException
	 */
	public CramRecordBlock read() throws IOException, CramFormatException {
		CramRecordBlock block = new CramRecordBlock();
		try {
			ensureExpectedBytes("BLOCKBEGIN".getBytes());
			block.setSequenceName(dis.readUTF());
		} catch (EOFException e) {
			return null;
		}
		block.setSequenceLength(dis.readInt());
		block.setFirstRecordPosition(dis.readLong());
		block.setRecordCount(dis.readLong());
		block.setReadLength(dis.readInt());
		block.setPositiveStrandBasePositionReversed(dis.readBoolean());
		block.setNegativeStrandBasePositionReversed(dis.readBoolean());
		block.setUnmappedReadQualityScoresIncluded(dis.readBoolean());
		block.setSubstitutionQualityScoresIncluded(dis.readBoolean());
		block.setMaskedQualityScoresIncluded(dis.readBoolean());
		block.losslessQualityScores = dis.readBoolean();

		ensureExpectedBytes("COMPRESSIONBEGIN".getBytes());
		block.setCompression(readCompression());
		ensureExpectedBytes("BLOCKEND".getBytes());
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

		compression.setStopBaseAlphabet(readByteArray());
		compression.setStopBaseFrequencies(readIntArray());

		compression.setStopScoreAlphabet(readByteArray());
		compression.setStopScoreFrequencies(readIntArray());

		compression.setReadFeatureAlphabet(readByteArray());
		compression.setReadFeatureFrequencies(readIntArray());

		compression.setReadLengthAlphabet(readIntArray());
		compression.setReadLengthFrequencies(readIntArray());

		compression.setReadAnnotationIndexes(readIntArray());
		compression.setReadAnnotationFrequencies(readIntArray());

		compression.setReadGroupIndexes(readIntArray());
		compression.setReadGroupFrequencies(readIntArray());

		compression.setMappingQualityAlphabet(readByteArray());
		compression.setMappingQualityFrequencies(readIntArray());

		compression.setHeapByteAlphabet(readByteArray());
		compression.setHeapByteFrequencies(readIntArray());

		compression.tagKeyAlphabet = new String[dis.readInt()];
		compression.tagKeyFrequency = new int[compression.tagKeyAlphabet.length];
		compression.tagByteFrequencyMap.clear();
		for (int i = 0; i < compression.tagKeyAlphabet.length; i++) {
			byte[] tagKeyBytes = new byte[4];
			dis.readFully(tagKeyBytes);
			String tagKey = new String(tagKeyBytes);
			compression.tagKeyAlphabet[i] = tagKey;
			compression.tagKeyFrequency[i] = dis.readInt();

			ByteFrequencies bf = new ByteFrequencies(readByteArray(), readIntArray());
			compression.tagByteFrequencyMap.put(tagKey, bf);

			IntFrequencies bl = new IntFrequencies(readIntArray(), readIntArray());
			compression.tagByteLengthMap.put(tagKey, bl);
		}

		compression.flagStats = new ByteFrequencies(readByteArray(), readIntArray());
		
//		compression.score2 = readDiByteFrequencies(dis) ;

		return compression;
	}

	private Encoding readEncoding() throws IOException {
		byte ordinal = dis.readByte();
		Encoding encoding = new Encoding();
		encoding.setAlgorithm(EncodingAlgorithm.values()[ordinal]);
		encoding.setParameters(dis.readUTF());
		return encoding;
	}
	
	private static DiByteFrequencies readDiByteFrequencies (DataInputStream dis) throws IOException {
		int len = dis.readInt() ;
		DiByteFrequencies f = new DiByteFrequencies() ;
		for (int i=0; i<len; i++) {
			byte value1 = dis.readByte() ;
			byte value2 = dis.readByte() ;
			f.add(value1, value2, dis.readInt()) ;
		}
		
		return f ;
	}
}

package uk.ac.ebi.ena.sra.cram;

import java.io.ByteArrayInputStream;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.zip.GZIPInputStream;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.SeekableStream;
import uk.ac.ebi.ena.sra.cram.CramIndexer.CountingInputStream;
import uk.ac.ebi.ena.sra.cram.format.CramHeader;
import uk.ac.ebi.ena.sra.cram.format.CramReadGroup;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.CramReferenceSequence;
import uk.ac.ebi.ena.sra.cram.format.ReadFeature;

public class Utils {

	public final static byte[] toBytes(int value) {
		final byte[] bytes = new byte[4];
		bytes[0] = (byte) (value >>> 24);
		bytes[1] = (byte) (value >>> 16);
		bytes[2] = (byte) (value >>> 8);
		bytes[3] = (byte) (value >>> 0);
		return bytes;
	}

	public final static byte[] toBytes(long value) {
		final byte[] bytes = new byte[8];
		for (int i = 0; i < 8; i++)
			bytes[i] = (byte) (value >>> (64 - 8 - i * 8));
		return bytes;
	}

	public static String toBitString(final byte[] b) {
		final char[] bits = new char[8 * b.length];
		for (int i = 0; i < b.length; i++) {
			final byte byteval = b[i];
			int bytei = i << 3;
			int mask = 0x1;
			for (int j = 7; j >= 0; j--) {
				final int bitval = byteval & mask;
				if (bitval == 0) {
					bits[bytei + j] = '0';
				} else {
					bits[bytei + j] = '1';
				}
				mask <<= 1;
			}
		}
		return String.valueOf(bits);
	}

	public static String toBitString(final int value) {
		return toBitString(toBytes(value));
	}

	public static byte[] transformSequence(byte[] bases, boolean compliment, boolean reverse) {
		byte[] result = new byte[bases.length];
		for (int i = 0; i < bases.length; i++) {
			byte base = bases[i];

			int index = reverse ? bases.length - i - 1 : i;

			result[index] = compliment ? complimentBase(base) : base;
		}
		return result;
	}

	public static final byte complimentBase(byte base) {
		switch (base) {
		case 'A':
			return 'T';
		case 'C':
			return 'G';
		case 'G':
			return 'C';
		case 'T':
			return 'A';
		case 'N':
			return 'N';

		default:
			throw new RuntimeException("Unkown base: " + base);
		}
	}

	public static int mostSignificantBit(final long value) {
		int i = 64;
		while (--i >= 0 && (((1L << i) & value)) == 0)
			;
		return i;
	}

	public static Byte[] autobox(byte[] array) {
		Byte[] newArray = new Byte[array.length];
		for (int i = 0; i < array.length; i++)
			newArray[i] = array[i];
		return newArray;
	}

	public static Integer[] autobox(int[] array) {
		Integer[] newArray = new Integer[array.length];
		for (int i = 0; i < array.length; i++)
			newArray[i] = array[i];
		return newArray;
	}

	public static void changeReadLength(SAMRecord record, int newLength) {
		if (newLength == record.getReadLength())
			return;
		if (newLength < 1 || newLength >= record.getReadLength())
			throw new IllegalArgumentException("Cannot change read length to " + newLength);

		List<CigarElement> newCigarElements = new ArrayList<CigarElement>();
		int len = 0;
		for (CigarElement ce : record.getCigar().getCigarElements()) {
			switch (ce.getOperator()) {
			case D:
				break;
			case S:
				// dump = true;
				// len -= ce.getLength();
				// break;
			case M:
			case I:
			case X:
				len += ce.getLength();
				break;

			default:
				throw new IllegalArgumentException("Unexpected cigar operator: " + ce.getOperator() + " in cigar "
						+ record.getCigarString());
			}

			if (len <= newLength) {
				newCigarElements.add(ce);
				continue;
			}
			CigarElement newCe = new CigarElement(ce.getLength() - (record.getReadLength() - newLength),
					ce.getOperator());
			if (newCe.getLength() > 0)
				newCigarElements.add(newCe);
			break;
		}

		byte[] newBases = new byte[newLength];
		System.arraycopy(record.getReadBases(), 0, newBases, 0, newLength);
		record.setReadBases(newBases);

		byte[] newScores = new byte[newLength];
		System.arraycopy(record.getBaseQualities(), 0, newScores, 0, newLength);

		record.setCigar(new Cigar(newCigarElements));
	}

	public static void reversePositionsInRead(CramRecord record) {
		if (record.getReadFeatures() == null || record.getReadFeatures().isEmpty())
			return;
		for (ReadFeature f : record.getReadFeatures())
			f.setPosition((int) (record.getReadLength() - f.getPosition() - 1));

		Collections.reverse(record.getReadFeatures());
	}

	public static byte[] getBasesFromReferenceFile(String referenceFilePath, String seqName, int from, int length) {
		ReferenceSequenceFile referenceSequenceFile = ReferenceSequenceFileFactory.getReferenceSequenceFile(new File(
				referenceFilePath));
		ReferenceSequence sequence = referenceSequenceFile.getSequence(seqName);
		byte[] bases = referenceSequenceFile.getSubsequenceAt(sequence.getName(), from, from + length).getBases();
		return bases;
	}

	public static void capitaliseAndCheckBases(byte[] bases, boolean strict) {
		for (int i = 0; i < bases.length; i++) {
			switch (bases[i]) {
			case 'A':
			case 'C':
			case 'G':
			case 'T':
			case 'N':
				break;
			case 'a':
				bases[i] = 'A';
				break;
			case 'c':
				bases[i] = 'C';
				break;
			case 'g':
				bases[i] = 'G';
				break;
			case 't':
				bases[i] = 'T';
				break;
			case 'n':
				bases[i] = 'N';
				break;

			default:
				if (strict)
					throw new RuntimeException("Illegal base at " + i + ": " + bases[i]);
				else
					bases[i] = 'N';
				break;
			}
		}
	}

	private static int readInt(InputStream in) throws IOException {
		int ch1 = in.read();
		int ch2 = in.read();
		int ch3 = in.read();
		int ch4 = in.read();
		if ((ch1 | ch2 | ch3 | ch4) < 0)
			throw new EOFException();
		return ((ch1 << 24) + (ch2 << 16) + (ch3 << 8) + (ch4 << 0));
	}

	public static final DataInputStream getNextChunk(DataInputStream dis) throws IOException {
		int compressedBlockSize;
		try {
			compressedBlockSize = readInt(dis);
			// compressedBlockSize = dis.readInt();
			// System.out.println("Compressed block sise: " +
			// compressedBlockSize);
		} catch (EOFException e) {
			return null;
		}
		byte[] compressedBlockData = new byte[compressedBlockSize];
		try {
			dis.readFully(compressedBlockData);
		} catch (EOFException e) {
			byte[] buf = new byte[6];
			System.arraycopy(compressedBlockData, 0, buf, 0, buf.length);
			System.err.println("Offensive data block start: " + Arrays.toString(buf));
			throw e;
		}
		GZIPInputStream gizIS = new GZIPInputStream(new ByteArrayInputStream(compressedBlockData));
		CountingInputStream uncompressedCIS = new CountingInputStream(gizIS);
		DataInputStream uncompressedDIS = new DataInputStream(uncompressedCIS);
		return uncompressedDIS;
	}

	public static SAMFileHeader cramHeader2SamHeader(CramHeader cramHeader) {
		SAMFileHeader samFileHeader = new SAMFileHeader();
		for (CramReferenceSequence crs : cramHeader.getReferenceSequences()) {
			SAMSequenceRecord samSequence = new SAMSequenceRecord(crs.getName(), crs.getLength());
			samFileHeader.addSequence(samSequence);
		}

		if (cramHeader.getReadGroups() != null)
			for (CramReadGroup crg : cramHeader.getReadGroups()) {
				if (crg.getId() == null)
					continue;
				SAMReadGroupRecord samReadGroupRecord = new SAMReadGroupRecord(crg.getId());
				samReadGroupRecord.setSample(crg.getSample());
				samFileHeader.addReadGroup(samReadGroupRecord);
			}

		return samFileHeader;
	}

	/**
	 * Reads and discards first bytes which must correspond to the CRAM magic.
	 * 
	 * @param is
	 * @return
	 * @throws IOException
	 */
	public static boolean isCRAM(InputStream is) throws IOException {
		byte[] magick = "CRAM".getBytes();
		for (byte b : magick)
			if (is.read() != b)
				return false;

		return true;
	}

	public static String peekStream(SeekableStream ss, long pos, int len) throws IOException {
		ss.mark(len);
		ss.seek(pos);
		byte[] buf = new byte[len];
		ss.read(buf);
		ss.reset();
		return Arrays.toString(buf);
	}
}

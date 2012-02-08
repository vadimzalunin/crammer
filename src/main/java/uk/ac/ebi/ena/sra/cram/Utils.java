package uk.ac.ebi.ena.sra.cram;

import java.io.ByteArrayInputStream;
import java.io.DataInputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.regex.Pattern;
import java.util.zip.GZIPInputStream;

import net.sf.picard.PicardException;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.picard.sam.SamPairUtil;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMReadGroupRecord;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.SAMTag;
import net.sf.samtools.util.SeekableStream;

import org.apache.log4j.Logger;

import uk.ac.ebi.ena.sra.cram.CramIndexer.CountingInputStream;
import uk.ac.ebi.ena.sra.cram.format.CramHeader;
import uk.ac.ebi.ena.sra.cram.format.CramReadGroup;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.CramReferenceSequence;
import uk.ac.ebi.ena.sra.cram.format.ReadFeature;

public class Utils {
	private static Logger log = Logger.getLogger(Utils.class);

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

	/**
	 * Copied from net.sf.picard.sam.SamPairUtil. This is a more permissive
	 * version of the method, which does not reset alignment start and reference
	 * for unmapped reads.
	 * 
	 * @param rec1
	 * @param rec2
	 * @param header
	 */
	public static void setLooseMateInfo(final SAMRecord rec1, final SAMRecord rec2, final SAMFileHeader header) {
		// If neither read is unmapped just set their mate info
		if (!rec1.getReadUnmappedFlag() && !rec2.getReadUnmappedFlag()) {

			rec1.setMateReferenceIndex(rec2.getReferenceIndex());
			rec1.setMateAlignmentStart(rec2.getAlignmentStart());
			rec1.setMateNegativeStrandFlag(rec2.getReadNegativeStrandFlag());
			rec1.setMateUnmappedFlag(false);
			rec1.setAttribute(SAMTag.MQ.name(), rec2.getMappingQuality());

			rec2.setMateReferenceIndex(rec1.getReferenceIndex());
			rec2.setMateAlignmentStart(rec1.getAlignmentStart());
			rec2.setMateNegativeStrandFlag(rec1.getReadNegativeStrandFlag());
			rec2.setMateUnmappedFlag(false);
			rec2.setAttribute(SAMTag.MQ.name(), rec1.getMappingQuality());
		}
		// Else if they're both unmapped set that straight
		else if (rec1.getReadUnmappedFlag() && rec2.getReadUnmappedFlag()) {
			rec1.setMateNegativeStrandFlag(rec2.getReadNegativeStrandFlag());
			rec1.setMateUnmappedFlag(true);
			rec1.setAttribute(SAMTag.MQ.name(), null);
			rec1.setInferredInsertSize(0);

			rec2.setMateNegativeStrandFlag(rec1.getReadNegativeStrandFlag());
			rec2.setMateUnmappedFlag(true);
			rec2.setAttribute(SAMTag.MQ.name(), null);
			rec2.setInferredInsertSize(0);
		}
		// And if only one is mapped copy it's coordinate information to the
		// mate
		else {
			final SAMRecord mapped = rec1.getReadUnmappedFlag() ? rec2 : rec1;
			final SAMRecord unmapped = rec1.getReadUnmappedFlag() ? rec1 : rec2;

			mapped.setMateReferenceIndex(unmapped.getReferenceIndex());
			mapped.setMateAlignmentStart(unmapped.getAlignmentStart());
			mapped.setMateNegativeStrandFlag(unmapped.getReadNegativeStrandFlag());
			mapped.setMateUnmappedFlag(true);
			mapped.setInferredInsertSize(0);

			unmapped.setMateReferenceIndex(mapped.getReferenceIndex());
			unmapped.setMateAlignmentStart(mapped.getAlignmentStart());
			unmapped.setMateNegativeStrandFlag(mapped.getReadNegativeStrandFlag());
			unmapped.setMateUnmappedFlag(false);
			unmapped.setInferredInsertSize(0);
		}

		final int insertSize = SamPairUtil.computeInsertSize(rec1, rec2);
		rec1.setInferredInsertSize(insertSize);
		rec2.setInferredInsertSize(-insertSize);
	}

	public static IndexedFastaSequenceFile createIndexedFastaSequenceFile(File file) throws CramException,
			FileNotFoundException {
		if (IndexedFastaSequenceFile.canCreateIndexedFastaReader(file)) {
			IndexedFastaSequenceFile ifsFile = new IndexedFastaSequenceFile(file);

			return ifsFile;
		} else
			throw new CramException(
					"Reference fasta file is not indexed or index file not found. Try executing 'samtools faidx "
							+ file.getAbsolutePath() + "'");
	}

	public static ReferenceSequence getReferenceSequenceOrNull(ReferenceSequenceFile rsFile, String name) {
		ReferenceSequence rs = null;
		try {
			return rsFile.getSequence(name);
		} catch (PicardException e) {
			return null;
		}
	}

	private static final Pattern chrPattern = Pattern.compile("chr.*", Pattern.CASE_INSENSITIVE);

	public static byte[] getBasesOrNull(ReferenceSequenceFile rsFile, String name, int start, int len) {
		ReferenceSequence rs = getReferenceSequenceOrNull(rsFile, name);
		if (rs == null && name.equals("M")) {
			rs = getReferenceSequenceOrNull(rsFile, "MT");
		}

		if (rs == null && name.equals("MT")) {
			rs = getReferenceSequenceOrNull(rsFile, "M");
		}

		boolean chrPatternMatch = chrPattern.matcher(name).matches();
		if (rs == null) {
			if (chrPatternMatch)
				rs = getReferenceSequenceOrNull(rsFile, name.substring(3));
			else
				rs = getReferenceSequenceOrNull(rsFile, "chr" + name);
		}

		if (len < 1)
			return rs.getBases();
		else
			return rsFile.getSubsequenceAt(rs.getName(), 1, len).getBases();
	}

	public static byte[] getReferenceSequenceBases(ReferenceSequenceFile referenceSequenceFile, String seqName)
			throws CramException {
		long time1 = System.currentTimeMillis();
		byte[] refBases = Utils.getBasesOrNull(referenceSequenceFile, seqName, 1, 0);
		if (refBases == null)
			throw new CramException("Reference sequence " + seqName + " not found in the fasta file "
					+ referenceSequenceFile.toString());
		
		long time2 = System.currentTimeMillis();
		log.debug(String.format("Reference sequence %s read in %.2f seconds.", seqName, (time2 - time1) / 1000f));
		
		Utils.capitaliseAndCheckBases(refBases, false);
		
		long time3 = System.currentTimeMillis();
		log.debug(String.format("Reference sequence normalized in %.2f seconds.", (time3 - time2) / 1000f));
		return refBases;
	}
}

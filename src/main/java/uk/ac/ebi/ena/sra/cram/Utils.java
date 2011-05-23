package uk.ac.ebi.ena.sra.cram;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMRecord;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.DeletionVariation;
import uk.ac.ebi.ena.sra.cram.format.InsertionVariation;
import uk.ac.ebi.ena.sra.cram.format.ReadFeature;
import uk.ac.ebi.ena.sra.cram.format.SubstitutionVariation;

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

	public static byte[] transformSequence(byte[] bases, boolean compliment,
			boolean reverse) {
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

	public static String toString(String readName, byte[] readBases,
			CramRecord cramRecord) {
		// ERR005143.135209 342 1 CCCGCCCACGGCTGCTCCAGCTGCTGCTGTAGCGAC [(2, 'S',
		// 'ct', None), (4, 'I', 'c', None), (35, 'D', 2, None)]
		StringBuilder sb = new StringBuilder();
		sb.append(readName);
		sb.append("\t");
		sb.append(cramRecord.getAlignmentStart());
		sb.append("\t");
		sb.append(cramRecord.isNegativeStrand() ? 1 : 0);
		sb.append("\t");
		sb.append(new String(readBases));

		sb.append("\t[");
		boolean firstVariation = true;
		for (ReadFeature v : sortVariationsByPosition(cramRecord)) {
			if (!firstVariation)
				sb.append(", ");
			sb.append("(");
			sb.append(v.getPosition() - 1);
			sb.append(", ");
			if (v instanceof SubstitutionVariation) {
				SubstitutionVariation sv = (SubstitutionVariation) v;
				sb.append("'S'");
				sb.append(", ");
				sb.append("'")
						.append(Character.toLowerCase((char) sv.getBase()))
						.append(Character.toLowerCase((char) sv
								.getRefernceBase())).append("'");
				sb.append(", ");
				if (sv.getQualityScore() > -1)
					sb.append(sv.getQualityScore()).append(")");
				else
					sb.append("None)");
			} else if (v instanceof InsertionVariation) {
				InsertionVariation iv = (InsertionVariation) v;
				sb.append("'I'");
				sb.append(", ");
				sb.append("'")
						.append(new String(iv.getSequence()).toLowerCase())
						.append("'");
				sb.append(", ");
				sb.append("None)");
			} else if (v instanceof DeletionVariation) {
				DeletionVariation dv = (DeletionVariation) v;
				sb.append("'D'");
				sb.append(", ");
				sb.append(dv.getLength());
				sb.append(", ");
				sb.append("None)");
			}

			firstVariation = false;
		}
		sb.append("]");

		return sb.toString();
	}

	public static List<ReadFeature> sortVariationsByPosition(CramRecord record) {
		List<ReadFeature> list = new ArrayList<ReadFeature>();
		if (record.getSubstitutionVariations() != null
				&& !record.getSubstitutionVariations().isEmpty())
			list.addAll(record.getSubstitutionVariations());

		if (record.getInsertionVariations() != null
				&& !record.getInsertionVariations().isEmpty())
			list.addAll(record.getInsertionVariations());

		if (record.getDeletionVariations() != null
				&& !record.getDeletionVariations().isEmpty())
			list.addAll(record.getDeletionVariations());

		Collections.sort(list, variationPositionComparator);

		return list;
	}

	private static Comparator<ReadFeature> variationPositionComparator = new Comparator<ReadFeature>() {

		@Override
		public int compare(ReadFeature o1, ReadFeature o2) {
			int result = o1.getPosition() - o2.getPosition();
			if (result != 0)
				return result;

			return o1.getOperator() - o2.getOperator();
		}
	};

	public static byte[] restoreBases(CramRecord record,
			SequenceBaseProvider provider, String seqName) throws IOException {
		byte[] bases = new byte[(int) record.getReadLength()];

		int posInRead = 1;
		long posInSeq = record.getAlignmentStart() - 1;
		List<ReadFeature> variations = sortVariationsByPosition(record);
		variations = record.getReadFeatures();
		for (ReadFeature v : variations) {
			for (; posInRead <= v.getPosition(); posInRead++)
				bases[posInRead - 1] = provider.getBaseAt(seqName, posInSeq++);

			switch (v.getOperator()) {
			case 'S':
				SubstitutionVariation sv = (SubstitutionVariation) v;
				// byte refBase = provider.getBaseAt(seqName, ++posInSeq) ;
				// byte base = sv.getBaseChange().getBaseForReference(refBase) ;
				// sv.setBase(base) ;
				// sv.setRefernceBase(refBase) ;
				bases[posInRead++ - 1] = sv.getBase();
				posInSeq++;
				break;
			case 'I':
				InsertionVariation iv = (InsertionVariation) v;
				for (int i = 0; i < iv.getSequence().length; i++)
					bases[posInRead++ - 1] = iv.getSequence()[i];

				break;
			case 'D':
				DeletionVariation dv = (DeletionVariation) v;
				posInSeq += dv.getLength();
				break;

			default:
				throw new RuntimeException("Uknown variation operator: "
						+ v.getOperator());
			}
		}
		for (; posInRead <= record.getReadLength(); posInRead++)
			bases[posInRead - 1] = provider.getBaseAt(seqName, posInSeq++);
		return bases;
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
			throw new IllegalArgumentException("Cannot change read length to "
					+ newLength);

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
				throw new IllegalArgumentException(
						"Unexpected cigar operator: " + ce.getOperator()
								+ " in cigar " + record.getCigarString());
			}

			if (len <= newLength) {
				newCigarElements.add(ce);
				continue;
			}
			CigarElement newCe = new CigarElement(ce.getLength()
					- (record.getReadLength() - newLength), ce.getOperator());
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
		if (record.getReadFeatures() == null
				|| record.getReadFeatures().isEmpty())
			return;
		for (ReadFeature f : record.getReadFeatures())
			f.setPosition((int) (record.getReadLength() - f.getPosition() - 1));

		Collections.reverse(record.getReadFeatures());
	}

	public static byte[] getBasesFromReferenceFile(String referenceFilePath,
			String seqName, int from, int length) {
		ReferenceSequenceFile referenceSequenceFile = ReferenceSequenceFileFactory
				.getReferenceSequenceFile(new File(referenceFilePath));
		ReferenceSequence sequence = referenceSequenceFile.getSequence(seqName);
		byte[] bases = referenceSequenceFile.getSubsequenceAt(
				sequence.getName(), from, from + length).getBases();
		return bases;
	}
}

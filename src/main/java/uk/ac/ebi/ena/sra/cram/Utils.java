package uk.ac.ebi.ena.sra.cram;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.DeletionVariation;
import uk.ac.ebi.ena.sra.cram.format.InsertionVariation;
import uk.ac.ebi.ena.sra.cram.format.SubstitutionVariation;
import uk.ac.ebi.ena.sra.cram.format.Variation;

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
		sb.append(cramRecord.getAlignmentStart() - 1);
		sb.append("\t");
		sb.append(cramRecord.isNegativeStrand() ? 1 : 0);
		sb.append("\t");
		sb.append(new String(readBases));

		sb.append("\t[");
		boolean firstVariation = true;
		for (Variation v : sortVariationsByPosition(cramRecord)) {
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
			} else if (v instanceof InsertionVariation) {
				InsertionVariation iv = (InsertionVariation) v;
				sb.append("'I'");
				sb.append(", ");
				sb.append("'")
						.append(new String(iv.getSequence()).toLowerCase())
						.append("'");
			} else if (v instanceof DeletionVariation) {
				DeletionVariation dv = (DeletionVariation) v;
				sb.append("'D'");
				sb.append(", ");
				sb.append(dv.getLength());
			}
			sb.append(", ");
			sb.append("None)");
			firstVariation = false;
		}
		sb.append("]");

		return sb.toString();
	}

	public static List<Variation> sortVariationsByPosition(CramRecord record) {
		List<Variation> list = new ArrayList<Variation>();
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

	private static Comparator<Variation> variationPositionComparator = new Comparator<Variation>() {

		@Override
		public int compare(Variation o1, Variation o2) {
			int result = o1.getPosition() - o2.getPosition();
			if (result != 0)
				return result;

			return o1.getOperator() - o2.getOperator();
		}
	};

	public static byte[] restoreBases(CramRecord record,
			SequenceBaseProvider provider, String seqName) throws IOException {
		byte[] bases = new byte[record.getReadLength()];

		int posInRead = 1;
		long posInSeq = record.getAlignmentStart() - 1;
		List<Variation> variations = sortVariationsByPosition(record);
		for (Variation v : variations) {
			for (; posInRead < v.getPosition(); posInRead++)
				bases[posInRead - 1] = provider.getBaseAt(seqName, posInSeq++);

			switch (v.getOperator()) {
			case 'S':
				SubstitutionVariation sv = (SubstitutionVariation) v;
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
}

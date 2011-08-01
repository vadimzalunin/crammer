package uk.ac.ebi.ena.sra.cram.format.text;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import uk.ac.ebi.ena.sra.cram.encoding.BaseChange;
import uk.ac.ebi.ena.sra.cram.format.BaseQualityScore;
import uk.ac.ebi.ena.sra.cram.format.DeletionVariation;
import uk.ac.ebi.ena.sra.cram.format.InsertBase;
import uk.ac.ebi.ena.sra.cram.format.InsertionVariation;
import uk.ac.ebi.ena.sra.cram.format.ReadBase;
import uk.ac.ebi.ena.sra.cram.format.ReadFeature;
import uk.ac.ebi.ena.sra.cram.format.SubstitutionVariation;

class DefaultReadFeaturesFormat implements ReadFeaturesFormat {
	public static final byte MATCH = 'M';
	public static final byte DELETION = 'D';
	public static final byte INSERTION = 'I';
	public static final byte INSERT_BASE = 'i';
	public static final byte SUBSTITUTION = 'S';
	public static final byte SUBSTITUTION_0 = 'z';
	public static final byte SUBSTITUTION_1 = 'y';
	public static final byte SUBSTITUTION_2 = 'x';
	public static final byte SUBSTITUTION_3 = 'w';
	public static final byte BASE_A = 'A';
	public static final byte BASE_C = 'C';
	public static final byte BASE_G = 'G';
	public static final byte BASE_T = 'T';
	public static final byte BASE_N = 'N';
	public static final byte QUALITY_SCORE = 'Q';
	public static final byte UNSET_QUALITY_SCORE = 32;

	public static final byte INSERT_A = 'a';
	public static final byte INSERT_C = 'c';
	public static final byte INSERT_G = 'g';
	public static final byte INSERT_T = 't';
	public static final byte INSERT_N = 'n';

	private boolean useOneLetterCodesForSubstitutions = true;

	public static void main(String[] args) {
		test ("M1", 1) ;
		test ("M2", 2) ;
		test ("M2D1M1", 3) ;
		test ("D1M2D1", 2) ;
		test ("a ", 1) ;
		test ("a c!", 2) ;
		test ("M1a c!M1", 4) ;
		test ("M1y!D1x a c!D1M1", 6) ;
	}

	private static void test(String s, int readLength) {
		DefaultReadFeaturesFormat format = new DefaultReadFeaturesFormat();
		List<ReadFeature> featureList = format.asReadFeatureList(s);
		StringBuilder sb = new StringBuilder();
		format.addFeaturesToStringBuilder(featureList, readLength, sb);
		if (s.equals(sb.toString()))
			System.out.println("ok: " + s);
		else {
			System.out.println("failed: " + s + "\t\t\t" + sb.toString());
			for (ReadFeature f:featureList)
				System.out.println(f.toString());
		}
	}

	@Override
	public void addFeaturesToStringBuilder(Iterable<ReadFeature> features,
			int readLength, StringBuilder sb) {
		int prevPos = 1;
		boolean qsRequired = false ;
		for (ReadFeature feature : features) {
			if (feature.getOperator() != BaseQualityScore.operator && qsRequired) {
				sb.appendCodePoint(UNSET_QUALITY_SCORE) ;
				qsRequired = false ;
			}
			if (feature.getPosition() > prevPos) {
				sb.appendCodePoint(MATCH).append(
						feature.getPosition() - prevPos);
				prevPos = feature.getPosition();
			}
			switch (feature.getOperator()) {
			case ReadBase.operator:
				ReadBase rb = (ReadBase) feature;
				sb.appendCodePoint(rb.getBase()).appendCodePoint(
						rb.getQualityScore() == -1 ? Byte.MAX_VALUE : rb.getQualityScore());
				prevPos++;
				break;
			case DeletionVariation.operator:
				DeletionVariation del = (DeletionVariation) feature;
				sb.appendCodePoint(DELETION).append(del.getLength());
				break;
			case InsertionVariation.operator:
				InsertionVariation ins = (InsertionVariation) feature;
				sb.appendCodePoint(INSERTION)
						.append(new String(ins.getSequence())).append(".");
				prevPos += ins.getSequence().length;
				break;
			case InsertBase.operator:
				InsertBase ibs = (InsertBase) feature ;
				sb.appendCodePoint(Character.toLowerCase(ibs.getBase())) ;
				qsRequired = true ;
				prevPos++;
				break ;
			case BaseQualityScore.operator:
				BaseQualityScore bqs = (BaseQualityScore) feature ;
				sb.appendCodePoint(bqs.getQualityScore()) ;
				qsRequired = false ;
				break ;
			case SubstitutionVariation.operator:
				SubstitutionVariation sub = (SubstitutionVariation) feature;
				if (useOneLetterCodesForSubstitutions) {
					switch (sub.getBaseChange().getChange()) {
					case 0:
						sb.appendCodePoint(SUBSTITUTION_0);
						break;
					case 1:
						sb.appendCodePoint(SUBSTITUTION_1);
						break;
					case 2:
						sb.appendCodePoint(SUBSTITUTION_2);
						break;
					case 3:
						sb.appendCodePoint(SUBSTITUTION_3);
						break;

					default:
						break;
					}
				} else
					sb.appendCodePoint(SUBSTITUTION).append(
							sub.getBaseChange().getChange());

//				sb.appendCodePoint(sub.getQualityScore());
				prevPos++;
				qsRequired = true ;
				break;

			default:
				throw new RuntimeException("Unkown read feature operator: '"
						+ (char) feature.getOperator() + "'");
			}
		}
		if (qsRequired) {
			sb.appendCodePoint(UNSET_QUALITY_SCORE) ;
			qsRequired = false ;
		}
		if (prevPos < readLength + 1)
			sb.appendCodePoint(MATCH).append(readLength - prevPos + 1);
	}

	@Override
	public List<ReadFeature> asReadFeatureList(String string) {
		List<ReadFeature> features = new ArrayList<ReadFeature>();
		byte[] bytes = string.getBytes();
		int pos = 1;
		for (int i = 0; i < bytes.length; i++) {
			switch (bytes[i]) {
			case MATCH:
				if (bytes.length <= i + 1)
					throw new IllegalArgumentException(
							"A positive integer is expected after 'M' operator in "
									+ string + " at " + i + "th position.");
				int value = readPositiveInt(bytes, i + 1);
				pos += value;
				i += 1 + Math.log10(value);
				break;
			case DELETION:
				if (bytes.length <= i + 1)
					throw new IllegalArgumentException(
							"A positive integer is expected after 'D' operator in "
									+ string + " at " + i + "th position.");
				int len = readPositiveInt(bytes, i + 1);
				DeletionVariation dv = new DeletionVariation();
				dv.setPosition(pos);
				dv.setLength(len);
				features.add(dv);
				i += 1 + Math.log10(len);
				break;
			case INSERTION:
				if (bytes.length <= i + 1)
					throw new IllegalArgumentException(
							"A positive integer is expected after 'D' operator in "
									+ string + " at " + i + "th position.");
				byte[] seq = readNonEmptyUntil(bytes, i + 1, (byte) '.');
				InsertionVariation iv = new InsertionVariation();
				iv.setPosition(pos);
				iv.setSequence(seq);
				features.add(iv);
				pos += seq.length;
				i += seq.length + 1;
				break;
			case INSERT_A:
			case INSERT_C:
			case INSERT_G:
			case INSERT_T:
			case INSERT_N:
				InsertBase ibv = new InsertBase();
				ibv.setPosition(pos);
				ibv.setBase((byte) Character.toUpperCase(bytes[i]));
				features.add(ibv);
				pos++;
				i++;

				if (bytes[i] != 32) {
					BaseQualityScore bqs = new BaseQualityScore(pos-1, bytes[i]);
					features.add(bqs);
				}
				break;
			case BASE_A:
			case BASE_C:
			case BASE_G:
			case BASE_T:
			case BASE_N:
				ReadBase rb = new ReadBase(pos, bytes[i], bytes[i + 1]);
				features.add(rb);
				pos++;
				i++;
				break;
			case SUBSTITUTION:
				if (bytes.length <= i + 1)
					throw new IllegalArgumentException(
							"A positive integer is expected after '"
									+ (char) SUBSTITUTION + "' operator in "
									+ string + " at " + i + "th position.");

				int change = readPositiveInt(bytes, i + 1);
				SubstitutionVariation sv = new SubstitutionVariation();
				sv.setPosition(pos);
				sv.setBaseChange(new BaseChange(change));
				features.add(sv);
				pos++;
				i++;
				break;
			case SUBSTITUTION_3:
				byte score = bytes[++i];
				SubstitutionVariation svW = new SubstitutionVariation();
				svW.setPosition(pos);
				svW.setBaseChange(new BaseChange(3));
				features.add(svW);

				if (score != 32) {
					BaseQualityScore bqs = new BaseQualityScore(pos, bytes[i]);
					features.add(bqs);
				}
				pos++;
				break;
			case SUBSTITUTION_2:
				score = bytes[++i];
				SubstitutionVariation svX = new SubstitutionVariation();
				svX.setPosition(pos);
				svX.setBaseChange(new BaseChange(2));
				features.add(svX);

				if (score != 32) {
					BaseQualityScore bqs = new BaseQualityScore(pos, bytes[i]);
					features.add(bqs);
				}
				pos++;
				break;
			case SUBSTITUTION_1:
				score = bytes[++i];
				SubstitutionVariation svY = new SubstitutionVariation();
				svY.setPosition(pos);
				svY.setBaseChange(new BaseChange(1));
				features.add(svY);

				if (score != 32) {
					BaseQualityScore bqs = new BaseQualityScore(pos, bytes[i]);
					features.add(bqs);
				}
				pos++;
				break;
			case SUBSTITUTION_0:
				score = bytes[++i];
				SubstitutionVariation svZ = new SubstitutionVariation();
				svZ.setPosition(pos);
				svZ.setBaseChange(new BaseChange(0));
				features.add(svZ);

				if (score != 32) {
					BaseQualityScore bqs = new BaseQualityScore(pos, bytes[i]);
					features.add(bqs);
				}
				pos++;
				break;
			default:
				throw new IllegalArgumentException("Unkown operator: "
						+ bytes[i] + " at " + i + " in: " + string);
			}
		}
		return features;
	}

	private static final int readPositiveInt(byte[] bytes, int from) {
		int value = 0;
		int maximumAllowedValue = 1024 * 1024;
		int pos = from;
		do {
			value = value * 10 + (bytes[pos] - '0');
		} while (value < maximumAllowedValue && ++pos < bytes.length
				&& bytes[pos] >= '0' && bytes[pos] <= '9');
		return value;
	}

	private static final byte[] readNonEmptyUntil(byte[] bytes, int from,
			byte stopByte) {
		int pos = from;

		while (++pos < bytes.length && stopByte != bytes[pos])
			;
		byte[] buf = new byte[pos - from];
		System.arraycopy(bytes, from, buf, 0, buf.length);
		return buf;
	}
}

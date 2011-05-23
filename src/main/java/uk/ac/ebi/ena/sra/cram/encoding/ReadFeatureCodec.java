package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.log4j.Logger;

import uk.ac.ebi.ena.sra.cram.format.DeletionVariation;
import uk.ac.ebi.ena.sra.cram.format.InsertionVariation;
import uk.ac.ebi.ena.sra.cram.format.ReadBase;
import uk.ac.ebi.ena.sra.cram.format.ReadFeature;
import uk.ac.ebi.ena.sra.cram.format.SubstitutionVariation;
import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;
import uk.ac.ebi.ena.sra.cram.io.NullBitOutputStream;

public class ReadFeatureCodec implements BitCodec<List<ReadFeature>> {
	public BitCodec<Long> inReadPosCodec;
	public BitCodec<ReadBase> readBaseCodec;
	public BitCodec<SubstitutionVariation> substitutionCodec;
	public BitCodec<InsertionVariation> insertionCodec;
	public BitCodec<DeletionVariation> deletionCodec;

	public BitCodec<Byte> featureOperationCodec;

	private long dumpInterval = 10000;
	private long counter = 0L;
	private long inReadPosLen = 0L;
	private long nLen = 0L;
	private long sLen = 0L;
	private long iLen = 0L;
	private long dLen = 0L;

	private long sCount = 0L;
	private long iCount = 0L;
	private long dCount = 0L;

	private static Logger log = Logger.getLogger(ReadFeatureCodec.class);

	private void dump() {
		log.debug(toString());
	}

	private static final String getCodecReport(BitCodec<?> codec) {
		if (codec instanceof MeasuringCodec) {
			MeasuringCodec<?> mc = (MeasuringCodec<?>) codec;
			return mc.toString();
		}
		return null;
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder("Read features codec report: \n");
		sb.append("inReadPosCodec: ").append(getCodecReport(inReadPosCodec))
				.append("\n");
		sb.append("readBaseCodec: ").append(getCodecReport(readBaseCodec))
				.append("\n");
		sb.append("featureOperationCodec: ")
				.append(getCodecReport(featureOperationCodec)).append("\n");
		sb.append("substitutionCodec: ")
				.append(getCodecReport(substitutionCodec)).append("\n");
		sb.append("insertionCodec: ").append(getCodecReport(insertionCodec))
				.append("\n");
		sb.append("deletionCodec: ").append(getCodecReport(deletionCodec))
				.append("\n");
		return sb.toString();
	}

	private void reset() {
		counter = 0L;
		inReadPosLen = 0L;
		nLen = 0L;
		sLen = 0L;
		iLen = 0L;
		dLen = 0L;

		sCount = 0L;
		iCount = 0L;
		dCount = 0L;
	}

	@Override
	public List<ReadFeature> read(BitInputStream bis) throws IOException {
		List<ReadFeature> list = new ArrayList<ReadFeature>();
		byte op;
		int prevPos = 1;
		while ((op = featureOperationCodec.read(bis)) != ReadFeature.STOP_OPERATOR) {
			ReadFeature feature = null;
			int pos = prevPos + inReadPosCodec.read(bis).intValue();
			prevPos = pos;
			switch (op) {
			case ReadBase.operator:
				ReadBase readBase = readBaseCodec.read(bis);
				readBase.setPosition(pos);
				feature = readBase;
				break;
			case SubstitutionVariation.operator:
				SubstitutionVariation sub = substitutionCodec.read(bis);
				sub.setPosition(pos);
				feature = sub;
				break;
			case InsertionVariation.operator:
				InsertionVariation ins = insertionCodec.read(bis);
				ins.setPosition(pos);
				feature = ins;
				break;
			case DeletionVariation.operator:
				DeletionVariation del = deletionCodec.read(bis);
				del.setPosition(pos);
				feature = del;
				break;

			default:
				throw new RuntimeException("Unknown read feature operator: "
						+ (char) op);
			}
			list.add(feature);
		}

		return list;
	}

	@Override
	public long write(BitOutputStream bos, List<ReadFeature> features)
			throws IOException {

		long len = 0L;
		int prevPos = 1;
		// long tmpLen = 0 ;
		for (ReadFeature feature : features) {
			len += featureOperationCodec.write(bos, feature.getOperator());

			inReadPosLen -= len;
			len += inReadPosCodec.write(bos, (long) feature.getPosition()
					- prevPos);
			inReadPosLen += len;

			prevPos = feature.getPosition();
			switch (feature.getOperator()) {
			case ReadBase.operator:
				nLen -= len;
				len += readBaseCodec.write(bos, (ReadBase) feature);
				nLen += len;
				break;
			case SubstitutionVariation.operator:
				sCount++;
				sLen -= len;
				// tmpLen = -len ;

				len += substitutionCodec.write(bos,
						(SubstitutionVariation) feature);
				sLen += len;

				// tmpLen += len ;
				// System.out.println("sub len: " + tmpLen);

				break;
			case InsertionVariation.operator:
				iCount++;
				iLen -= len;
				len += insertionCodec.write(bos, (InsertionVariation) feature);
				iLen += len;
				break;
			case DeletionVariation.operator:
				dCount++;
				dLen -= len;
				len += deletionCodec.write(bos, (DeletionVariation) feature);
				dLen += len;
				break;

			default:
				throw new RuntimeException("Unknown read feature operator: "
						+ (char) feature.getOperator());
			}
			if (++counter >= dumpInterval) {
				dump();
				reset();
			}
		}
		featureOperationCodec.write(bos, ReadFeature.STOP_OPERATOR);

		return len;
	}

	@Override
	public long numberOfBits(List<ReadFeature> features) {
		try {
			return write(NullBitOutputStream.INSTANCE, features);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	public long getDumpInterval() {
		return dumpInterval;
	}

	public void setDumpInterval(long dumpInterval) {
		this.dumpInterval = dumpInterval;
	}

}

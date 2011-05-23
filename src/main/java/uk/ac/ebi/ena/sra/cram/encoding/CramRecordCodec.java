package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.IOException;
import java.util.List;

import org.apache.log4j.Logger;

import uk.ac.ebi.ena.sra.cram.SequenceBaseProvider;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.ReadFeature;
import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;
import uk.ac.ebi.ena.sra.cram.io.NullBitOutputStream;

public class CramRecordCodec implements BitCodec<CramRecord> {
	public BitCodec<Long> inSeqPosCodec;
	public BitCodec<Long> recordsToNextFragmentCodec;
	public BitCodec<Long> readlengthCodec;
	public BitCodec<List<ReadFeature>> variationsCodec;
	public SequenceBaseProvider sequenceBaseProvider;
	public String sequenceName;
	public long prevPosInSeq = 1L;
	public long defaultReadLength = 0L;

	private long dumpInterval = 100000;
	private long inSeqPosLen = 0L;
	private long variationsLen = 0L;
	private long readLenLen = 0L;
	private long recordCounter = 0L;

	private static Logger log = Logger.getLogger(CramRecordCodec.class);

	private void dumpLengths() {
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
		StringBuilder sb = new StringBuilder("Cram record codec report: \n");
		sb.append("inSeqPosCodec: ").append(getCodecReport(inSeqPosCodec))
				.append("\n");
		sb.append("readlengthCodec: ").append(getCodecReport(readlengthCodec))
				.append("\n");
		sb.append("recordsToNextFragmentCodec: ")
				.append(getCodecReport(recordsToNextFragmentCodec))
				.append("\n");
		sb.append("variationsCodec: ").append(getCodecReport(variationsCodec))
				.append("\n");
		return sb.toString();
	}

	@Override
	public CramRecord read(BitInputStream bis) throws IOException {
		CramRecord record = new CramRecord();

		boolean readMapped = bis.readBit();
		record.setReadMapped(readMapped);
		if (readMapped) {
			long position = prevPosInSeq + inSeqPosCodec.read(bis);
			prevPosInSeq = position;
			record.setAlignmentStart(position);

			boolean imperfectMatch = bis.readBit();
			record.setPerfectMatch(!imperfectMatch);
			if (imperfectMatch) {
				List<ReadFeature> features = variationsCodec.read(bis);
				record.setReadFeatures(features);
				// populateFeatures(record, features);
			}
		}

		record.setNegativeStrand(bis.readBit());
		record.setLastFragment(bis.readBit());
		if (!record.isLastFragment())
			record.setRecordsToNextFragment(recordsToNextFragmentCodec
					.read(bis));

		if (bis.readBit())
			record.setReadLength(readlengthCodec.read(bis));
		else
			record.setReadLength(defaultReadLength);

		return record;
	}

	// private final void populateFeatures(CramRecord record,
	// Collection<ReadFeature> features) throws IOException {
	// for (ReadFeature feature : features) {
	// switch (feature.getOperator()) {
	// case 'N':
	// if (record.getReadBases() == null)
	// record.setReadBases(new ArrayList<ReadBase>());
	// ReadBase readBase = (ReadBase) feature;
	// record.getReadBases().add(readBase);
	// break;
	// case 'S':
	// if (record.getSubstitutionVariations() == null)
	// record.setSubstitutionVariations(new ArrayList<SubstitutionVariation>());
	//
	// SubstitutionVariation v = (SubstitutionVariation) feature;
	// record.getSubstitutionVariations().add(v);
	// byte refBase = sequenceBaseProvider.getBaseAt(sequenceName,
	// record.getAlignmentStart() + v.getPosition() - 2);
	//
	// byte base = v.getBaseChange().getBaseForReference(refBase);
	//
	// v.setRefernceBase(refBase);
	// v.setBase(base);
	// break;
	// case 'I':
	// if (record.getInsertionVariations() == null)
	// record.setInsertionVariations(new ArrayList<InsertionVariation>());
	// record.getInsertionVariations().add(
	// (InsertionVariation) feature);
	// break;
	// case 'D':
	// if (record.getDeletionVariations() == null)
	// record.setDeletionVariations(new ArrayList<DeletionVariation>());
	// record.getDeletionVariations().add((DeletionVariation) feature);
	// break;
	//
	// default:
	// throw new RuntimeException("Unknown read feature operator: "
	// + (char) feature.getOperator());
	// }
	// }
	// }

	@Override
	public long write(BitOutputStream bos, CramRecord record)
			throws IOException {
		long len = 0L;
		if (record.isReadMapped()) {
			bos.write(true);
			len++;
			inSeqPosLen -= len;
			if (record.getAlignmentStart() - prevPosInSeq < 0) {
				log.error("Negative relative position in sequence: prev="
						+ prevPosInSeq);
				log.error(record.toString());
			}
			len += inSeqPosCodec.write(bos, record.getAlignmentStart()
					- prevPosInSeq);
			inSeqPosLen += len;

			prevPosInSeq = record.getAlignmentStart();
			if (!record.isPerfectMatch()) {
				bos.write(true);
				len++;
				variationsLen -= len;
				List<ReadFeature> vars = record.getReadFeatures();
				len += variationsCodec.write(bos, vars);
				variationsLen += len;
			} else
				bos.write(false);
		} else {
			throw new RuntimeException("Unmapped reads are not supported.");
		}

		bos.write(record.isNegativeStrand());
		len++;
		bos.write(record.isLastFragment());
		len++;
		if (!record.isLastFragment()) {
			len += recordsToNextFragmentCodec.write(bos,
					record.getRecordsToNextFragment());
			throw new RuntimeException(
					"Distance to next fragment not supported.");
		}

		if (record.getReadLength() != defaultReadLength) {
			bos.write(true);
			readLenLen -= len;
			len += readlengthCodec.write(bos, record.getReadLength());
			readLenLen += len;
		} else
			bos.write(false);
		len++;

		recordCounter++;
		if (recordCounter >= dumpInterval) {
			dumpLengths();
		}
		return len;
	}

	@Override
	public long numberOfBits(CramRecord record) {
		try {
			return write(NullBitOutputStream.INSTANCE, record);
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

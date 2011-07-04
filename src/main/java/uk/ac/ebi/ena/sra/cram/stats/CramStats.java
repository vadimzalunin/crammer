package uk.ac.ebi.ena.sra.cram.stats;

import java.util.Collection;
import java.util.Iterator;

import org.apache.commons.collections.bag.TreeBag;
import org.apache.commons.math.stat.Frequency;
import org.apache.log4j.Logger;

import uk.ac.ebi.ena.sra.cram.format.CramCompression;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.CramRecordBlock;
import uk.ac.ebi.ena.sra.cram.format.DeletionVariation;
import uk.ac.ebi.ena.sra.cram.format.Encoding;
import uk.ac.ebi.ena.sra.cram.format.InsertionVariation;
import uk.ac.ebi.ena.sra.cram.format.ReadBase;
import uk.ac.ebi.ena.sra.cram.format.ReadFeature;
import uk.ac.ebi.ena.sra.cram.format.SubstitutionVariation;
import uk.ac.ebi.ena.sra.cram.format.compression.CramCompressionException;
import uk.ac.ebi.ena.sra.cram.format.compression.EncodingAlgorithm;
import uk.ac.ebi.ena.sra.cram.format.compression.NumberCodecStub;

public class CramStats {
	private CramRecord prevRecord;
	private Frequency inSeqPosFreq = new Frequency();
	private Frequency inReadPosFreq = new Frequency();
	private Frequency qualityScoreFreq = new Frequency();
	private Frequency basesFreq = new Frequency();
	private Frequency readLengthFreq = new Frequency();
	private Frequency delLengthFreq = new Frequency();
	private Frequency readFeatureFreq = new Frequency();
	private Frequency distanceToNextFragmentFreq = new Frequency();

	private long nofSubstituions = 0L;
	private long nofInsertions = 0L;
	private long nofDeletions = 0L;

	private long baseCount = 0L;
	private long recordCount = 0L;

	private long firstRecordPosition = -1L;

	private TreeBag inserts = new TreeBag();

	private static Logger log = Logger.getLogger(CramStats.class);

	public void addRecord(CramRecord record) {
		recordCount++;
		baseCount += record.getReadLength();
		readLengthFreq.addValue(record.getReadLength());
		
		if (!record.isLastFragment())
			distanceToNextFragmentFreq.addValue(record
					.getRecordsToNextFragment());

		if (record.getAlignmentStart() > 0) {
			// throw new IllegalArgumentException(
			// "Records position in reference is negative: "
			// + record.getAlignmentStart());

			if (firstRecordPosition < 0)
				firstRecordPosition = record.getAlignmentStart();

			if (prevRecord == null) {
				inSeqPosFreq.addValue(0);
			} else {
				inSeqPosFreq.addValue(record.getAlignmentStart()
						- prevRecord.getAlignmentStart());
			}
			prevRecord = record;
		}

		if (!record.isReadMapped()) {
			byte[] bases = record.getReadBases();
			byte[] scores = record.getQualityScores();
			for (int i = 0; i < record.getReadLength(); i++) {
				basesFreq.addValue(bases[i]);
				qualityScoreFreq.addValue(scores[i]);
			}

			basesFreq.addValue((byte) '$');
			qualityScoreFreq.addValue((byte) -1);

			return;
		}

		Collection<ReadFeature> variationsByPosition = record.getReadFeatures();
		int prevInReadPos = 1;
		for (ReadFeature f : variationsByPosition) {
			readFeatureFreq.addValue(f.getOperator());

			inReadPosFreq.addValue(f.getPosition() - prevInReadPos);
			switch (f.getOperator()) {
			case SubstitutionVariation.operator:
				SubstitutionVariation sv = (SubstitutionVariation) f;
				if (sv.getBase() != 0)
					basesFreq.addValue(sv.getBase());
				qualityScoreFreq.addValue(sv.getQualityScore());
				nofSubstituions++;
				break;
			case ReadBase.operator:
				ReadBase rb = (ReadBase) f;
				qualityScoreFreq.addValue(rb.getQualityScore());
				break;
			case InsertionVariation.operator:
				InsertionVariation iv = (InsertionVariation) f;
				for (byte base : iv.getSequence())
					basesFreq.addValue(base);
				nofInsertions++;
				// inserts.add(new String (iv.getSequence())) ;
				break;
			case DeletionVariation.operator:
				DeletionVariation dv = (DeletionVariation) f;
				nofDeletions++;
				delLengthFreq.addValue(dv.getLength());
				break;

			default:
				break;
			}
			prevInReadPos = f.getPosition();
		}
		if (!variationsByPosition.isEmpty())
			readFeatureFreq.addValue(ReadFeature.STOP_OPERATOR);

	}

	public void adjustBlock(CramRecordBlock block)
			throws CramCompressionException {
		block.setFirstRecordPosition(firstRecordPosition);
		block.setRecordCount(recordCount);

		CramCompression compression = block.getCompression();
		if (compression == null) {
			block.setCompression(new CramCompression());
			compression = block.getCompression();
		}

		ValueFrequencyHolder holder = getValueFrequencies(basesFreq);
		compression.setBaseAlphabet(holder.values);
		compression.setBaseFrequencies(holder.frequencies);

		holder = getValueFrequencies(qualityScoreFreq);
		compression.setScoreAlphabet(holder.values);
		compression.setScoreFrequencies(holder.frequencies);

		if (readLengthFreq.getUniqueCount() == 1)
			block.setReadLength(((Long) readLengthFreq.valuesIterator().next())
					.intValue());
		else
			block.setReadLength(0);

		holder = getValueFrequencies(readLengthFreq);
		compression.setReadLengthAlphabet(holder.intValues);
		compression.setReadLengthFrequencies(holder.frequencies);

		block.getCompression().setInSeqPosEncoding(getEncoding(inSeqPosFreq));
		block.getCompression().setInReadPosEncoding(getEncoding(inReadPosFreq));
		block.getCompression().setDelLengthEncoding(getEncoding(delLengthFreq));
		block.getCompression().setReadLengthEncoding(
				getEncoding(readLengthFreq));
		block.getCompression().setRecordsToNextFragmentEncoding(
				getEncoding(distanceToNextFragmentFreq));

		holder = getValueFrequencies(readFeatureFreq);
		compression.setReadFeatureAlphabet(holder.values);
		compression.setReadFeatureFrequencies(holder.frequencies);

		log.debug("Read feature frequencies: ");
		log.debug(readFeatureFreq.toString());
		log.debug("In sequence position frequencies: ");
		log.debug(inSeqPosFreq.toString());
		log.debug("In read position frequencies: ");
		log.debug(inReadPosFreq.toString());
		log.debug("Quality score frequencies: ");
		log.debug(qualityScoreFreq.toString());
		log.debug("Insert base frequencies: ");
		log.debug(basesFreq.toString());
		log.debug("Read length frequencies: ");
		log.debug(readLengthFreq.toString());
		log.debug("Deletion length frequencies: ");
		log.debug(delLengthFreq.toString());

		log.info("Nubmer of substitutions: " + nofSubstituions);
		log.info("Nubmer of insertions: " + nofInsertions);
		log.info("Nubmer of deletions: " + nofDeletions);

		log.info("Bases: " + baseCount);
		log.info("Records: " + recordCount);

		// Iterator iterator = inserts.uniqueSet().iterator() ;
		// while (iterator.hasNext()) {
		// String insert = (String) iterator.next() ;
		// int count = inserts.getCount(insert) ;
		// System.out.printf("%d\t%s\n", count, insert);
		// }
	}

	private static final Encoding getEncoding(Frequency f)
			throws CramCompressionException {
		NumberCodecOptimiser inSeqPosOpt = new NumberCodecOptimiser();
		Iterator<Comparable<?>> valuesIterator = f.valuesIterator();
		int i = 0;
		long valueCounter = 0;
		long maxInSeqPos = Long.MIN_VALUE;
		while (valuesIterator.hasNext()) {
			Comparable<?> next = valuesIterator.next();
			Long value = (Long) next;
			if (maxInSeqPos < value)
				maxInSeqPos = value;
			i++;
			long count = f.getCount(next);
			valueCounter += count;
			inSeqPosOpt.addValue(value, count);
		}
		long binaryBitsRequiredPerPosition;
		if (maxInSeqPos < 1)
			binaryBitsRequiredPerPosition = 0;
		else
			binaryBitsRequiredPerPosition = (long) Math.ceil(Math
					.log(maxInSeqPos + 1) / Math.log(2));
		long totalBinaryBitsRequired = binaryBitsRequiredPerPosition
				* valueCounter;
		if (inSeqPosOpt.getMinLength() > totalBinaryBitsRequired) {
			return new Encoding(EncodingAlgorithm.BETA, "0,"
					+ binaryBitsRequiredPerPosition);
		} else {
			NumberCodecStub inSeqPosStub = inSeqPosOpt.getMinLengthStub();
			return new Encoding(inSeqPosStub.getEncoding(),
					inSeqPosStub.getStringRepresentation());
		}
	}

	private static class ValueFrequencyHolder {
		byte[] values;
		int[] intValues;
		int[] frequencies;

		public ValueFrequencyHolder(int size) {
			values = new byte[size];
			intValues = new int[size];
			frequencies = new int[size];
		}
	}

	private static final ValueFrequencyHolder getValueFrequencies(Frequency f) {
		ValueFrequencyHolder holder = new ValueFrequencyHolder(
				f.getUniqueCount());
		Iterator<Comparable<?>> valuesIterator = f.valuesIterator();
		int i = 0;
		while (valuesIterator.hasNext()) {
			Comparable<?> next = valuesIterator.next();
			holder.values[i] = ((Long) next).byteValue();
			holder.intValues[i] = ((Long) next).intValue();
			holder.frequencies[i] = (int) f.getCount(next);
			i++;
		}

		return holder;
	}

	public long getBaseCount() {
		return baseCount;
	}
}

package uk.ac.ebi.ena.sra.cram.stats;

import java.io.PrintStream;
import java.util.Collection;
import java.util.Iterator;
import java.util.Set;

import org.apache.commons.collections.bag.HashBag;
import org.apache.commons.math.stat.Frequency;

import uk.ac.ebi.ena.sra.compression.huffman.HuffmanCode;
import uk.ac.ebi.ena.sra.compression.huffman.HuffmanTree;
import uk.ac.ebi.ena.sra.cram.encoding.HuffmanCodec;
import uk.ac.ebi.ena.sra.cram.format.BaseQualityScore;
import uk.ac.ebi.ena.sra.cram.format.CramCompression;
import uk.ac.ebi.ena.sra.cram.format.CramHeader;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.CramRecordBlock;
import uk.ac.ebi.ena.sra.cram.format.DeletionVariation;
import uk.ac.ebi.ena.sra.cram.format.Encoding;
import uk.ac.ebi.ena.sra.cram.format.InsertBase;
import uk.ac.ebi.ena.sra.cram.format.InsertionVariation;
import uk.ac.ebi.ena.sra.cram.format.ReadAnnotation;
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
	private Frequency readGroupIndexFreq = new Frequency();
	private Frequency mappingQualityFreq = new Frequency();

	private long nofSubstituions = 0L;
	private long nofInsertions = 0L;
	private long nofInsertedBases = 0L;
	private long nofDeletions = 0L;

	private long baseCount = 0L;
	private long recordCount = 0L;

	private long firstRecordPosition = -1L;

	private HashBag readAnnoKeyBag = new HashBag();

	private final CramHeader header;
	private final PrintStream statsPS;
	private final float qualityBudget;

	public CramStats(CramHeader header, PrintStream statsPS, float qualityBudget) {
		this.header = header;
		this.statsPS = statsPS;
		this.qualityBudget = qualityBudget;
	}

	public CramStats(CramHeader header, PrintStream statsPS) {
		this(header, statsPS, 0);
	}

	public void addRecord(CramRecord record) {
		recordCount++;
		baseCount += record.getReadLength();
		readLengthFreq.addValue(record.getReadLength());

		if (record.getAnnotations() != null) {
			for (ReadAnnotation a : record.getAnnotations())
				readAnnoKeyBag.add(a);
		}

		readGroupIndexFreq.addValue(record.getReadGroupID());

		if (!record.isLastFragment())
			distanceToNextFragmentFreq.addValue(record
					.getRecordsToNextFragment());

		if (record.getAlignmentStart() > 0) {
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
				if (bases != null)
					basesFreq.addValue(bases[i]);
				if (scores != null && scores.length != 0)
					qualityScoreFreq.addValue(scores[i]);
			}

			basesFreq.addValue((byte) '$');
			qualityScoreFreq.addValue((byte) -1);

			return;
		} else
			mappingQualityFreq.addValue(record.getMappingQuality());

		Collection<ReadFeature> variationsByPosition = record.getReadFeatures();
		if (variationsByPosition != null) {
			int prevInReadPos = 0;
			for (ReadFeature f : variationsByPosition) {
				readFeatureFreq.addValue(f.getOperator());

				inReadPosFreq.addValue(f.getPosition() - prevInReadPos);
				switch (f.getOperator()) {
				case SubstitutionVariation.operator:
					SubstitutionVariation sv = (SubstitutionVariation) f;
					if (sv.getBase() != 0)
						basesFreq.addValue(sv.getBase());
					nofSubstituions++;
					break;
				case ReadBase.operator:
					ReadBase rb = (ReadBase) f;
					basesFreq.addValue(rb.getBase());
					qualityScoreFreq.addValue(rb.getQualityScore());
					break;
				case InsertionVariation.operator:
					InsertionVariation iv = (InsertionVariation) f;
					for (byte base : iv.getSequence())
						basesFreq.addValue(base);
					nofInsertions++;
					nofInsertedBases += iv.getSequence().length;
					basesFreq.addValue((byte) '$');
					// inserts.add(new String (iv.getSequence())) ;
					break;
				case DeletionVariation.operator:
					DeletionVariation dv = (DeletionVariation) f;
					nofDeletions++;
					delLengthFreq.addValue(dv.getLength());
					break;
				case InsertBase.operator:
					InsertBase ib = (InsertBase) f;
					basesFreq.addValue(ib.getBase());
					nofInsertions++;
					nofInsertedBases++;
					break;
				case BaseQualityScore.operator:
					BaseQualityScore bqs = (BaseQualityScore) f;
					qualityScoreFreq.addValue(bqs.getQualityScore());
					break;

				default:
					break;
				}
				prevInReadPos = f.getPosition();
			}
			if (!variationsByPosition.isEmpty())
				readFeatureFreq.addValue(ReadFeature.STOP_OPERATOR);
		}
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

		holder = getValueFrequencies(mappingQualityFreq);
		compression.setMappingQualityAlphabet(holder.values);
		compression.setMappingQualityFrequencies(holder.frequencies);

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

		if (header != null && header.getReadAnnotations() != null) {
			Set uniqueSet = readAnnoKeyBag.uniqueSet();
			int[] raIndexes = new int[uniqueSet.size()];
			int[] raFreqs = new int[raIndexes.length];
			int i = 0;
			for (Object o : uniqueSet) {
				ReadAnnotation ra = (ReadAnnotation) o;
				int count = readAnnoKeyBag.getCount(ra);
				if (count == 0)
					continue;

				int index = header.getReadAnnotations().indexOf(ra);
				if (index < 0)
					throw new CramCompressionException(
							"Annotation not found in the dictionary: "
									+ ra.getKey());

				raIndexes[i] = index;
				raFreqs[i] = count;
				i++;
			}
			compression.setReadAnnotationIndexes(raIndexes);
			compression.setReadAnnotationFrequencies(raFreqs);
		} else {
			compression.setReadAnnotationIndexes(new int[0]);
			compression.setReadAnnotationFrequencies(new int[0]);
		}

		holder = getValueFrequencies(readGroupIndexFreq);
		compression.setReadGroupIndexes(holder.intValues);
		compression.setReadGroupFrequencies(holder.frequencies);

		if (statsPS != null) {
			statsPS.println("Read feature frequencies: ");
			statsPS.println(readFeatureFreq.toString());
			statsPS.println("Base frequencies: ");
			statsPS.println(basesFreq.toString());
			statsPS.println("Quality score frequencies: ");
			statsPS.println(qualityScoreFreq.toString());
			statsPS.println("In read position frequencies: ");
			statsPS.println(inReadPosFreq.toString());
			statsPS.println("Read length frequencies: ");
			statsPS.println(readLengthFreq.toString());
			statsPS.println("Deletion length frequencies: ");
			statsPS.println(delLengthFreq.toString());
			statsPS.println("In sequence position frequencies: ");
			statsPS.println(inSeqPosFreq.toString());
			statsPS.println("Records to next fragment frequencies: ");
			statsPS.println(distanceToNextFragmentFreq.toString());
			statsPS.println("Mapping quality frequencies: ");
			statsPS.println(mappingQualityFreq.toString());

			statsPS.println("Nubmer of substitutions: " + nofSubstituions);
			statsPS.println("Nubmer of insertions: " + nofInsertions);
			statsPS.println("Nubmer of inserted bases: " + nofInsertedBases);
			statsPS.println("Nubmer of deletions: " + nofDeletions);

			statsPS.println("Bases: " + baseCount);
			statsPS.println("Records: " + recordCount);
		}
	}

	private static final Encoding getEncoding(Frequency f)
			throws CramCompressionException {
		NumberCodecOptimiser optimiser = new NumberCodecOptimiser();
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
			optimiser.addValue(value, count);
		}

		// long huffmanBits = getHuffmanBitLengthEstimate(f);
		// if (optimiser.getMinLength() > huffmanBits) {
		// System.out.printf("Huffman is better: %d vs %d\n",
		// optimiser.getMinLength(), huffmanBits);
		// } else
		// System.out.printf("Huffman is worse: %d vs %d\n",
		// optimiser.getMinLength(), huffmanBits);

		long binaryBitsRequiredPerPosition;
		if (maxInSeqPos < 1)
			binaryBitsRequiredPerPosition = 0;
		else
			binaryBitsRequiredPerPosition = (long) Math.ceil(Math
					.log(maxInSeqPos + 1) / Math.log(2));
		long totalBinaryBitsRequired = binaryBitsRequiredPerPosition
				* valueCounter;

		if (optimiser.getMinLength() > totalBinaryBitsRequired) {
			return new Encoding(EncodingAlgorithm.BETA, "0,"
					+ binaryBitsRequiredPerPosition);
		} else {
			NumberCodecStub inSeqPosStub = optimiser.getMinLengthStub();
			return new Encoding(inSeqPosStub.getEncoding(),
					inSeqPosStub.getStringRepresentation());
		}
	}

	private static long getHuffmanBitLengthEstimate(Frequency f) {
		Iterator<Comparable<?>> valuesIterator = f.valuesIterator();
		int i = 0;
		Long[] values = new Long[f.getUniqueCount()];
		int[] freqs = new int[values.length];
		while (valuesIterator.hasNext()) {
			Comparable<?> next = valuesIterator.next();
			values[i] = ((Long) next).longValue();
			freqs[i] = (int) f.getCount(next);
			i++;
		}
		HuffmanTree<Long> tree = HuffmanCode.buildTree(freqs, values);

		valuesIterator = f.valuesIterator();
		HuffmanCodec<Long> codec = new HuffmanCodec<Long>(tree);
		long totalBits = 0;
		while (valuesIterator.hasNext()) {
			Comparable<?> next = valuesIterator.next();
			long bits = codec.numberOfBits(((Long) next).longValue());
			totalBits += bits * (int) f.getCount(next);
		}

		// length of int array plus length of int array in bits:
		int treeLen = ((4 + values.length * 4) + (4 + freqs.length * 4)) * 8;
		return totalBits + treeLen;
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

package uk.ac.ebi.ena.sra.cram.impl;

import java.util.List;

import uk.ac.ebi.ena.sra.compression.huffman.HuffmanCode;
import uk.ac.ebi.ena.sra.compression.huffman.HuffmanTree;
import uk.ac.ebi.ena.sra.cram.SequenceBaseProvider;
import uk.ac.ebi.ena.sra.cram.Utils;
import uk.ac.ebi.ena.sra.cram.encoding.BaseChange;
import uk.ac.ebi.ena.sra.cram.encoding.BaseChangeCodec;
import uk.ac.ebi.ena.sra.cram.encoding.BaseSequenceCodec;
import uk.ac.ebi.ena.sra.cram.encoding.BitCodec;
import uk.ac.ebi.ena.sra.cram.encoding.CramRecordCodec;
import uk.ac.ebi.ena.sra.cram.encoding.DeletionVariationCodec;
import uk.ac.ebi.ena.sra.cram.encoding.HuffmanCodec;
import uk.ac.ebi.ena.sra.cram.encoding.InsertionVariationCodec;
import uk.ac.ebi.ena.sra.cram.encoding.MeasuringCodec;
import uk.ac.ebi.ena.sra.cram.encoding.ReadBaseCodec;
import uk.ac.ebi.ena.sra.cram.encoding.ReadFeatureCodec;
import uk.ac.ebi.ena.sra.cram.encoding.SubstitutionVariationCodec;
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
import uk.ac.ebi.ena.sra.cram.format.compression.NumberCodecFactory;
import uk.ac.ebi.ena.sra.cram.format.compression.NumberCodecStub;

public class RecordCodecFactory {
	private long dumpRecordInterval = Long.MAX_VALUE;
	private long dumpReadFeaturesInterval = Long.MAX_VALUE;

	public BitCodec<CramRecord> createRecordCodec(CramRecordBlock block,
			SequenceBaseProvider referenceBaseProvider)
			throws CramCompressionException {
		CramCompression compression = block.getCompression();
		// given a bunch of block info and compression info create a codec:
		CramRecordCodec recordCodec = new CramRecordCodec();
		recordCodec.setDumpInterval(dumpRecordInterval);
		recordCodec.prevPosInSeq = block.getFirstRecordPosition();
		recordCodec.inSeqPosCodec = new MeasuringCodec<Long>(
				createStub(compression.getInSeqPosEncoding()), "Refpos codec");

		// hack:
		Long[] readLengthAlphabet = new Long[compression
				.getReadLengthAlphabet().length];
		for (int i = 0; i < readLengthAlphabet.length; i++)
			readLengthAlphabet[i] = new Long(
					compression.getReadLengthAlphabet()[i]);
		HuffmanTree<Long> readlengthTree = HuffmanCode.buildTree(
				compression.getReadLengthFrequencies(), readLengthAlphabet);
		HuffmanCodec<Long> readlengthCodec = new HuffmanCodec<Long>(
				readlengthTree);
		// createStub(compression.getReadLengthEncoding());
		recordCodec.readlengthCodec = new MeasuringCodec<Long>(readlengthCodec,
				"Read length codec");
		recordCodec.recordsToNextFragmentCodec = new MeasuringCodec<Long>(
				createStub(compression.getRecordsToNextFragmentEncoding()),
				"Rank codec");
		recordCodec.sequenceBaseProvider = referenceBaseProvider;
		recordCodec.defaultReadLength = block.getReadLength();

		NumberCodecStub inReadPosCodecStub = createStub(compression
				.getInReadPosEncoding());
		HuffmanTree<Byte> qualityScoreTree = HuffmanCode.buildTree(
				compression.getScoreFrequencies(),
				Utils.autobox(compression.getScoreAlphabet()));
		BitCodec<Byte> qualityScoreCodec = new MeasuringCodec<Byte>(
				new HuffmanCodec<Byte>(qualityScoreTree), "Quality score codec");

		ReadFeatureCodec readFearureCodec = new ReadFeatureCodec();
		readFearureCodec.setDumpInterval(dumpReadFeaturesInterval);
		readFearureCodec.inReadPosCodec = new MeasuringCodec<Long>(
				inReadPosCodecStub, "Position in read");
		HuffmanTree<Byte> featureOperatorTree = HuffmanCode.buildTree(
				compression.getReadFeatureFrequencies(),
				Utils.autobox(compression.getReadFeatureAlphabet()));
		HuffmanCodec<Byte> featureOperatorCodec = new HuffmanCodec<Byte>(
				featureOperatorTree);
		readFearureCodec.featureOperationCodec = new MeasuringCodec<Byte>(
				featureOperatorCodec, "Read feature operators");

		ReadBaseCodec readBaseCodec = new ReadBaseCodec();
		readBaseCodec.qualityScoreCodec = qualityScoreCodec;
		readFearureCodec.readBaseCodec = new MeasuringCodec<ReadBase>(
				readBaseCodec, "Read base codec");

		DeletionVariationCodec deletionCodec = new DeletionVariationCodec();
		deletionCodec.dellengthPosCodec = new MeasuringCodec<Long>(
				createStub(compression.getDelLengthEncoding()),
				"Deletion length codec");
		readFearureCodec.deletionCodec = new MeasuringCodec<DeletionVariation>(
				deletionCodec, "Deletion codec");

		SubstitutionVariationCodec substitutionCodec = new SubstitutionVariationCodec();
		substitutionCodec.baseChangeCodec = new MeasuringCodec<BaseChange>(
				new BaseChangeCodec(), "Base change codec");
		substitutionCodec.qualityScoreCodec = qualityScoreCodec;
		readFearureCodec.substitutionCodec = new MeasuringCodec<SubstitutionVariation>(
				substitutionCodec, "Substitution codec");

		InsertionVariationCodec insertionCodec = new InsertionVariationCodec();
		insertionCodec.insertBasesCodec = new BaseSequenceCodec(
				BaseSequenceCodec.BaseCodecType.RAISED, "ACGTSN".getBytes());
		readFearureCodec.insertionCodec = new MeasuringCodec<InsertionVariation>(
				insertionCodec, "Insertion codec");

		recordCodec.variationsCodec = new MeasuringCodec<List<ReadFeature>>(
				readFearureCodec, "Variations codec");

		return recordCodec;
	}

	private static NumberCodecStub createStub(Encoding encoding)
			throws CramCompressionException {
		NumberCodecStub stub = NumberCodecFactory.createStub(encoding
				.getAlgorithm());
		stub.initFromString(encoding.getParameters());
		return stub;
	}

	public long getDumpRecordInterval() {
		return dumpRecordInterval;
	}

	public void setDumpRecordInterval(long dumpRecordInterval) {
		this.dumpRecordInterval = dumpRecordInterval;
	}

	public long getDumpReadFeaturesInterval() {
		return dumpReadFeaturesInterval;
	}

	public void setDumpReadFeaturesInterval(long dumpReadFeaturesInterval) {
		this.dumpReadFeaturesInterval = dumpReadFeaturesInterval;
	}
}

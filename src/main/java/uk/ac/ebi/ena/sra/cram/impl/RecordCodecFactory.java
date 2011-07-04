package uk.ac.ebi.ena.sra.cram.impl;

import java.util.Enumeration;
import java.util.List;

import javax.swing.tree.DefaultMutableTreeNode;

import uk.ac.ebi.ena.sra.compression.huffman.HuffmanCode;
import uk.ac.ebi.ena.sra.compression.huffman.HuffmanTree;
import uk.ac.ebi.ena.sra.cram.SequenceBaseProvider;
import uk.ac.ebi.ena.sra.cram.Utils;
import uk.ac.ebi.ena.sra.cram.encoding.BaseChange;
import uk.ac.ebi.ena.sra.cram.encoding.BaseChangeCodec;
import uk.ac.ebi.ena.sra.cram.encoding.BaseSequenceCodec;
import uk.ac.ebi.ena.sra.cram.encoding.BitCodec;
import uk.ac.ebi.ena.sra.cram.encoding.ByteArrayHuffmanCodec;
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

	public DefaultMutableTreeNode buildCodecTree(CramRecordBlock block,
			SequenceBaseProvider referenceBaseProvider)
			throws CramCompressionException {
		DefaultMutableTreeNode root = new DefaultMutableTreeNode();

		CramRecordCodec recordCodec = new CramRecordCodec();
		root.setUserObject(new MeasuringCodec<CramRecord>(recordCodec,
				"Cram record codec"));

		CramCompression compression = block.getCompression();
		// given a bunch of block info and compression info create a codec:
		recordCodec.prevPosInSeq = block.getFirstRecordPosition();
		MeasuringCodec<Long> refPosMeasuringCodec = new MeasuringCodec<Long>(
				createStub(compression.getInSeqPosEncoding()), "Refpos codec");
		recordCodec.inSeqPosCodec = refPosMeasuringCodec;
		root.add(new DefaultMutableTreeNode(refPosMeasuringCodec));

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
		MeasuringCodec<Long> readlengthMeasuringCodec = new MeasuringCodec<Long>(
				readlengthCodec, "Read length codec");
		recordCodec.readlengthCodec = readlengthMeasuringCodec;
		root.add(new DefaultMutableTreeNode(readlengthMeasuringCodec));

		MeasuringCodec<Long> rankCodec = new MeasuringCodec<Long>(
				createStub(compression.getRecordsToNextFragmentEncoding()),
				"Rank codec");
		recordCodec.recordsToNextFragmentCodec = rankCodec;
		root.add(new DefaultMutableTreeNode(rankCodec));

		recordCodec.sequenceBaseProvider = referenceBaseProvider;
		recordCodec.defaultReadLength = block.getReadLength();

		// read feature list node:

		DefaultMutableTreeNode rflNode = new DefaultMutableTreeNode();
		root.add(rflNode);

		NumberCodecStub inReadPosCodecStub = createStub(compression
				.getInReadPosEncoding());

		HuffmanTree<Byte> qualityScoreTree = HuffmanCode.buildTree(
				compression.getScoreFrequencies(),
				Utils.autobox(compression.getScoreAlphabet()));
		MeasuringCodec<Byte> qsMeasuringCodec = new MeasuringCodec<Byte>(
				new HuffmanCodec<Byte>(qualityScoreTree), "Quality score codec");
		BitCodec<Byte> qualityScoreCodec = qsMeasuringCodec;
		rflNode.add(new DefaultMutableTreeNode(qsMeasuringCodec));

		ReadFeatureCodec readFearureCodec = new ReadFeatureCodec();
		MeasuringCodec<Long> posInReadMeasuringCodec = new MeasuringCodec<Long>(
				inReadPosCodecStub, "Position in read");
		readFearureCodec.inReadPosCodec = posInReadMeasuringCodec;
		rflNode.add(new DefaultMutableTreeNode(posInReadMeasuringCodec));

		HuffmanTree<Byte> featureOperatorTree = HuffmanCode.buildTree(
				compression.getReadFeatureFrequencies(),
				Utils.autobox(compression.getReadFeatureAlphabet()));
		HuffmanCodec<Byte> featureOperatorCodec = new HuffmanCodec<Byte>(
				featureOperatorTree);
		MeasuringCodec<Byte> rfopsMeasuringCodec = new MeasuringCodec<Byte>(
				featureOperatorCodec, "Read feature operators");
		readFearureCodec.featureOperationCodec = rfopsMeasuringCodec;
		rflNode.add(new DefaultMutableTreeNode(rfopsMeasuringCodec));

		// read base node:
		DefaultMutableTreeNode rbNode = new DefaultMutableTreeNode();
		rflNode.add(rbNode);

		ReadBaseCodec readBaseCodec = new ReadBaseCodec();
		MeasuringCodec<Byte> rbqsMeasuringCodec = new MeasuringCodec<Byte>(
				qualityScoreCodec, "RBQS codec");
		readBaseCodec.qualityScoreCodec = rbqsMeasuringCodec;
		rbNode.add(new DefaultMutableTreeNode(rbqsMeasuringCodec));
		MeasuringCodec<ReadBase> readBaseMeasuringCodec = new MeasuringCodec<ReadBase>(
				readBaseCodec, "Read base codec");
		readFearureCodec.readBaseCodec = readBaseMeasuringCodec;
		rbNode.setUserObject(readBaseMeasuringCodec);

		// deletion node:
		DefaultMutableTreeNode delNode = new DefaultMutableTreeNode();
		rflNode.add(delNode);

		DeletionVariationCodec deletionCodec = new DeletionVariationCodec();
		MeasuringCodec<Long> delLenMeasuringCodec = new MeasuringCodec<Long>(
				createStub(compression.getDelLengthEncoding()),
				"Deletion length codec");
		deletionCodec.dellengthPosCodec = delLenMeasuringCodec;
		delNode.add(new DefaultMutableTreeNode(delLenMeasuringCodec));
		MeasuringCodec<DeletionVariation> delMeasuringCodec = new MeasuringCodec<DeletionVariation>(
				deletionCodec, "Deletion codec");
		readFearureCodec.deletionCodec = delMeasuringCodec;
		delNode.setUserObject(delMeasuringCodec);

		// substitution node:
		DefaultMutableTreeNode subNode = new DefaultMutableTreeNode();
		rflNode.add(subNode);

		SubstitutionVariationCodec substitutionCodec = new SubstitutionVariationCodec();
		MeasuringCodec<BaseChange> baseChangeMeasuringCodec = new MeasuringCodec<BaseChange>(
				new BaseChangeCodec(), "Base change codec");
		substitutionCodec.baseChangeCodec = baseChangeMeasuringCodec;
		subNode.add(new DefaultMutableTreeNode(baseChangeMeasuringCodec));

		MeasuringCodec<Byte> subqsMeasuringCodec = new MeasuringCodec<Byte>(
				qualityScoreCodec, "Sub QS codec");
		substitutionCodec.qualityScoreCodec = subqsMeasuringCodec;
		subNode.add(new DefaultMutableTreeNode(subqsMeasuringCodec));
		MeasuringCodec<SubstitutionVariation> subMeasuringCodec = new MeasuringCodec<SubstitutionVariation>(
				substitutionCodec, "Substitution codec");
		readFearureCodec.substitutionCodec = subMeasuringCodec;
		subNode.setUserObject(subMeasuringCodec);

		// insertion node:
		DefaultMutableTreeNode insNode = new DefaultMutableTreeNode();
		rflNode.add(insNode);

		InsertionVariationCodec insertionCodec = new InsertionVariationCodec();
		insertionCodec.insertBasesCodec = new BaseSequenceCodec(
				BaseSequenceCodec.BaseCodecType.RAISED, "ACGTSN".getBytes());
		MeasuringCodec<InsertionVariation> insMeasuringCodec = new MeasuringCodec<InsertionVariation>(
				insertionCodec, "Insertion codec");
		readFearureCodec.insertionCodec = insMeasuringCodec;
		insNode.setUserObject(insMeasuringCodec);

		MeasuringCodec<List<ReadFeature>> varMeasuringCodec = new MeasuringCodec<List<ReadFeature>>(
				readFearureCodec, "Variations codec");
		recordCodec.variationsCodec = varMeasuringCodec;
		rflNode.setUserObject(varMeasuringCodec);

		// unampped bases:
		ByteArrayHuffmanCodec unmappedBasesCodec = new ByteArrayHuffmanCodec(
				compression.getBaseAlphabet(),
				compression.getBaseFrequencies(), (byte) '$');
		MeasuringCodec<byte[]> unmappedBasesMeasuringCodec = new MeasuringCodec<byte[]>(
				unmappedBasesCodec, "Unmapped bases codec");
		recordCodec.basesCodec = unmappedBasesMeasuringCodec;
		root.add(new DefaultMutableTreeNode(unmappedBasesMeasuringCodec));

		// unampped quality scores:
		ByteArrayHuffmanCodec unmappedScoresCodec = new ByteArrayHuffmanCodec(
				compression.getScoreAlphabet(),
				compression.getScoreFrequencies(), (byte) -1);
		MeasuringCodec<byte[]> unmappedScoresMeasuringCodec = new MeasuringCodec<byte[]>(
				unmappedScoresCodec, "Unmapped quality scores codec");
		recordCodec.qualitiesCodec = unmappedScoresMeasuringCodec;
		root.add(new DefaultMutableTreeNode(unmappedScoresMeasuringCodec));

		return root;
	}

	public void dump(DefaultMutableTreeNode node) {
		MeasuringCodec<?> codec = (MeasuringCodec<?>) node.getUserObject();
		for (int i = 0; i < node.getLevel(); i++)
			System.err.print("\t");

		System.err.println(codec.toString());
		if (!node.isLeaf()) {
			Enumeration children = node.children();
			while (children.hasMoreElements()) {
				DefaultMutableTreeNode childNode = (DefaultMutableTreeNode) children
						.nextElement();
				dump(childNode);
			}
		}
	}

	public BitCodec<CramRecord> createRecordCodec(CramRecordBlock block,
			SequenceBaseProvider referenceBaseProvider)
			throws CramCompressionException {
		CramCompression compression = block.getCompression();
		// given a bunch of block info and compression info create a codec:
		CramRecordCodec recordCodec = new CramRecordCodec();
		recordCodec.prevPosInSeq = block.getFirstRecordPosition();
		recordCodec.inSeqPosCodec = createStub(compression
				.getInSeqPosEncoding());

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
		recordCodec.readlengthCodec = readlengthCodec;
		recordCodec.recordsToNextFragmentCodec = createStub(compression
				.getRecordsToNextFragmentEncoding());
		recordCodec.sequenceBaseProvider = referenceBaseProvider;
		recordCodec.defaultReadLength = block.getReadLength();

		NumberCodecStub inReadPosCodecStub = createStub(compression
				.getInReadPosEncoding());
		HuffmanTree<Byte> qualityScoreTree = HuffmanCode.buildTree(
				compression.getScoreFrequencies(),
				Utils.autobox(compression.getScoreAlphabet()));
		BitCodec<Byte> qualityScoreCodec = new HuffmanCodec<Byte>(
				qualityScoreTree);

		ReadFeatureCodec readFearureCodec = new ReadFeatureCodec();
		readFearureCodec.inReadPosCodec = inReadPosCodecStub;
		HuffmanTree<Byte> featureOperatorTree = HuffmanCode.buildTree(
				compression.getReadFeatureFrequencies(),
				Utils.autobox(compression.getReadFeatureAlphabet()));
		HuffmanCodec<Byte> featureOperatorCodec = new HuffmanCodec<Byte>(
				featureOperatorTree);
		readFearureCodec.featureOperationCodec = featureOperatorCodec;

		ReadBaseCodec readBaseCodec = new ReadBaseCodec();
		readBaseCodec.qualityScoreCodec = qualityScoreCodec;
		readFearureCodec.readBaseCodec = readBaseCodec;

		DeletionVariationCodec deletionCodec = new DeletionVariationCodec();
		deletionCodec.dellengthPosCodec = createStub(compression
				.getDelLengthEncoding());
		readFearureCodec.deletionCodec = deletionCodec;

		SubstitutionVariationCodec substitutionCodec = new SubstitutionVariationCodec();
		substitutionCodec.baseChangeCodec = new BaseChangeCodec();
		substitutionCodec.qualityScoreCodec = qualityScoreCodec;
		readFearureCodec.substitutionCodec = substitutionCodec;

		InsertionVariationCodec insertionCodec = new InsertionVariationCodec();
		insertionCodec.insertBasesCodec = new BaseSequenceCodec(
				BaseSequenceCodec.BaseCodecType.RAISED, "ACGTSN".getBytes());
		readFearureCodec.insertionCodec = insertionCodec;

		recordCodec.variationsCodec = readFearureCodec;

		// unampped bases:
		ByteArrayHuffmanCodec unmappedBasesCodec = new ByteArrayHuffmanCodec(
				compression.getBaseAlphabet(),
				compression.getBaseFrequencies(), (byte) '$');
		MeasuringCodec<byte[]> unmappedBasesMeasuringCodec = new MeasuringCodec<byte[]>(
				unmappedBasesCodec, "Unmapped bases codec");
		recordCodec.basesCodec = unmappedBasesMeasuringCodec;

		// unampped quality scores:
		ByteArrayHuffmanCodec unmappedScoresCodec = new ByteArrayHuffmanCodec(
				compression.getScoreAlphabet(),
				compression.getScoreFrequencies(), (byte) -1);
		MeasuringCodec<byte[]> unmappedScoresMeasuringCodec = new MeasuringCodec<byte[]>(
				unmappedScoresCodec, "Unmapped quality scores codec");
		recordCodec.qualitiesCodec = unmappedScoresMeasuringCodec;

		return recordCodec;
	}

	private static NumberCodecStub createStub(Encoding encoding)
			throws CramCompressionException {
		NumberCodecStub stub = NumberCodecFactory.createStub(encoding
				.getAlgorithm());
		stub.initFromString(encoding.getParameters());
		return stub;
	}

}

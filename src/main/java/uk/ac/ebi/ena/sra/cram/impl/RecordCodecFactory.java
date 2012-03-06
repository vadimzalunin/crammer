package uk.ac.ebi.ena.sra.cram.impl;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.Enumeration;
import java.util.List;
import java.util.TreeMap;

import javax.swing.tree.DefaultMutableTreeNode;

import uk.ac.ebi.ena.sra.compression.huffman.HuffmanCode;
import uk.ac.ebi.ena.sra.compression.huffman.HuffmanTree;
import uk.ac.ebi.ena.sra.cram.SequenceBaseProvider;
import uk.ac.ebi.ena.sra.cram.Utils;
import uk.ac.ebi.ena.sra.cram.encoding.BaseChange;
import uk.ac.ebi.ena.sra.cram.encoding.BaseChangeCodec;
import uk.ac.ebi.ena.sra.cram.encoding.BitCodec;
import uk.ac.ebi.ena.sra.cram.encoding.ByteArrayHuffmanCodec;
import uk.ac.ebi.ena.sra.cram.encoding.CramRecordCodec;
import uk.ac.ebi.ena.sra.cram.encoding.DeletionVariationCodec;
import uk.ac.ebi.ena.sra.cram.encoding.FixedLengthByteArrayHuffmanCodec;
import uk.ac.ebi.ena.sra.cram.encoding.HuffmanByteCodec;
import uk.ac.ebi.ena.sra.cram.encoding.HuffmanCodec;
import uk.ac.ebi.ena.sra.cram.encoding.InsertionVariationCodec;
import uk.ac.ebi.ena.sra.cram.encoding.MeasuringCodec;
import uk.ac.ebi.ena.sra.cram.encoding.ReadAnnotationCodec;
import uk.ac.ebi.ena.sra.cram.encoding.ReadBaseCodec;
import uk.ac.ebi.ena.sra.cram.encoding.ReadFeatureCodec;
import uk.ac.ebi.ena.sra.cram.encoding.SingleValueBitCodec;
import uk.ac.ebi.ena.sra.cram.encoding.SubstitutionVariationCodec;
import uk.ac.ebi.ena.sra.cram.encoding.VariableLengthByteArrayHuffmanCodec;
import uk.ac.ebi.ena.sra.cram.format.BaseQualityScore;
import uk.ac.ebi.ena.sra.cram.format.ByteFrequencies;
import uk.ac.ebi.ena.sra.cram.format.CramCompression;
import uk.ac.ebi.ena.sra.cram.format.CramHeader;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.CramRecordBlock;
import uk.ac.ebi.ena.sra.cram.format.DeletionVariation;
import uk.ac.ebi.ena.sra.cram.format.Encoding;
import uk.ac.ebi.ena.sra.cram.format.InsertBase;
import uk.ac.ebi.ena.sra.cram.format.InsertionVariation;
import uk.ac.ebi.ena.sra.cram.format.IntFrequencies;
import uk.ac.ebi.ena.sra.cram.format.ReadAnnotation;
import uk.ac.ebi.ena.sra.cram.format.ReadBase;
import uk.ac.ebi.ena.sra.cram.format.ReadFeature;
import uk.ac.ebi.ena.sra.cram.format.SubstitutionVariation;
import uk.ac.ebi.ena.sra.cram.format.compression.CramCompressionException;
import uk.ac.ebi.ena.sra.cram.format.compression.EncodingAlgorithm;
import uk.ac.ebi.ena.sra.cram.format.compression.NumberCodecFactory;
import uk.ac.ebi.ena.sra.cram.format.compression.NumberCodecStub;
import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;

public class RecordCodecFactory {

	public DefaultMutableTreeNode buildCodecTree(CramHeader header, CramRecordBlock block,
			SequenceBaseProvider referenceBaseProvider) throws CramCompressionException {
		DefaultMutableTreeNode root = new DefaultMutableTreeNode();

		CramRecordCodec recordCodec = new CramRecordCodec();
		root.setUserObject(new MeasuringCodec<CramRecord>(recordCodec, "Cram record codec"));

		CramCompression compression = block.getCompression();
		// given a bunch of block info and compression info create a codec:
		recordCodec.prevPosInSeq = block.getFirstRecordPosition();
		MeasuringCodec<Long> refPosMeasuringCodec = new MeasuringCodec<Long>(
				createStub(compression.getInSeqPosEncoding()), "Refpos codec");
		recordCodec.inSeqPosCodec = refPosMeasuringCodec;
		root.add(new DefaultMutableTreeNode(refPosMeasuringCodec));

		// hack:
		Long[] readLengthAlphabet = new Long[compression.getReadLengthAlphabet().length];
		for (int i = 0; i < readLengthAlphabet.length; i++)
			readLengthAlphabet[i] = new Long(compression.getReadLengthAlphabet()[i]);
		HuffmanTree<Long> readlengthTree = HuffmanCode.buildTree(compression.getReadLengthFrequencies(),
				readLengthAlphabet);
		HuffmanCodec<Long> readlengthCodec = new HuffmanCodec<Long>(readlengthTree);
		// createStub(compression.getReadLengthEncoding());
		MeasuringCodec<Long> readlengthMeasuringCodec = new MeasuringCodec<Long>(readlengthCodec, "Read length codec");
		recordCodec.readlengthCodec = readlengthMeasuringCodec;
		root.add(new DefaultMutableTreeNode(readlengthMeasuringCodec));

		MeasuringCodec<Long> rankCodec = new MeasuringCodec<Long>(
				createStub(compression.getRecordsToNextFragmentEncoding()), "Rank codec");
		recordCodec.recordsToNextFragmentCodec = rankCodec;
		root.add(new DefaultMutableTreeNode(rankCodec));

		recordCodec.sequenceBaseProvider = referenceBaseProvider;
		recordCodec.defaultReadLength = block.getReadLength();

		// read feature list node:

		DefaultMutableTreeNode rflNode = new DefaultMutableTreeNode();
		root.add(rflNode);

		NumberCodecStub inReadPosCodecStub = createStub(compression.getInReadPosEncoding());

		HuffmanTree<Byte> baseTree = HuffmanCode.buildTree(compression.getBaseFrequencies(),
				Utils.autobox(compression.getBaseAlphabet()));
		MeasuringCodec<Byte> baseMeasuringCodec = new MeasuringCodec<Byte>(new HuffmanByteCodec(baseTree), "Base codec");
		final BitCodec<Byte> baseCodec = baseMeasuringCodec;
		rflNode.add(new DefaultMutableTreeNode(baseCodec));

		HuffmanTree<Byte> qualityScoreTree = HuffmanCode.buildTree(compression.getScoreFrequencies(),
				Utils.autobox(compression.getScoreAlphabet()));
		// MeasuringCodec<Byte> qsMeasuringCodec = new MeasuringCodec<Byte>(new
		// PPMCodec(3),
		// "Quality score codec");
		MeasuringCodec<Byte> qsMeasuringCodec = new MeasuringCodec<Byte>(new HuffmanByteCodec(qualityScoreTree),
				"Quality score codec");
		final BitCodec<Byte> qualityScoreCodec = qsMeasuringCodec;
		rflNode.add(new DefaultMutableTreeNode(qsMeasuringCodec));

		ReadFeatureCodec readFearureCodec = new ReadFeatureCodec();
		MeasuringCodec<Long> posInReadMeasuringCodec = new MeasuringCodec<Long>(inReadPosCodecStub, "Position in read");
		readFearureCodec.inReadPosCodec = posInReadMeasuringCodec;
		rflNode.add(new DefaultMutableTreeNode(posInReadMeasuringCodec));

		HuffmanTree<Byte> featureOperatorTree = HuffmanCode.buildTree(compression.getReadFeatureFrequencies(),
				Utils.autobox(compression.getReadFeatureAlphabet()));
		HuffmanByteCodec featureOperatorCodec = new HuffmanByteCodec(featureOperatorTree);
		MeasuringCodec<Byte> rfopsMeasuringCodec = new MeasuringCodec<Byte>(featureOperatorCodec,
				"Read feature operators");
		readFearureCodec.featureOperationCodec = rfopsMeasuringCodec;
		rflNode.add(new DefaultMutableTreeNode(rfopsMeasuringCodec));

		// read base node:
		DefaultMutableTreeNode rbNode = new DefaultMutableTreeNode();
		rflNode.add(rbNode);

		ReadBaseCodec readBaseCodec = new ReadBaseCodec();
		MeasuringCodec<Byte> rbBaseMeasuringCodec = new MeasuringCodec<Byte>(baseCodec, "RBQS codec");
		readBaseCodec.baseCodec = rbBaseMeasuringCodec;
		rbNode.add(new DefaultMutableTreeNode(rbBaseMeasuringCodec));

		MeasuringCodec<Byte> rbqsMeasuringCodec = new MeasuringCodec<Byte>(qualityScoreCodec, "RBQS codec");
		readBaseCodec.qualityScoreCodec = rbqsMeasuringCodec;
		rbNode.add(new DefaultMutableTreeNode(rbqsMeasuringCodec));
		MeasuringCodec<ReadBase> readBaseMeasuringCodec = new MeasuringCodec<ReadBase>(readBaseCodec, "Read base codec");
		readFearureCodec.readBaseCodec = readBaseMeasuringCodec;
		rbNode.setUserObject(readBaseMeasuringCodec);

		// deletion node:
		DefaultMutableTreeNode delNode = new DefaultMutableTreeNode();
		rflNode.add(delNode);

		DeletionVariationCodec deletionCodec = new DeletionVariationCodec();
		MeasuringCodec<Long> delLenMeasuringCodec = new MeasuringCodec<Long>(
				createStub(compression.getDelLengthEncoding()), "Deletion length codec");
		deletionCodec.dellengthPosCodec = delLenMeasuringCodec;
		// delNode.add(new DefaultMutableTreeNode(delLenMeasuringCodec));
		MeasuringCodec<DeletionVariation> delMeasuringCodec = new MeasuringCodec<DeletionVariation>(deletionCodec,
				"Deletion codec");
		readFearureCodec.deletionCodec = delMeasuringCodec;
		// delNode.setUserObject(delMeasuringCodec);
		delNode.setUserObject(delLenMeasuringCodec);

		// substitution node:
		DefaultMutableTreeNode subNode = new DefaultMutableTreeNode();
		rflNode.add(subNode);

		SubstitutionVariationCodec substitutionCodec = new SubstitutionVariationCodec();
		MeasuringCodec<BaseChange> baseChangeMeasuringCodec = new MeasuringCodec<BaseChange>(new BaseChangeCodec(),
				"Base change codec");
		substitutionCodec.baseChangeCodec = baseChangeMeasuringCodec;
		// subNode.add(new DefaultMutableTreeNode(baseChangeMeasuringCodec));

		MeasuringCodec<SubstitutionVariation> subMeasuringCodec = new MeasuringCodec<SubstitutionVariation>(
				substitutionCodec, "Substitution codec");
		readFearureCodec.substitutionCodec = subMeasuringCodec;
		// subNode.setUserObject(subMeasuringCodec);
		subNode.setUserObject(baseChangeMeasuringCodec);

		// insertion node:
		DefaultMutableTreeNode insNode = new DefaultMutableTreeNode();
		rflNode.add(insNode);

		InsertionVariationCodec insertionCodec = new InsertionVariationCodec();
		insertionCodec.insertBasesCodec = new ByteArrayHuffmanCodec(compression.getBaseAlphabet(),
				compression.getStopBaseFrequencies(), (byte) '$');
		MeasuringCodec<InsertionVariation> insMeasuringCodec = new MeasuringCodec<InsertionVariation>(insertionCodec,
				"Insertion codec");
		readFearureCodec.insertionCodec = insMeasuringCodec;
		insNode.setUserObject(insMeasuringCodec);

		MeasuringCodec<List<ReadFeature>> varMeasuringCodec = new MeasuringCodec<List<ReadFeature>>(readFearureCodec,
				"Variations codec");
		recordCodec.variationsCodec = varMeasuringCodec;
		rflNode.setUserObject(varMeasuringCodec);

		BitCodec<InsertBase> insertBaseCodec = new BitCodec<InsertBase>() {

			@Override
			public InsertBase read(BitInputStream bis) throws IOException {
				InsertBase ib = new InsertBase();
				ib.setBase(baseCodec.read(bis));
				return ib;
			}

			@Override
			public long write(BitOutputStream bos, InsertBase ib) throws IOException {
				return baseCodec.write(bos, ib.getBase());
			}

			@Override
			public long numberOfBits(InsertBase ib) {
				return baseCodec.numberOfBits(ib.getBase());
			}
		};
		readFearureCodec.insertBaseCodec = new MeasuringCodec<InsertBase>(insertBaseCodec, "Insert base codec");
		rflNode.add(new DefaultMutableTreeNode(readFearureCodec.insertBaseCodec));

		BitCodec<BaseQualityScore> baseQSCodec = new BitCodec<BaseQualityScore>() {

			@Override
			public BaseQualityScore read(BitInputStream bis) throws IOException {
				BaseQualityScore bqs = new BaseQualityScore(-1, qualityScoreCodec.read(bis));
				return bqs;
			}

			@Override
			public long write(BitOutputStream bos, BaseQualityScore bqs) throws IOException {
				return qualityScoreCodec.write(bos, bqs.getQualityScore());
			}

			@Override
			public long numberOfBits(BaseQualityScore bqs) {
				return qualityScoreCodec.numberOfBits(bqs.getQualityScore());
			}
		};
		readFearureCodec.baseQSCodec = new MeasuringCodec<BaseQualityScore>(baseQSCodec, "Base QS codec");
		rflNode.add(new DefaultMutableTreeNode(readFearureCodec.baseQSCodec));

		if (header != null && header.getReadAnnotations() != null && !header.getReadAnnotations().isEmpty()) {
			ReadAnnotation[] anns = new ReadAnnotation[compression.getReadAnnotationIndexes().length];
			for (int i = 0; i < anns.length; i++)
				anns[i] = header.getReadAnnotations().get(i);

			BitCodec<ReadAnnotation> readAnnoCodec = new ReadAnnotationCodec(anns,
					compression.getReadAnnotationFrequencies());
			MeasuringCodec<ReadAnnotation> readAnnoMesuringCodec = new MeasuringCodec<ReadAnnotation>(readAnnoCodec,
					"Read anno codec");
			recordCodec.readAnnoCodec = readAnnoMesuringCodec;
			root.add(new DefaultMutableTreeNode(readAnnoMesuringCodec));
		}

		if (compression.getReadGroupIndexes() == null || compression.getReadGroupIndexes().length < 2) {
			SingleValueBitCodec<Integer> singleValueBitCodec = new SingleValueBitCodec<Integer>();
			singleValueBitCodec.setValue(compression.getReadGroupIndexes()[0]);
			recordCodec.readGroupCodec = singleValueBitCodec;
		} else {
			HuffmanTree<Integer> readGroupIndexTree = HuffmanCode.buildTree(compression.getReadGroupFrequencies(),
					Utils.autobox(compression.getReadGroupIndexes()));
			HuffmanCodec<Integer> readGroupIndexCodec = new HuffmanCodec<Integer>(readGroupIndexTree);
			MeasuringCodec<Integer> measuringReadGroupIndexCodec = new MeasuringCodec<Integer>(readGroupIndexCodec,
					"Read group index codec");
			recordCodec.readGroupCodec = measuringReadGroupIndexCodec;
			root.add(new DefaultMutableTreeNode(measuringReadGroupIndexCodec));
		}

		HuffmanTree<Byte> mappingQualityTree = HuffmanCode.buildTree(compression.getMappingQualityFrequencies(),
				Utils.autobox(compression.getMappingQualityAlphabet()));
		HuffmanByteCodec mappingQualityCodec = new HuffmanByteCodec(mappingQualityTree);
		recordCodec.mappingQualityCodec = new MeasuringCodec<Byte>(mappingQualityCodec, "Mapping quality codec");
		root.add(new DefaultMutableTreeNode(recordCodec.mappingQualityCodec));

		recordCodec.storeMappedQualityScores = block.losslessQualityScores;

		recordCodec.baseCodec = new MeasuringCodec<Byte>(baseCodec, "Record base codec");
		recordCodec.qualityCodec = new MeasuringCodec<Byte>(qualityScoreCodec, "Record quality codec");

		HuffmanTree<Byte> heapByteTree = HuffmanCode.buildTree(compression.getHeapByteFrequencies(),
				Utils.autobox(compression.getHeapByteAlphabet()));
		HuffmanByteCodec heapByteCodec = new HuffmanByteCodec(heapByteTree);
		recordCodec.heapByteCodec = new MeasuringCodec<Byte>(heapByteCodec, "Heap bytes codec");
		root.add(new DefaultMutableTreeNode(recordCodec.heapByteCodec));

		if (compression.tagKeyAlphabet != null && compression.tagKeyAlphabet.length > 0) {
			HuffmanTree<String> tree = HuffmanCode.buildTree(compression.tagKeyFrequency, compression.tagKeyAlphabet);
			recordCodec.tagKeyAndTypeCodec = new MeasuringCodec<String>(new HuffmanCodec<String>(tree), "Tag key codec");
			DefaultMutableTreeNode tagCodecsNode = new DefaultMutableTreeNode(recordCodec.tagKeyAndTypeCodec);
			root.add(tagCodecsNode);

			recordCodec.tagCodecMap = new TreeMap<String, BitCodec<byte[]>>();
			for (int i = 0; i < compression.tagKeyAlphabet.length; i++) {
				String tagKey = compression.tagKeyAlphabet[i];
				ByteFrequencies byteFreqs = compression.tagByteFrequencyMap.get(tagKey);

				IntFrequencies lenFreqs = compression.tagByteLengthMap.get(tagKey);

				int[] lenValues = lenFreqs.getValues();

				BitCodec<byte[]> codec = null;
				if (lenValues.length == 1)
					codec = new FixedLengthByteArrayHuffmanCodec(byteFreqs.getValues(), byteFreqs.getFrequencies(),
							lenValues[0]);
				else {
					codec = new VariableLengthByteArrayHuffmanCodec(byteFreqs.getValues(), byteFreqs.getFrequencies(),
							lenFreqs.getValues(), lenFreqs.getFrequencies());
				}
				MeasuringCodec<byte[]> measuringCodec = new MeasuringCodec<byte[]>(codec, tagKey);
				recordCodec.tagCodecMap.put(tagKey, measuringCodec);
				tagCodecsNode.add(new DefaultMutableTreeNode(measuringCodec));
			}

		}
		HuffmanTree<Byte> flagsTree = HuffmanCode.buildTree(compression.flagStats.getFrequencies(),
				Utils.autobox(compression.flagStats.getValues()));
		HuffmanByteCodec flagsCodec = new HuffmanByteCodec(flagsTree);
		recordCodec.flagsCodec = new MeasuringCodec<Byte>(flagsCodec, "Read flags codec");
		root.add(new DefaultMutableTreeNode(recordCodec.flagsCodec));

		return root;
	}

	public void dump(DefaultMutableTreeNode node, long totalBits) {
		MeasuringCodec<?> codec = (MeasuringCodec<?>) node.getUserObject();
		if (totalBits == 0)
			totalBits = codec.getWrittenBits();

		if (codec.getWrittenBits() == 0 && node.isLeaf())
			return;
		for (int i = 0; i < node.getLevel(); i++)
			System.err.print("\t");

		System.err.println(codec.toString());
		if (!node.isLeaf()) {
			Enumeration children = node.children();
			while (children.hasMoreElements()) {
				DefaultMutableTreeNode childNode = (DefaultMutableTreeNode) children.nextElement();
				dump(childNode, totalBits);
			}
		}
	}

	public Collection<MeasuringCodec> listAllCodecs(DefaultMutableTreeNode node) {
		ArrayList list = Collections.list(node.depthFirstEnumeration());
		ArrayList<MeasuringCodec> codecs = new ArrayList<MeasuringCodec>();
		for (Object o : list)
			codecs.add((MeasuringCodec) ((DefaultMutableTreeNode) o).getUserObject());
		Collections.sort(codecs, byBitsComparator);
		return codecs;
	}

	// private static Comparator<DefaultMutableTreeNode> byBitsComparator = new
	// Comparator<DefaultMutableTreeNode>() {
	//
	// @Override
	// public int compare(DefaultMutableTreeNode o1, DefaultMutableTreeNode o2)
	// {
	// MeasuringCodec c1 = (MeasuringCodec) o1.getUserObject();
	// MeasuringCodec c2 = (MeasuringCodec) o2.getUserObject();
	// return (int) (c1.getWrittenBits() - c2.getWrittenBits());
	// }
	// };

	private static Comparator<MeasuringCodec> byBitsComparator = new Comparator<MeasuringCodec>() {

		@Override
		public int compare(MeasuringCodec o1, MeasuringCodec o2) {
			return -(int) (o1.getWrittenBits() - o2.getWrittenBits());
		}
	};

	public BitCodec<CramRecord> createRecordCodec(CramHeader header, CramRecordBlock block,
			SequenceBaseProvider referenceBaseProvider) throws CramCompressionException {
		CramCompression compression = block.getCompression();
		// given a bunch of block info and compression info create a codec:
		CramRecordCodec recordCodec = new CramRecordCodec();
		recordCodec.prevPosInSeq = block.getFirstRecordPosition();
		recordCodec.inSeqPosCodec = createStub(compression.getInSeqPosEncoding());

		// hack:
		Long[] readLengthAlphabet = new Long[compression.getReadLengthAlphabet().length];
		for (int i = 0; i < readLengthAlphabet.length; i++)
			readLengthAlphabet[i] = new Long(compression.getReadLengthAlphabet()[i]);
		HuffmanTree<Long> readlengthTree = HuffmanCode.buildTree(compression.getReadLengthFrequencies(),
				readLengthAlphabet);
		HuffmanCodec<Long> readlengthCodec = new HuffmanCodec<Long>(readlengthTree);
		// createStub(compression.getReadLengthEncoding());
		recordCodec.readlengthCodec = readlengthCodec;
		recordCodec.recordsToNextFragmentCodec = createStub(compression.getRecordsToNextFragmentEncoding());
		recordCodec.sequenceBaseProvider = referenceBaseProvider;
		recordCodec.defaultReadLength = block.getReadLength();

		NumberCodecStub inReadPosCodecStub = createStub(compression.getInReadPosEncoding());

		HuffmanTree<Byte> baseTree = HuffmanCode.buildTree(compression.getBaseFrequencies(),
				Utils.autobox(compression.getBaseAlphabet()));
		final BitCodec<Byte> baseCodec = new HuffmanByteCodec(baseTree);

		HuffmanTree<Byte> qualityScoreTree = HuffmanCode.buildTree(compression.getScoreFrequencies(),
				Utils.autobox(compression.getScoreAlphabet()));
		final BitCodec<Byte> qualityScoreCodec = new HuffmanByteCodec(qualityScoreTree);

		ReadFeatureCodec readFearureCodec = new ReadFeatureCodec();
		readFearureCodec.inReadPosCodec = inReadPosCodecStub;
		HuffmanTree<Byte> featureOperatorTree = HuffmanCode.buildTree(compression.getReadFeatureFrequencies(),
				Utils.autobox(compression.getReadFeatureAlphabet()));
		HuffmanByteCodec featureOperatorCodec = new HuffmanByteCodec(featureOperatorTree);
		readFearureCodec.featureOperationCodec = featureOperatorCodec;

		ReadBaseCodec readBaseCodec = new ReadBaseCodec();
		readBaseCodec.baseCodec = baseCodec;
		// if (block.isMaskedQualityScoresIncluded())
		readBaseCodec.qualityScoreCodec = qualityScoreCodec;
		// else
		// readBaseCodec.qualityScoreCodec = new NullBitCodec<Byte>();
		readFearureCodec.readBaseCodec = readBaseCodec;

		DeletionVariationCodec deletionCodec = new DeletionVariationCodec();
		deletionCodec.dellengthPosCodec = createStub(compression.getDelLengthEncoding());
		readFearureCodec.deletionCodec = deletionCodec;

		SubstitutionVariationCodec substitutionCodec = new SubstitutionVariationCodec();
		substitutionCodec.baseChangeCodec = new BaseChangeCodec();
		readFearureCodec.substitutionCodec = substitutionCodec;

		InsertionVariationCodec insertionCodec = new InsertionVariationCodec();
		// insertionCodec.insertBasesCodec = new BaseSequenceCodec(
		// BaseSequenceCodec.BaseCodecType.RAISED, "ACGTSN".getBytes());
		insertionCodec.insertBasesCodec = new ByteArrayHuffmanCodec(compression.getBaseAlphabet(),
				compression.getStopBaseFrequencies(), (byte) '$');
		readFearureCodec.insertionCodec = insertionCodec;

		recordCodec.variationsCodec = readFearureCodec;

		// // unampped bases:
		// final ByteArrayHuffmanCodec unmappedBasesCodec = new
		// ByteArrayHuffmanCodec(compression.getBaseAlphabet(),
		// compression.getBaseFrequencies(), (byte) '$');
		// recordCodec.basesCodec = unmappedBasesCodec;
		//
		// // unampped quality scores:
		// if (block.isUnmappedReadQualityScoresIncluded()) {
		// ByteArrayHuffmanCodec unmappedScoresCodec = new
		// ByteArrayHuffmanCodec(compression.getScoreAlphabet(),
		// compression.getScoreFrequencies(), (byte) -1);
		// recordCodec.qualitiesCodec = unmappedScoresCodec;
		// } else
		// recordCodec.qualitiesCodec = new NullBitCodec<byte[]>();

		BitCodec<InsertBase> insertBaseCodec = new BitCodec<InsertBase>() {

			@Override
			public InsertBase read(BitInputStream bis) throws IOException {
				InsertBase ib = new InsertBase();
				ib.setBase(baseCodec.read(bis));
				return ib;
			}

			@Override
			public long write(BitOutputStream bos, InsertBase ib) throws IOException {
				return baseCodec.write(bos, ib.getBase());
			}

			@Override
			public long numberOfBits(InsertBase ib) {
				return baseCodec.numberOfBits(ib.getBase());
			}
		};
		readFearureCodec.insertBaseCodec = insertBaseCodec;

		BitCodec<BaseQualityScore> baseQSCodec = new BitCodec<BaseQualityScore>() {

			@Override
			public BaseQualityScore read(BitInputStream bis) throws IOException {
				BaseQualityScore bqs = new BaseQualityScore(-1, qualityScoreCodec.read(bis));
				return bqs;
			}

			@Override
			public long write(BitOutputStream bos, BaseQualityScore bqs) throws IOException {
				return qualityScoreCodec.write(bos, bqs.getQualityScore());
			}

			@Override
			public long numberOfBits(BaseQualityScore bqs) {
				return qualityScoreCodec.numberOfBits(bqs.getQualityScore());
			}
		};
		readFearureCodec.baseQSCodec = baseQSCodec;

		if (header != null && header.getReadAnnotations() != null && !header.getReadAnnotations().isEmpty()) {
			ReadAnnotation[] anns = new ReadAnnotation[compression.getReadAnnotationIndexes().length];
			for (int i = 0; i < anns.length; i++)
				anns[i] = header.getReadAnnotations().get(i);

			BitCodec<ReadAnnotation> readAnnoCodec = new ReadAnnotationCodec(anns,
					compression.getReadAnnotationFrequencies());
			recordCodec.readAnnoCodec = readAnnoCodec;
		}

		if (compression.getReadGroupIndexes() == null || compression.getReadGroupIndexes().length < 2) {
			SingleValueBitCodec<Integer> singleValueBitCodec = new SingleValueBitCodec<Integer>();
			singleValueBitCodec.setValue(compression.getReadGroupIndexes()[0]);
			recordCodec.readGroupCodec = singleValueBitCodec;
		} else {
			HuffmanTree<Integer> readGroupIndexTree = HuffmanCode.buildTree(compression.getReadGroupFrequencies(),
					Utils.autobox(compression.getReadGroupIndexes()));
			HuffmanCodec<Integer> readGroupIndexCodec = new HuffmanCodec<Integer>(readGroupIndexTree);
			recordCodec.readGroupCodec = readGroupIndexCodec;
		}

		HuffmanTree<Byte> mappingQualityTree = HuffmanCode.buildTree(compression.getMappingQualityFrequencies(),
				Utils.autobox(compression.getMappingQualityAlphabet()));
		HuffmanByteCodec fmappingQualityCodec = new HuffmanByteCodec(mappingQualityTree);
		recordCodec.mappingQualityCodec = fmappingQualityCodec;

		recordCodec.storeMappedQualityScores = block.losslessQualityScores;

		recordCodec.baseCodec = baseCodec;
		recordCodec.qualityCodec = qualityScoreCodec;

		HuffmanTree<Byte> heapByteTree = HuffmanCode.buildTree(compression.getHeapByteFrequencies(),
				Utils.autobox(compression.getHeapByteAlphabet()));
		recordCodec.heapByteCodec = new HuffmanByteCodec(heapByteTree);

		if (compression.tagKeyAlphabet != null && compression.tagKeyAlphabet.length > 0) {
			recordCodec.tagCodecMap = new TreeMap<String, BitCodec<byte[]>>();
			for (int i = 0; i < compression.tagKeyAlphabet.length; i++) {
				String tagKey = compression.tagKeyAlphabet[i];
				ByteFrequencies byteFreqs = compression.tagByteFrequencyMap.get(tagKey);
				IntFrequencies lenFreqs = compression.tagByteLengthMap.get(tagKey);

				int[] lenValues = lenFreqs.getValues();

				BitCodec<byte[]> codec = null;
				if (lenValues.length == 1)
					codec = new FixedLengthByteArrayHuffmanCodec(byteFreqs.getValues(), byteFreqs.getFrequencies(),
							lenValues[0]);
				else {
					codec = new VariableLengthByteArrayHuffmanCodec(byteFreqs.getValues(), byteFreqs.getFrequencies(),
							lenFreqs.getValues(), lenFreqs.getFrequencies());
				}
				recordCodec.tagCodecMap.put(tagKey, codec);
			}

			HuffmanTree<String> tree = HuffmanCode.buildTree(compression.tagKeyFrequency, compression.tagKeyAlphabet);
			recordCodec.tagKeyAndTypeCodec = new HuffmanCodec<String>(tree);
		}

		HuffmanTree<Byte> flagsTree = HuffmanCode.buildTree(compression.flagStats.getFrequencies(),
				Utils.autobox(compression.flagStats.getValues()));
		recordCodec.flagsCodec = new HuffmanByteCodec(flagsTree);

		return recordCodec;
	}

	private static NumberCodecStub createStub(Encoding encoding) throws CramCompressionException {
		if (encoding == null || encoding.getAlgorithm() == EncodingAlgorithm.NULL)
			return NumberCodecFactory.createStub(EncodingAlgorithm.NULL);

		NumberCodecStub stub = NumberCodecFactory.createStub(encoding.getAlgorithm());
		stub.initFromString(encoding.getParameters());
		return stub;
	}

}

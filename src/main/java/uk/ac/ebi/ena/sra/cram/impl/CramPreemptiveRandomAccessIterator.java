package uk.ac.ebi.ena.sra.cram.impl;

import java.io.DataInputStream;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.SeekableFileStream;
import net.sf.samtools.util.SeekableStream;
import uk.ac.ebi.ena.sra.cram.CramException;
import uk.ac.ebi.ena.sra.cram.SequenceBaseProvider;
import uk.ac.ebi.ena.sra.cram.Utils;
import uk.ac.ebi.ena.sra.cram.encoding.BitCodec;
import uk.ac.ebi.ena.sra.cram.encoding.CramRecordCodec;
import uk.ac.ebi.ena.sra.cram.encoding.MeasuringCodec;
import uk.ac.ebi.ena.sra.cram.format.CramFormatException;
import uk.ac.ebi.ena.sra.cram.format.CramHeader;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.CramRecordBlock;
import uk.ac.ebi.ena.sra.cram.format.compression.CramCompressionException;
import uk.ac.ebi.ena.sra.cram.index.CramIndex;
import uk.ac.ebi.ena.sra.cram.index.RecordPointer;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitInputStream;

public class CramPreemptiveRandomAccessIterator implements CloseableIterator<CramRecord> {

	private SeekableStream seekableStream;
	private DataInputStream dis;
	private CramHeader header;
	private CramRecordBlock block;
	private boolean eof = false;

	private long recordCounter = 0;

	private RecordCodecFactory recordCodecFactory = new RecordCodecFactory();
	private SequenceBaseProvider referenceBaseProvider;
	private BitCodec<CramRecord> recordCodec;
	private DefaultBitInputStream bis;
	private RestoreBases restoreBases;
	private RestoreQualityScores restoreScores;
	private final RecordPointer start;
	private DataInputStream blockDIS;

	private CramRecord nextRecord;
	private ReferenceSequenceFile referenceSequenceFile;
	private String refSequenceName;

	private long alStartAddjustment = 0;

	public static void main(String[] args) throws CramException, IOException {
		File cramFile = new File("c:/temp/HG00096.mapped.illumina.mosaik.GBR.exome.20110411.chr20.cram");
		File cramIndexFile = new File(cramFile.getAbsolutePath() + ".crai");
		CramIndex index = CramIndex.fromFile(cramIndexFile);

		File refFile = new File("c:/temp/chr20.fa");
		ReferenceSequenceFile referenceSequenceFile = Utils.createIndexedFastaSequenceFile(refFile);

		// ReferenceSequence nextSequence = referenceSequenceFile
		// .getSequence("20");
		// byte[] refBases = referenceSequenceFile.getSubsequenceAt(
		// nextSequence.getName(), 1, nextSequence.length()).getBases();
		// Utils.capitaliseAndCheckBases(refBases, false);
		// ByteArraySequenceBaseProvider provider = new
		// ByteArraySequenceBaseProvider(
		// refBases);

		Random random = new Random();
		long time1 = System.currentTimeMillis();
		for (int i = 0; i < 1000; i++) {
			long queryStart = random.nextInt(62905065);
			getRecordsStartingFrom(cramFile, index, referenceSequenceFile, "20", queryStart, 1000);
		}
		long time2 = System.currentTimeMillis();
		System.out.println(time2 - time1);
	}

	private static List<CramRecord> getRecordsStartingFrom(File cramFile, CramIndex index,
			ReferenceSequenceFile referenceSequenceFile, String seqName, long queryStart, int nofRecords)
			throws IOException, CramException {

		SeekableFileStream stream = new SeekableFileStream(cramFile);

		List<CramRecord> records = new ArrayList<CramRecord>(nofRecords);

		RecordPointer pointer = index.findRecordPointerAt(seqName, queryStart);
		// System.out.println(pointer);
		CramPreemptiveRandomAccessIterator iterator = new CramPreemptiveRandomAccessIterator(stream,
				referenceSequenceFile, pointer);

		int counter = 0;
		while (iterator.hasNext()) {
			CramRecord record = iterator.next();
			if (record.getAlignmentStart() < queryStart)
				continue;
			records.add(record);
			if (++counter >= nofRecords)
				break;
		}

		return records;
	}

	public CramPreemptiveRandomAccessIterator(SeekableStream stream, ReferenceSequenceFile referenceSequenceFile,
			RecordPointer start) throws IOException, CramException {
		if (!Utils.isCRAM(stream))
			throw new RuntimeException("Not a valid CRAM format.");

		this.seekableStream = stream;
		this.referenceSequenceFile = referenceSequenceFile;
		this.start = start;

		dis = new DataInputStream(stream);
		readHeader();
		seekableStream.seek(start.getBlockStart());
		readNextBlock();
		if (start.getBitOffset() == 0)
			blockDIS.skip(start.getByteOffset());
		else
			blockDIS.skip(start.getByteOffset() - 1);
		bis.readBits(start.getBitOffset());
		recordCounter = start.getRecordNumber();
		CramRecordCodec codec = null;

		alStartAddjustment = start.getAlignmentStart() - block.getFirstRecordPosition();

		advance();

		// ugly hack:
		nextRecord.setAlignmentStart(start.getAlignmentStart());
		if (nextRecord.isReadMapped())
			restoreBases.restoreReadBases(nextRecord);
		if (recordCodec instanceof MeasuringCodec)
			codec = (CramRecordCodec) ((MeasuringCodec) recordCodec).getDelegate();
		else
			codec = (CramRecordCodec) recordCodec;
		codec.prevPosInSeq = start.getAlignmentStart();
	}

	public CramHeader getCramHeader() {
		return header;
	}

	private void readHeader() throws IOException {
		header = CramHeaderIO.read(Utils.getNextChunk(dis));
	}

	private void readNextBlock() throws IOException, CramException {
		blockDIS = Utils.getNextChunk(this.dis);
		if (blockDIS == null) {
			eof = true;
			return;
		}
		CramRecordBlockReader crbReader = new CramRecordBlockReader(blockDIS);
		block = crbReader.read();
		recordCounter = 0;

		if (block == null)
			eof = true;
		else {
			if (refSequenceName == null || !block.getSequenceName().equals(refSequenceName)) {
				refSequenceName = block.getSequenceName();
				byte[] refBases = Utils.getReferenceSequenceBases(referenceSequenceFile, refSequenceName);
				referenceBaseProvider = new ByteArraySequenceBaseProvider(refBases);

			}

			recordCodec = recordCodecFactory.createRecordCodec(header, block, referenceBaseProvider);
			restoreBases = new RestoreBases(referenceBaseProvider, block.getSequenceName());
			restoreScores = new RestoreQualityScores();
			bis = new DefaultBitInputStream(blockDIS);
			alStartAddjustment = 0;
		}
	}

	@Override
	public boolean hasNext() {
		return nextRecord != null;
	}

	@Override
	public CramRecord next() {
		final CramRecord result = nextRecord;
		advance();
		return result;
	}

	@Override
	public void remove() {
		throw new RuntimeException("remove unsupported in CramRecordIterator.");
	}

	private boolean blockHasMoreRecords() {
		return block == null || recordCounter < block.getRecordCount();
	}

	private void advance() {
		if (eof) {
			nextRecord = null;
			return;
		}

		try {
			if (!blockHasMoreRecords())
				readNextBlock();

			if (eof) {
				nextRecord = null;
				return;
			}

			nextRecord = getNextRecord();

		} catch (IOException exc) {
			throw new RuntimeException(exc.getMessage(), exc);
		} catch (CramFormatException e) {
			throw new RuntimeException(e.getMessage(), e);
		} catch (CramCompressionException e) {
			throw new RuntimeException(e.getMessage(), e);
		} catch (CramException e) {
			throw new RuntimeException(e.getMessage(), e);
		}
	}

	/**
	 * Read the next record from the input stream.
	 */
	private CramRecord getNextRecord() throws IOException {
		CramRecord record = recordCodec.read(bis);
		record.setSequenceName(block.getSequenceName());
		// record.setAlignmentStart(record.getAlignmentStart() +
		// alStartAddjustment);
		recordCounter++;

		if (record.isReadMapped())
			restoreBases.restoreReadBases(record);
		restoreScores.restoreQualityScores(record);

		return record;
	}

	/**
	 * @return The record that will be return by the next call to next()
	 */
	protected CramRecord peek() {
		return nextRecord;
	}

	@Override
	public void close() {
		// TODO Auto-generated method stub

	}

}

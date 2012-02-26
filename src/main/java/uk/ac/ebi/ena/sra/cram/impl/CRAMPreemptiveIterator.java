package uk.ac.ebi.ena.sra.cram.impl;

import java.io.DataInputStream;
import java.io.IOException;
import java.io.InputStream;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.samtools.util.CloseableIterator;

import org.apache.log4j.Logger;

import uk.ac.ebi.ena.sra.cram.CramException;
import uk.ac.ebi.ena.sra.cram.SequenceBaseProvider;
import uk.ac.ebi.ena.sra.cram.Utils;
import uk.ac.ebi.ena.sra.cram.encoding.BitCodec;
import uk.ac.ebi.ena.sra.cram.format.CramFormatException;
import uk.ac.ebi.ena.sra.cram.format.CramHeader;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.CramRecordBlock;
import uk.ac.ebi.ena.sra.cram.format.compression.CramCompressionException;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitInputStream;

public class CRAMPreemptiveIterator implements CloseableIterator<CramRecord> {
	private static Logger log = Logger.getLogger(CRAMPreemptiveIterator.class);
	private DataInputStream is;
	private CramHeader header;
	private CramRecordBlock block;
	private boolean eof = false;

	private long recordCounter = 0;

	private RecordCodecFactory recordCodecFactory = new RecordCodecFactory();
	private SequenceBaseProvider referenceBaseProvider;
	private BitCodec<CramRecord> recordCodec;
	private DefaultBitInputStream bis;
	private RestoreBases restoreBases;
	private ReferenceSequenceFile referenceSequenceFile;
	private String refSequenceName;
	private RestoreQualityScores restoreScores;

	private CramRecord mNextRecord;

	public CRAMPreemptiveIterator(InputStream is, ReferenceSequenceFile referenceSequenceFile, CramHeader header)
			throws IOException, CramException {
		if (!Utils.isCRAM(is))
			throw new RuntimeException("Not a valid CRAM format.");

		if (is instanceof DataInputStream)
			this.is = (DataInputStream) is;
		else
			this.is = new DataInputStream(is);
		this.referenceSequenceFile = referenceSequenceFile;
		if (header == null)
			readHeader();

		readNextBlock();
		advance();
	}

	public CramHeader getCramHeader() {
		return header;
	}

	private void readHeader() throws IOException {
		header = CramHeaderIO.read(Utils.getNextChunk(is));
	}

	private void readNextBlock() throws IOException, CramException {
		DataInputStream dis = Utils.getNextChunk(is);
		if (dis == null) {
			eof = true;
			return;
		}
		CramRecordBlockReader crbReader = new CramRecordBlockReader(dis);
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
			bis = new DefaultBitInputStream(dis);
		}
	}

	public void close() {

	}

	public boolean hasNext() {
		return (mNextRecord != null);
	}

	public CramRecord next() {
		final CramRecord result = mNextRecord;
		try {
			advance();
		} catch (CramException e) {
			throw new RuntimeException(e);
		}
		return result;
	}

	public void remove() {
		throw new UnsupportedOperationException("Not supported: remove");
	}

	private boolean blockHasMoreRecords() {
		return block == null || recordCounter < block.getRecordCount();
	}

	private void advance() throws CramException {
		if (eof) {
			mNextRecord = null;
			return;
		}

		try {
			if (!blockHasMoreRecords())
				readNextBlock();

			if (eof || block == null) {
				mNextRecord = null;
				return;
			}

			mNextRecord = getNextRecord();

		} catch (IOException exc) {
			throw new RuntimeException(exc.getMessage(), exc);
		} catch (CramFormatException e) {
			throw new RuntimeException(e.getMessage(), e);
		} catch (CramCompressionException e) {
			throw new RuntimeException(e.getMessage(), e);
		}
	}

	/**
	 * Read the next record from the input stream.
	 */
	private CramRecord getNextRecord() throws IOException {
		CramRecord record = recordCodec.read(bis);
		record.setSequenceName(block.getSequenceName());
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
		return mNextRecord;
	}
}

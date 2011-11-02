package uk.ac.ebi.ena.sra.cram.impl;

import java.io.DataInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.zip.GZIPInputStream;

import net.sf.picard.PicardException;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.samtools.util.CloseableIterator;
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
	private ReferenceSequence referenceSequence;

	private CramRecord mNextRecord;

	public CRAMPreemptiveIterator(InputStream is, ReferenceSequenceFile referenceSequenceFile, CramHeader header)
			throws IOException, CramFormatException, CramCompressionException {
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

	private void readNextBlock() throws CramFormatException, IOException, CramCompressionException {
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
			if (referenceSequence == null || !block.getSequenceName().equals(referenceSequence.getName())) {
					referenceSequence = referenceSequenceFile.getSequence(block.getSequenceName());

				try {
					referenceSequence = referenceSequenceFile.getSubsequenceAt(block.getSequenceName(), 1,
							referenceSequence.length());
				} catch (PicardException e) {
					System.err.println("Reference sequence length: " + referenceSequence.length());
					System.err.println("offensive block: ");
					System.err.println(block.toString());
					throw new RuntimeException(e);
				}
				referenceBaseProvider = new ByteArraySequenceBaseProvider(referenceSequence.getBases());

			}

			recordCodec = recordCodecFactory.createRecordCodec(header, block, referenceBaseProvider);
			restoreBases = new RestoreBases(referenceBaseProvider, block.getSequenceName());
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
		advance();
		return result;
	}

	public void remove() {
		throw new UnsupportedOperationException("Not supported: remove");
	}

	private boolean blockHasMoreRecords() {
		return block == null || recordCounter < block.getRecordCount();
	}

	private void advance() {
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

		return record;
	}

	/**
	 * @return The record that will be return by the next call to next()
	 */
	protected CramRecord peek() {
		return mNextRecord;
	}
}

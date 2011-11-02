package uk.ac.ebi.ena.sra.cram.impl;

import java.io.DataInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.util.zip.GZIPInputStream;

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

public class CramIterator implements CloseableIterator<CramRecord> {

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

	public CramIterator(InputStream is,
			ReferenceSequenceFile referenceSequenceFile) throws IOException {
		this(is, referenceSequenceFile, null);
	}

	public CramIterator(InputStream is,
			ReferenceSequenceFile referenceSequenceFile, CramHeader header)
			throws IOException {
		if (!Utils.isCRAM(is)) throw new RuntimeException("Not a valid CRAM format.") ;
		
		if (is instanceof DataInputStream)
			this.is = (DataInputStream) is;
		else
			this.is = new DataInputStream(is);
		this.referenceSequenceFile = referenceSequenceFile;
		if (header == null)
			readHeader();
	}
	
	public CramHeader getCramHeader () {
		return header ;
	}

	private static InputStream uncompressNextChunk(InputStream is)
			throws IOException {
		int len = 0;
		for (int i = 0; i < 4; i++)
			len = (len << 8) | is.read();

		InputStream limitedIS = new LimitedInputStream(is, len);
		return new GZIPInputStream(limitedIS);
	}

	private static class LimitedInputStream extends InputStream {
		private InputStream delegate;
		private long limit = 0;

		public LimitedInputStream(InputStream delegate, long limit) {
			super();
			this.delegate = delegate;
			this.limit = limit;
		}

		@Override
		public int read() throws IOException {
			if (limit < 1)
				return -1;
			int result = delegate.read();
			limit--;
			return result;
		}

		@Override
		public int read(byte[] b, int off, int len) throws IOException {
			if (limit - len < 0)
				len = (int) limit;

			int result = super.read(b, off, len);
			limit -= result;
			return result;
		}
	}

	private void readHeader() throws IOException {
		header = CramHeaderIO.read(Utils.getNextChunk(is));
	}

	private void readNextBlock() throws CramFormatException, IOException,
			CramCompressionException {
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
			if (referenceSequence == null
					|| !block.getSequenceName().equals(
							referenceSequence.getName())) {
				if (referenceSequence == null)
					referenceSequence = referenceSequenceFile.getSequence(block
							.getSequenceName());

				referenceSequence = referenceSequenceFile.getSubsequenceAt(
						block.getSequenceName(), 1, referenceSequence.length());
				referenceBaseProvider = new ByteArraySequenceBaseProvider(
						referenceSequence.getBases());

			}

			recordCodec = recordCodecFactory.createRecordCodec(header, block,
					referenceBaseProvider);
			restoreBases = new RestoreBases(referenceBaseProvider,
					block.getSequenceName());
			bis = new DefaultBitInputStream(dis);
		}
	}

	public CramRecord skipToAlignmentStart(long alStart) {
		while (hasNext()) {
			CramRecord record = next();
			if (record.getAlignmentStart() >= alStart)
				return record;
		}
		return null;
	}

	public CramRecord skipToAlignmentEnd(long alEnd) {
		while (hasNext()) {
			CramRecord record = next();
			if (record.getAlignmentStart() + record.getReadLength() >= alEnd)
				return record;
		}
		return null;
	}

	@Override
	public boolean hasNext() {
		do {
			if (eof)
				return false;

			if (block == null) {
				try {
					readNextBlock();
				} catch (Exception e) {
					throw new RuntimeException(e);
				}

				if (eof)
					return false;
			}
		} while (recordCounter >= block.getRecordCount());
		return true;
	}

	@Override
	public CramRecord next() {
		if (eof)
			throw new RuntimeException("End of stream.");

		if (block == null || recordCounter >= block.getRecordCount())
			throw new RuntimeException("Iterator depleted or invalid state.");

		CramRecord record;
		try {
			record = recordCodec.read(bis);
			record.setSequenceName(block.getSequenceName()) ;
			recordCounter++;
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

		if (record.isReadMapped())
			try {
				restoreBases.restoreReadBases(record);
			} catch (IOException e) {
				throw new RuntimeException(e);
			}

		if (recordCounter >= block.getRecordCount())
			block = null;

		return record;
	}

	@Override
	public void remove() {
		throw new RuntimeException("remove unsupported in CramRecordIterator.");
	}

	@Override
	public void close() {
	}

}

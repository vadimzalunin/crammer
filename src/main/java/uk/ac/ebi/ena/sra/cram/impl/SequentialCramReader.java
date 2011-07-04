package uk.ac.ebi.ena.sra.cram.impl;

import java.io.DataInputStream;
import java.io.IOException;

import uk.ac.ebi.ena.sra.cram.SequenceBaseProvider;
import uk.ac.ebi.ena.sra.cram.encoding.BitCodec;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.CramRecordBlock;
import uk.ac.ebi.ena.sra.cram.format.compression.CramCompressionException;
import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitInputStream;

public class SequentialCramReader {

	private CramRecordBlockReader blockReader;
	private CramRecordBlock block;

	private BitCodec<CramRecord> recordCodec;
	private BitInputStream bis;

	private long awaitingRecords = -1L;

	private DataInputStream dis;
	private RecordCodecFactory recordCodecFactory = new RecordCodecFactory();
	private SequenceBaseProvider referenceBaseProvider;

	private RestoreBases restoreBases;

	public SequentialCramReader(DataInputStream dis,
			SequenceBaseProvider referenceBaseProvider) {
		this.dis = dis;
		this.referenceBaseProvider = referenceBaseProvider;
		blockReader = new CramRecordBlockReader(dis);
	}

	public CramRecordBlock readBlock() throws IOException,
			CramCompressionException {
		if (awaitingRecords > 0L)
			throw new RuntimeException("Pending records found. ");
		block = blockReader.read();
		if (block == null) {
			awaitingRecords = -1;
			referenceBaseProvider = null;
			recordCodec = null;
			restoreBases = null;
			return null;
		}
		awaitingRecords = block.getRecordCount();
		bis = new DefaultBitInputStream(dis);

		if (referenceBaseProvider != null) {
			recordCodec = recordCodecFactory.createRecordCodec(block,
					referenceBaseProvider);
			restoreBases = new RestoreBases(referenceBaseProvider,
					block.getSequenceName());
		}
		return block;
	}

	public CramRecord readRecord() throws IOException {
		if (block == null)
			throw new RuntimeException("Read block before reading records.");
		if (awaitingRecords < 1)
			throw new RuntimeException("No more records in this block.");
		CramRecord record = recordCodec.read(bis);
		if (record.isReadMapped())
			restoreBases.restoreReadBases(record);
		awaitingRecords--;
		return record;
	}

	public SequenceBaseProvider getReferenceBaseProvider() {
		return referenceBaseProvider;
	}

	public void setReferenceBaseProvider(
			SequenceBaseProvider referenceBaseProvider)
			throws CramCompressionException {
		this.referenceBaseProvider = referenceBaseProvider;
		if (block != null) {
			recordCodec = recordCodecFactory.createRecordCodec(block,
					referenceBaseProvider);
			restoreBases = new RestoreBases(referenceBaseProvider,
					block.getSequenceName());
		}
	}
}

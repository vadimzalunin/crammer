package uk.ac.ebi.ena.sra.cram.impl;

import java.io.IOException;
import java.io.OutputStream;

import uk.ac.ebi.ena.sra.cram.SequenceBaseProvider;
import uk.ac.ebi.ena.sra.cram.encoding.BitCodec;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.CramRecordBlock;
import uk.ac.ebi.ena.sra.cram.format.compression.CramCompressionException;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitOutputStream;

public class SequentialCramWriter {
	private RecordCodecFactory recordCodecFactory = new RecordCodecFactory();
	private CramRecordBlockWriter blockWriter;
	private BitOutputStream bos;
	private BitCodec<CramRecord> recordCodec;
	private OutputStream os ;
	private SequenceBaseProvider referenceBaseProvider ;

	public SequentialCramWriter(OutputStream os, SequenceBaseProvider referenceBaseProvider) {
		this.os = os ;
		this.referenceBaseProvider = referenceBaseProvider ;
		bos = new DefaultBitOutputStream(os);
		blockWriter = new CramRecordBlockWriter(os);
	}

	public void write(CramRecordBlock block) throws IOException, CramCompressionException {
		bos.flush();
		recordCodec = recordCodecFactory.createRecordCodec(block, referenceBaseProvider);
		blockWriter.write(block);
	}

	public void write(CramRecord record) throws IOException {
		recordCodec.write(bos, record);
	}

	public void flush() throws IOException {
		bos.flush() ;
		os.flush() ;
	}
	
}

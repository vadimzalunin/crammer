package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.DataOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.Iterator;

import uk.ac.ebi.ena.sra.cram.CramWriter;
import uk.ac.ebi.ena.sra.cram.format.CramFormatException;
import uk.ac.ebi.ena.sra.cram.format.CramHeader;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.CramRecordBlock;
import uk.ac.ebi.ena.sra.cram.impl.CramRecordWriter;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitOutputStream;

public class JavaBinaryCramWriter implements CramWriter {
	private CramHeader cram;
	private DataOutputStream dos;
	private CramRecordWriter writer;
	private CramRecordBlock currentbBlock;
	private Iterator<CramRecordBlock> blockIterator;
	private int recordCounter;

	public JavaBinaryCramWriter(DataOutputStream dos) throws IOException,
			CramFormatException {
		this.dos = dos;
		writer = new CramRecordWriter(new DefaultBitOutputStream(dos));
	}

	@Override
	public void setCram(CramHeader cram) throws IOException {
		if (this.cram != null)
			throw new RuntimeException("CRAM header is already written.");

		this.cram = cram;
		ObjectOutputStream oos = new ObjectOutputStream(dos);
		oos.writeObject(cram);
		oos.flush();

		blockIterator = cram.getBlocks().iterator();
		currentbBlock = blockIterator.next();
		recordCounter = 0;
	}

	@Override
	public void addRecord(CramRecord record) throws IOException {
		writer.writeCramRecord(record);
		if (recordCounter++ >= currentbBlock.getRecordCount()) {
			if (!blockIterator.hasNext())
				throw new RuntimeException("Record is beyond last block.");
			currentbBlock = blockIterator.next();
			recordCounter = 0;
		}
	}

	@Override
	public void close() throws IOException {
		writer.flush();
	}

}

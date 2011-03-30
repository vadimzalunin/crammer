package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.DataInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.Iterator;

import uk.ac.ebi.ena.sra.cram.CramReader;
import uk.ac.ebi.ena.sra.cram.SequenceBaseProvider;
import uk.ac.ebi.ena.sra.cram.format.CramFormatException;
import uk.ac.ebi.ena.sra.cram.format.CramHeader;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.CramRecordBlock;
import uk.ac.ebi.ena.sra.cram.impl.CramRecordIterator;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitInputStream;

public class JavaBinaryCramReader implements CramReader {
	private SequenceBaseProvider sequenceBaseProvider;
	private CramRecordIterator iterator;
	private CramHeader cram;
	private long recordCounter;
	private CramRecordBlock currentBlock;
	private Iterator<CramRecordBlock> blockIterator;

	public JavaBinaryCramReader(SequenceBaseProvider sequenceBaseProvider,
			DataInputStream dis) throws IOException, CramFormatException {
		this.sequenceBaseProvider = sequenceBaseProvider;
		readCram(dis);
	}

	private void readCram(DataInputStream dis) throws IOException,
			CramFormatException {
		ObjectInputStream ois = new ObjectInputStream(dis);
		try {
			Object object = ois.readObject();
			if (object == null)
				throw new CramFormatException(
						"Failed to read any Cram objects.");

			if (!(object instanceof CramHeader))
				throw new CramFormatException(
						"Deserialized object is not Cram: "
								+ object.getClass().getName());

			cram = (CramHeader) object;
			currentBlock = cram.getBlocks().iterator().next();
			blockIterator = cram.getBlocks().iterator();
			recordCounter = 0;

			iterator = new CramRecordIterator(new DefaultBitInputStream(dis),
					sequenceBaseProvider);
			iterator.setSequenceName(currentBlock.getSequenceName());
		} catch (ClassNotFoundException e) {
			throw new CramFormatException(
					"Expecting a java serialised Cram object.");
		}
	}

	@Override
	public boolean hasNext() {
		if (recordCounter > currentBlock.getRecordCount()
				&& !blockIterator.hasNext())
			return false;
		return iterator.hasNext();
	}

	@Override
	public CramRecord next() {
		recordCounter++;
		if (recordCounter > currentBlock.getRecordCount()) {
			currentBlock = blockIterator.next();
			recordCounter = 0;
			iterator.setSequenceName(currentBlock.getSequenceName());
		}
		CramRecord record = iterator.next();

		// fill out block values:
		record.setReadLength(currentBlock.getReadLength());
		return record;
	}

	@Override
	public void remove() {
		throw new RuntimeException(
				"Changes to CRAM record iterator are not allowed.");
	}

	@Override
	public CramHeader getCram() {
		return cram;
	}

	public CramRecordBlock getCurrentBlock() {
		return currentBlock;
	}

}

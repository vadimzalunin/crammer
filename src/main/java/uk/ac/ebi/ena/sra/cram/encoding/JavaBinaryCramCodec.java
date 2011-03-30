package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.List;

import uk.ac.ebi.ena.sra.cram.SequenceBaseProvider;
import uk.ac.ebi.ena.sra.cram.format.CramFormatException;
import uk.ac.ebi.ena.sra.cram.format.CramHeader;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.CramRecordBlock;
import uk.ac.ebi.ena.sra.cram.impl.CramRecordIterator;
import uk.ac.ebi.ena.sra.cram.impl.CramRecordWriter;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitInputStream;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitOutputStream;

public class JavaBinaryCramCodec implements Codec<CramHeader> {
	private SequenceBaseProvider sequenceBaseProvider;

	// ugly hack:
	private int maxRecords = -1;

	public JavaBinaryCramCodec(SequenceBaseProvider sequenceBaseProvider) {
		super();
		this.sequenceBaseProvider = sequenceBaseProvider;
	}

	@Override
	public CramHeader read(DataInputStream dis) throws IOException,
			CramFormatException {
		ObjectInputStream ois = new ObjectInputStream(dis);
		CramHeader cram = null;
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

			CramRecordIterator reader = new CramRecordIterator(
					new DefaultBitInputStream(dis), sequenceBaseProvider);
			long recordCounter = 0;
			for (CramRecordBlock block : cram.getBlocks()) {
				List<CramRecord> blockRecords = new ArrayList<CramRecord>();
				reader.setSequenceName(block.getSequenceName());
				for (int i = 0; i < block.getRecordCount(); i++) {
					if (!reader.hasNext())
						throw new CramFormatException(
								"Incomplete alignmet sections, trailing records not found.");
					CramRecord record = reader.next();
					record.setReadLength(block.getReadLength());
					blockRecords.add(record);
					if (maxRecords > -1 && recordCounter++ >= maxRecords)
						break;
				}

				block.setRecords(blockRecords);
				if (maxRecords > -1 && recordCounter++ >= maxRecords)
					break;
			}
		} catch (ClassNotFoundException e) {
			throw new CramFormatException(
					"Expecting a java serialised Cram object.");
		}
		return cram;
	}

	@Override
	public void write(DataOutputStream dos, CramHeader cram) throws IOException {
		ObjectOutputStream oos = new ObjectOutputStream(dos);
		oos.writeObject(cram);
		oos.flush();
		CramRecordWriter writer = new CramRecordWriter(
				new DefaultBitOutputStream(dos));
		for (CramRecordBlock block : cram.getBlocks()) 
			for (CramRecord record : block.getRecords()) 
				writer.writeCramRecord(record);
			
		writer.flush();
	}

	public int getMaxRecords() {
		return maxRecords;
	}

	public void setMaxRecords(int maxRecords) {
		this.maxRecords = maxRecords;
	}

}

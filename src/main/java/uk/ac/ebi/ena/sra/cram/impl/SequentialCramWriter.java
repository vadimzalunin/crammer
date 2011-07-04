package uk.ac.ebi.ena.sra.cram.impl;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;

import javax.swing.tree.DefaultMutableTreeNode;

import uk.ac.ebi.ena.sra.cram.CramException;
import uk.ac.ebi.ena.sra.cram.SequenceBaseProvider;
import uk.ac.ebi.ena.sra.cram.encoding.BitCodec;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.CramRecordBlock;
import uk.ac.ebi.ena.sra.cram.format.compression.CramCompressionException;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitInputStream;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitOutputStream;

public class SequentialCramWriter {
	private RecordCodecFactory recordCodecFactory = new RecordCodecFactory();
	private CramRecordBlockWriter blockWriter;
	private BitOutputStream bos;
	private BitCodec<CramRecord> recordCodec;
	private OutputStream os;
	private SequenceBaseProvider referenceBaseProvider;

	private ExposedByteArrayOutputStream checkBAOS = new ExposedByteArrayOutputStream(
			1024);
	private BitOutputStream checkBOS = new DefaultBitOutputStream(checkBAOS);
	private ByteArrayInputStream checkBAIS = new ByteArrayInputStream(
			checkBAOS.getBuffer());
	private DefaultBitInputStream checkBIS = new DefaultBitInputStream(
			checkBAIS);
	private BitCodec<CramRecord> checkRecordReadCodec;
	private BitCodec<CramRecord> checkRecordWriteCodec;
	private RestoreBases restoreBases;

	private DefaultMutableTreeNode rootNode;

	public SequentialCramWriter(OutputStream os,
			SequenceBaseProvider referenceBaseProvider) {
		this.os = os;
		this.referenceBaseProvider = referenceBaseProvider;
		bos = new DefaultBitOutputStream(os);
		blockWriter = new CramRecordBlockWriter(os);
	}

	public long write(CramRecordBlock block) throws IOException,
			CramCompressionException {
		bos.flush();

		restoreBases = new RestoreBases(referenceBaseProvider,
				block.getSequenceName());

		rootNode = recordCodecFactory.buildCodecTree(block,
				referenceBaseProvider);
		checkRecordReadCodec = (BitCodec<CramRecord>) rootNode.getUserObject();

		rootNode = recordCodecFactory.buildCodecTree(block,
				referenceBaseProvider);
		checkRecordWriteCodec = (BitCodec<CramRecord>) rootNode.getUserObject();

		rootNode = recordCodecFactory.buildCodecTree(block,
				referenceBaseProvider);
		recordCodec = (BitCodec<CramRecord>) rootNode.getUserObject();

		return blockWriter.write(block);
	}

	public long write(CramRecord record) throws IOException {
		return recordCodec.write(bos, record);
	}

	public long writeAndCheck(CramRecord record) throws IOException,
			CramException {
		checkBAOS.reset();
		checkRecordWriteCodec.write(checkBOS, record);
		checkBOS.flush();

		checkBAIS.reset();
		checkBIS.reset();
		CramRecord checkRecord = checkRecordReadCodec.read(checkBIS);
		restoreBases.restoreReadBases(checkRecord);

		if (!checkRecord.equals(record)) {
			System.err.println(record.toString());
			System.err.println(checkRecord.toString());
			checkBAIS.reset();
			checkRecord = checkRecordReadCodec.read(checkBIS);

			throw new CramException("Round trip check failed.");
		}

		return write(record);
	}

	public void flush() throws IOException {
		bos.flush();
		os.flush();
	}

	public void dump() {
		recordCodecFactory.dump(rootNode);
	}

	public static class ExposedByteArrayOutputStream extends
			ByteArrayOutputStream {

		public ExposedByteArrayOutputStream(int size) {
			super(size);
		}

		public byte[] getBuffer() {
			return buf;
		}
	}
}

package uk.ac.ebi.ena.sra.cram.impl;

import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Collection;

import org.apache.log4j.Logger;

import uk.ac.ebi.ena.sra.cram.CramException;
import uk.ac.ebi.ena.sra.cram.SequenceBaseProvider;
import uk.ac.ebi.ena.sra.cram.format.CramHeader;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.CramRecordBlock;
import uk.ac.ebi.ena.sra.cram.format.CramReferenceSequence;
import uk.ac.ebi.ena.sra.cram.stats.CramStats;

public class CramWriter {
	private OutputStream os;
	private SequentialCramWriter writer;
	private SequenceBaseProvider provider;
	private CramStats stats;
	private CramRecordBlock block;
	private int maxRecordsPerBlock = 100000;
	private Collection<CramReferenceSequence> sequences;
	private boolean roundTripCheck;
	private long recordsBitLength;

	private static Logger log = Logger.getLogger(CramWriter.class);

	public CramWriter(OutputStream os, SequenceBaseProvider provider,
			Collection<CramReferenceSequence> sequences,
			boolean roundTripCheck, int maxBlockSize) {
		this.os = os;
		this.provider = provider;
		this.sequences = sequences;
		this.roundTripCheck = roundTripCheck;
		maxRecordsPerBlock = maxBlockSize;
	}

	public void dump() {
		writer.dump();
	}

	public void init() throws IOException {
		CramHeader header = new CramHeader();
		header.setVersion("0.3");
		header.setReferenceSequences(sequences);
		CramHeaderIO.write(header, os);
		stats = new CramStats();
		recordsBitLength = 0;
	}

	private CramReferenceSequence findSequenceByName(String name) {
		for (CramReferenceSequence sequence : sequences) {
			if (name.equals(sequence.getName()))
				return sequence;
		}

		return null;
	}

	public void startSequence(String name, byte[] bases) throws IOException,
			CramException {
		CramReferenceSequence sequence = findSequenceByName(name);
		if (sequence == null)
			throw new IllegalArgumentException(
					"Unknown refernce sequence name: " + name);

		recordsBitLength += purgeBlock(block);

		provider = new ByteArraySequenceBaseProvider(bases);
		writer = new SequentialCramWriter(os, provider);

		block.setSequenceName(sequence.getName());
		block.setSequenceLength(sequence.getLength());
		stats = new CramStats();
	}

	public void addRecord(CramRecord record) throws IOException, CramException {
		stats.addRecord(record);
		block.getRecords().add(record);

		if (block.getRecords().size() >= maxRecordsPerBlock) {
			recordsBitLength += purgeBlock(block);
		}
	}

	private long purgeBlock(CramRecordBlock block) throws IOException,
			CramException {
		long len = 0;
		if (block != null && block.getRecords() != null
				&& !block.getRecords().isEmpty()) {
			stats.adjustBlock(block);
			block.setRecordCount(block.getRecords().size());
			log.debug(block.toString());
			len = write(block);
			log.debug(String.format("Block purged: %d\t%d\t%.2f\t%.4f",
					block.getRecordCount(), len,
					(float) len / block.getRecordCount(),
					(float) len / stats.getBaseCount()));
		}

		stats = new CramStats();
		this.block = new CramRecordBlock();
		if (block != null) {
			this.block.setSequenceName(block.getSequenceName());
			this.block.setSequenceLength(block.getSequenceLength());
		}
		this.block
				.setRecords(new ArrayList<CramRecord>(maxRecordsPerBlock + 1));

		return len;
	}

	private long write(CramRecordBlock block) throws IOException, CramException {
		long len = 0;
		len = writer.write(block);
		len = 0;

		if (roundTripCheck)
			for (CramRecord record : block.getRecords())
				len += writer.writeAndCheck(record);
		else
			for (CramRecord record : block.getRecords())
				len += writer.write(record);

		writer.flush();

		return len;
	}

	public void close() throws IOException, CramException {
		recordsBitLength += purgeBlock(block);
	}

	public long getRecordsBitLength() {
		return recordsBitLength;
	}
}
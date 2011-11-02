package uk.ac.ebi.ena.sra.cram.impl;

import java.io.DataOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPOutputStream;

import org.apache.log4j.Logger;

import uk.ac.ebi.ena.sra.cram.CramException;
import uk.ac.ebi.ena.sra.cram.SequenceBaseProvider;
import uk.ac.ebi.ena.sra.cram.format.CramHeader;
import uk.ac.ebi.ena.sra.cram.format.CramReadGroup;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.CramRecordBlock;
import uk.ac.ebi.ena.sra.cram.format.CramReferenceSequence;
import uk.ac.ebi.ena.sra.cram.format.ReadAnnotation;
import uk.ac.ebi.ena.sra.cram.io.ExposedByteArrayOutputStream;
import uk.ac.ebi.ena.sra.cram.stats.CramStats;

public class CramWriter {
	private OutputStream os;
	private ExposedByteArrayOutputStream writerOS;
	private SequentialCramWriter writer;
	private SequenceBaseProvider provider;
	private CramStats stats;
	private CramRecordBlock block;
	private int maxRecordsPerBlock = 1000000;
	private List<CramReferenceSequence> sequences;
	private boolean roundTripCheck;
	private long bitsWritten;

	private static Logger log = Logger.getLogger(CramWriter.class);
	private final boolean captureUnammpedQualityScortes;
	private final boolean captureSubstituionQualityScore;
	private final boolean captureMaskedQualityScores;
	private boolean autodump;
	private CramHeader header;
	private final List<ReadAnnotation> readAnnotations;
	private final PrintStream statsPS;
	private final List<CramReadGroup> cramReadGroups;

	public CramWriter(OutputStream os, SequenceBaseProvider provider,
			List<CramReferenceSequence> sequences, boolean roundTripCheck,
			int maxRecordsPerBlock, boolean captureUnammpedQualityScortes,
			boolean captureSubstituionQualityScore,
			boolean captureMaskedQualityScores,
			List<ReadAnnotation> readAnnotations, PrintStream statsPS,
			List<CramReadGroup> cramReadGroups) {
		this.os = os;
		this.provider = provider;
		this.sequences = sequences;
		this.roundTripCheck = roundTripCheck;
		this.maxRecordsPerBlock = maxRecordsPerBlock;
		this.captureUnammpedQualityScortes = captureUnammpedQualityScortes;
		this.captureSubstituionQualityScore = captureSubstituionQualityScore;
		this.captureMaskedQualityScores = captureMaskedQualityScores;
		this.readAnnotations = readAnnotations;
		this.statsPS = statsPS;
		this.cramReadGroups = cramReadGroups;
	}

	public void dump() {
		writer.dump();
	}

	private void flushWriterOS() throws IOException {
		// stream for in-memory compressed data:
		ExposedByteArrayOutputStream compressedOS = new ExposedByteArrayOutputStream(
				1024);
		GZIPOutputStream gos = new GZIPOutputStream(compressedOS);

		// compress data into memory:
		gos.write(writerOS.getBuffer(), 0, writerOS.size());
		gos.close();

		writerOS.reset();

		// write compresseed data size:
		DataOutputStream dos = new DataOutputStream(os);
		dos.writeInt(compressedOS.size());

		// write compressed data:
		dos.write(compressedOS.getBuffer(), 0, compressedOS.size());
		dos.flush();
	}

	public void init() throws IOException {
		// magick:
		os.write("CRAM".getBytes()) ;
		
		writerOS = new ExposedByteArrayOutputStream(1024 * 1024 * 10);

		header = new CramHeader();
		header.setVersion("0.5");
		header.setReferenceSequences(sequences);
		header.setReadAnnotations(readAnnotations);
		header.setReadGroups(cramReadGroups);
		CramHeaderIO.write(header, writerOS);
		flushWriterOS();
		stats = new CramStats(header, statsPS);
		bitsWritten = 0;
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
					"Unknown reference sequence name: " + name);

		bitsWritten += purgeBlock(block);

		provider = new ByteArraySequenceBaseProvider(bases);
		writer = new SequentialCramWriter(writerOS, provider, header);

		block.setSequenceName(sequence.getName());
		block.setSequenceLength(sequence.getLength());
		stats = new CramStats(header, statsPS);
	}

	public void addRecord(CramRecord record) throws IOException, CramException {
		stats.addRecord(record);
		block.getRecords().add(record);

		if (block.getRecords().size() >= maxRecordsPerBlock)
			bitsWritten += purgeBlock(block);
	}

	private long purgeBlock(CramRecordBlock block) throws IOException,
			CramException {
		long len = 0;
		if (block != null && block.getRecords() != null
				&& !block.getRecords().isEmpty()) {
			if (statsPS != null)
				statsPS.printf("Sequence name: %s; ref length=%d\n",
						block.getSequenceName(), block.getSequenceLength());
			stats.adjustBlock(block);
			block.setRecordCount(block.getRecords().size());
			log.debug(block.toString());
			len = write(block);
			log.debug(String.format("Block purged: %d\t%d\t%.2f\t%.4f",
					block.getRecordCount(), len,
					(float) len / block.getRecordCount(),
					(float) len / stats.getBaseCount()));
		}

		stats = new CramStats(header, statsPS);
		this.block = new CramRecordBlock();
		if (block != null) {
			this.block.setSequenceName(block.getSequenceName());
			this.block.setSequenceLength(block.getSequenceLength());
		}
		this.block
				.setRecords(new ArrayList<CramRecord>(maxRecordsPerBlock + 1));
		this.block
				.setUnmappedReadQualityScoresIncluded(captureUnammpedQualityScortes);
		this.block
				.setSubstitutionQualityScoresIncluded(captureSubstituionQualityScore);
		this.block.setMaskedQualityScoresIncluded(captureMaskedQualityScores);

		return len;
	}

	private long write(CramRecordBlock block) throws IOException, CramException {
		long len = 0;
		len = writer.write(block);

		if (roundTripCheck)
			for (CramRecord record : block.getRecords())
				len += writer.writeAndCheck(record);
		else
			for (CramRecord record : block.getRecords())
				len += writer.write(record);

		writer.flush();
		flushWriterOS();

		if (autodump)
			dump();

		return len;
	}

	public void close() throws IOException, CramException {
		bitsWritten += purgeBlock(block);
	}

	public long getBitsWritten() {
		return bitsWritten;
	}

	public boolean isAutodump() {
		return autodump;
	}

	public void setAutodump(boolean autodump) {
		this.autodump = autodump;
	}
}
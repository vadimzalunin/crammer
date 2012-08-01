/*******************************************************************************
 * Copyright 2012 EMBL-EBI, Hinxton outstation
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *   http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 ******************************************************************************/
package uk.ac.ebi.ena.sra.cram.impl;

import java.io.DataOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.zip.GZIPOutputStream;

import net.sf.samtools.SAMRecord;

import org.apache.log4j.Logger;

import uk.ac.ebi.ena.sra.cram.CramException;
import uk.ac.ebi.ena.sra.cram.SequenceBaseProvider;
import uk.ac.ebi.ena.sra.cram.format.CramCompression;
import uk.ac.ebi.ena.sra.cram.format.CramHeader;
import uk.ac.ebi.ena.sra.cram.format.CramHeaderRecord;
import uk.ac.ebi.ena.sra.cram.format.CramReadGroup;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.CramRecordBlock;
import uk.ac.ebi.ena.sra.cram.format.CramReferenceSequence;
import uk.ac.ebi.ena.sra.cram.format.ReadAnnotation;
import uk.ac.ebi.ena.sra.cram.impl.RecordCodecFactory.CodecStats;
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
	private final boolean captureAllQS;

	private long blockCreationTime = -1;
	private long beyondHorizon = 0;
	private long extraChromosomePairs = 0;
	private final List<CramHeaderRecord> headerRecords;

	private CodecStats codecStats;
	private final boolean preserveReadNames;

	public CramWriter(OutputStream os, SequenceBaseProvider provider, List<CramReferenceSequence> sequences,
			boolean roundTripCheck, int maxRecordsPerBlock, boolean captureUnammpedQualityScortes,
			boolean captureSubstituionQualityScore, boolean captureMaskedQualityScores,
			List<ReadAnnotation> readAnnotations, PrintStream statsPS, List<CramReadGroup> cramReadGroups,
			boolean captureAllQS, List<CramHeaderRecord> headerRecords, boolean preserveReadNames) {
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
		this.captureAllQS = captureAllQS;
		this.headerRecords = headerRecords;
		this.preserveReadNames = preserveReadNames;
	}

	public void dump() {
		writer.dump();
	}

	private void flushWriterOS() throws IOException {
		// stream for in-memory compressed data:
		ExposedByteArrayOutputStream compressedOS = new ExposedByteArrayOutputStream(1024);
		// BZip2CompressorOutputStream gos = new
		// BZip2CompressorOutputStream(compressedOS);
		GZIPOutputStream gos = new GZIPOutputStream(compressedOS);

		// compress data into memory:
		gos.write(writerOS.getBuffer(), 0, writerOS.size());
		gos.close();

		writerOS.reset();

		// write compressed data size:
		DataOutputStream dos = new DataOutputStream(os);
		dos.writeInt(compressedOS.size());

		// write compressed data:
		dos.write(compressedOS.getBuffer(), 0, compressedOS.size());
		dos.flush();
	}

	public void init() throws IOException {
		// magick:
		os.write("CRAM".getBytes());

		writerOS = new ExposedByteArrayOutputStream(1024 * 1024 * 10);

		header = new CramHeader();
		header.setRecords(headerRecords);
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

	public void startSequence(String name, byte[] bases) throws IOException, CramException {
		if (SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(name)) {

			bitsWritten += purgeBlock(block);

			provider = new ByteArraySequenceBaseProvider(new byte[] {});
			writer = new SequentialCramWriter(writerOS, provider, header);

			block.setSequenceName(SAMRecord.NO_ALIGNMENT_REFERENCE_NAME);
			block.setSequenceLength(0);
			stats = new CramStats(header, statsPS, captureAllQS ? 1F : 0F);
		} else {
			CramReferenceSequence sequence = findSequenceByName(name);
			if (sequence == null)
				throw new IllegalArgumentException("Unknown reference sequence name: " + name);

			bitsWritten += purgeBlock(block);

			provider = new ByteArraySequenceBaseProvider(bases);
			writer = new SequentialCramWriter(writerOS, provider, header);

			block.setSequenceName(sequence.getName());
			block.setSequenceLength(sequence.getLength());
			stats = new CramStats(header, statsPS, captureAllQS ? 1F : 0F);
		}
	}

	public void addRecord(CramRecord record) throws IOException, CramException {
		if (record.next != null || record.previous != null) {
			CramRecord mate = record.next == null ? record.previous : record.next;
			if (block.getSequenceName().equals(mate.getSequenceName()))
				beyondHorizon++;
			else
				extraChromosomePairs++;
		}
		stats.addRecord(record);
		block.getRecords().add(record);

		if (block.getRecords().size() >= maxRecordsPerBlock)
			bitsWritten += purgeBlock(block);
	}

	private long purgeBlock(CramRecordBlock block) throws IOException, CramException {
		long len = 0;
		if (block != null && block.getRecords() != null && !block.getRecords().isEmpty()) {
			if (statsPS != null)
				statsPS.printf("Sequence name: %s; ref length=%d\n", block.getSequenceName(), block.getSequenceLength());
			stats.adjustBlock(block);
			block.setRecordCount(block.getRecords().size());
			log.debug(block.toString());
			len = write(block);
			log.info(String.format("Block purged: %s\t%d\t%d\t%d\t%.2f\t%.4f\t%.3fs\t%d\t%d", block.getSequenceName(),
					block.getRecordCount(), writer.gzippedBlockHeaderBytes, len, (float) len / block.getRecordCount(),
					(float) len / stats.getBaseCount(), (System.currentTimeMillis() - blockCreationTime) / 1000f,
					beyondHorizon, extraChromosomePairs));
		}

		beyondHorizon = 0;
		extraChromosomePairs = 0;
		blockCreationTime = System.currentTimeMillis();
		stats = new CramStats(header, statsPS, captureAllQS ? 1F : 0F);
		this.block = new CramRecordBlock();
		if (block != null) {
			this.block.setSequenceName(block.getSequenceName());
			this.block.setSequenceLength(block.getSequenceLength());
		}
		this.block.setRecords(new ArrayList<CramRecord>(maxRecordsPerBlock + 1));
		this.block.setUnmappedReadQualityScoresIncluded(captureUnammpedQualityScortes);
		this.block.setSubstitutionQualityScoresIncluded(captureSubstituionQualityScore);
		this.block.setMaskedQualityScoresIncluded(captureMaskedQualityScores);
		this.block.preserveReadNames = preserveReadNames ;
		if (this.block.getCompression() == null)
			this.block.setCompression(new CramCompression());
		this.block.losslessQualityScores = captureAllQS;

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

		if (codecStats == null)
			codecStats = writer.getCodecStats();
		else
			codecStats.add(writer.getCodecStats());

		return len;
	}

	public void close() throws IOException, CramException {
		bitsWritten += purgeBlock(block);
	}

	public CodecStats getCodecStats() {
		return codecStats ;
	}

	public boolean isAutodump() {
		return autodump;
	}

	public void setAutodump(boolean autodump) {
		this.autodump = autodump;
	}
}

package net.sf.samtools;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import uk.ac.ebi.ena.sra.cram.CramException;
import uk.ac.ebi.ena.sra.cram.Utils;
import uk.ac.ebi.ena.sra.cram.bam.Sam2CramRecordFactory;
import uk.ac.ebi.ena.sra.cram.bam.Sam2CramRecordFactory.TREAT_TYPE;
import uk.ac.ebi.ena.sra.cram.format.CramReadGroup;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.CramReferenceSequence;
import uk.ac.ebi.ena.sra.cram.impl.CramWriter;
import uk.ac.ebi.ena.sra.cram.mask.RefMaskUtils;
import uk.ac.ebi.ena.sra.cram.spot.PairedTemplateAssembler;

public class CRAMFileWriter implements SAMFileWriter {
	private ReferenceSequenceFile referenceSequenceFile;
	private OutputStream os;
	private SAMFileHeader samHeader;
	private CramWriter cramWriter;
	private PairedTemplateAssembler assembler;
	private Sam2CramRecordFactory cramRecordFactory;

	private Map<String, Integer> readGroupIdToIndexMap = new TreeMap<String, Integer>();
	private InputStream refSnpPosInputStream;
	private RefMaskUtils.RefMask refPileMasks;

	private String currentSequenceName;
	private LinkedList<SAMRecord> recordBuffer = new LinkedList<SAMRecord>();
	private byte[] refBases;

	public CRAMFileWriter(ReferenceSequenceFile rsf, OutputStream os,
			SAMFileHeader samHeader) throws IOException {
		this.referenceSequenceFile = rsf;
		this.os = os;
		this.samHeader = samHeader;

		init();
	}

	private void init() throws IOException {
		List<CramReferenceSequence> sequences;
		sequences = new ArrayList<CramReferenceSequence>();
		for (SAMSequenceRecord seq : samHeader.getSequenceDictionary()
				.getSequences()) {
			CramReferenceSequence cramSeq = new CramReferenceSequence(
					seq.getSequenceName(), seq.getSequenceLength());
			sequences.add(cramSeq);
		}

		List<CramReadGroup> cramReadGroups = new ArrayList<CramReadGroup>(
				samHeader.getReadGroups().size() + 1);
		cramReadGroups.add(new CramReadGroup(null));
		for (SAMReadGroupRecord rgr : samHeader.getReadGroups()) {
			readGroupIdToIndexMap.put(rgr.getReadGroupId(),
					readGroupIdToIndexMap.size() + 1);
			cramReadGroups.add(new CramReadGroup(rgr.getId(), rgr.getSample()));
		}

		cramWriter = new CramWriter(os, null, sequences, false, 100000, true,
				true, false, null, null, cramReadGroups);
		cramWriter.init();

		cramRecordFactory = new Sam2CramRecordFactory();
		cramRecordFactory.setCaptureInsertScores(true);
		cramRecordFactory.setCaptureSubtitutionScores(true);
		cramRecordFactory.setCaptureUnmappedScores(true);
		cramRecordFactory.setUncategorisedQualityScoreCutoff(0);
		cramRecordFactory.setTreatSoftClipsAs(TREAT_TYPE.IGNORE);
	}

	private static byte[] refPos2RefSNPs(InputStream is, int refLen,
			byte onValue) throws IOException {
		if (is == null)
			return null;
		byte[] refSNPs = new byte[refLen];
		BufferedReader r = new BufferedReader(new InputStreamReader(is));
		String line = null;
		while ((line = r.readLine()) != null) {
			int pos = Integer.valueOf(line);
			refSNPs[pos] = onValue;
		}
		r.close();
		return refSNPs;
	}

	private void flush() throws IOException, CramException {
		while (!recordBuffer.isEmpty())
			assembleTemplateAndWrite(recordBuffer.pollFirst());

		if (assembler != null) {
			SAMRecord assembledRecord = null;
			while ((assembledRecord = assembler.fetchNextSAMRecord()) != null) {
				writeSAMRecord(assembledRecord,
						assembler.distanceToNextFragment());
			}
		}
	}

	private void startSequence(String sequenceName) throws IOException,
			CramException {
		flush();

		assembler = new PairedTemplateAssembler();

		ReferenceSequence sequence = referenceSequenceFile
				.getSequence(sequenceName);
		refBases = referenceSequenceFile.getSubsequenceAt(sequence.getName(),
				1, sequence.length()).getBases();
		Utils.capitaliseAndCheckBases(refBases, false);

		cramWriter.startSequence(sequenceName, refBases);
		cramRecordFactory.setRefBases(refBases);

		byte[] refSNPs = refPos2RefSNPs(refSnpPosInputStream, refBases.length,
				(byte) '+');
		cramRecordFactory.setRefSNPs(refSNPs);

		refPileMasks = new RefMaskUtils.RefMask(refBases.length, 2);
		cramRecordFactory.setRefPile(refPileMasks);

		currentSequenceName = sequenceName;
	}

	private void writeSAMRecord(SAMRecord record, int distanceToNextFragment)
			throws IOException, CramException {
		CramRecord cramRecord = buildCramRecord(record);
		if (record.getReadGroup() != null) {
			String readGroupId = record.getReadGroup().getReadGroupId();
			Integer readGroupIndex = readGroupIdToIndexMap.get(readGroupId);
			cramRecord.setReadGroupID(readGroupIndex);
		}

		if (distanceToNextFragment > 0) {
			cramRecord.setLastFragment(false);
			cramRecord.setRecordsToNextFragment(distanceToNextFragment);
		} else
			cramRecord.setLastFragment(true);

		cramWriter.addRecord(cramRecord);
	}

	private CramRecord buildCramRecord(SAMRecord samRecord) {
		return cramRecordFactory.createCramRecord(samRecord);
	}

	private void assembleTemplateAndWrite(SAMRecord samRecord)
			throws IOException, CramException {
		assembler.addSAMRecord(samRecord);
		SAMRecord assembledRecord = null;
		while ((assembledRecord = assembler.nextSAMRecord()) != null) {
			writeSAMRecord(assembledRecord, assembler.distanceToNextFragment());
		}
	}

	private void addRefMask(SAMRecord record, byte[] refBases,
			RefMaskUtils.RefMask refMask) {
		int refStartInBlock;
		int readStartInBlock;
		int refStart;
		int readStart;
		byte readBase;
		byte refBase;
		for (AlignmentBlock block : record.getAlignmentBlocks()) {
			refStartInBlock = block.getReferenceStart();
			readStartInBlock = block.getReadStart();
			byte[] readBases = record.getReadBases();
			for (int i = 0; i < block.getLength(); i++) {
				refStart = refStartInBlock + i - 1;
				readStart = readStartInBlock + i - 1;
				readBase = readBases[readStart];
				refBase = refBases[refStart];
				refMask.addReadBase(refStart, readBase, refBase);
			}
		}
	}

	@Override
	public void addAlignment(SAMRecord samRecord) {
		try {
			if (currentSequenceName == null
					|| !currentSequenceName
							.equals(samRecord.getReferenceName())) {
				startSequence(samRecord.getReferenceName());
			}

			SAMRecord tempRecord = null;
			int alEnd = 0;

			recordBuffer.add(samRecord);

			if (refPileMasks != null)
				addRefMask(samRecord, refBases, refPileMasks);

			while (!recordBuffer.isEmpty()) {
				tempRecord = recordBuffer.peekFirst();
				if (tempRecord.getReadUnmappedFlag())
					alEnd = tempRecord.getAlignmentStart()
							+ tempRecord.getReadLength();
				else
					alEnd = tempRecord.getAlignmentEnd();

				if (alEnd < samRecord.getAlignmentStart()) {
					assembleTemplateAndWrite(recordBuffer.pollFirst());
				} else
					break;
			}

		} catch (IOException e) {
			throw new RuntimeException(e);
		} catch (CramException e) {
			throw new RuntimeException(e);
		}
	}

	@Override
	public SAMFileHeader getFileHeader() {
		return samHeader;
	}

	@Override
	public void close() {
		try {
			flush () ;
			cramWriter.close() ;
		} catch (IOException e) {
			throw new RuntimeException(e) ;
		} catch (CramException e) {
			throw new RuntimeException(e) ;
		}
	}

}

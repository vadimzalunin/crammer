package uk.ac.ebi.ena.sra.cram.spot;

import java.util.Map;
import java.util.Queue;
import java.util.TreeMap;
import java.util.concurrent.LinkedBlockingQueue;

import net.sf.samtools.SAMRecord;

public class PairedTemplateAssembler {
	public static final int POINTEE_DISTANCE_NOT_SET = -2;
	public static final int DISTANCE_NOT_SET = -1;

	private Queue<Record> recordQueue = new LinkedBlockingQueue<PairedTemplateAssembler.Record>();
	private Map<String, Record> recordsByNameMap = new TreeMap<String, PairedTemplateAssembler.Record>();

	private static class Record {
		SAMRecord samRecord;
		SAMRecord mateRecord;
		int nofRecordsToNextFragment = DISTANCE_NOT_SET;
		int index = -1;
		boolean needsAssembly = true;
	}

	private int counter = 0;
	private int alignmentHorizon;
	private int recordHorizon;

	private int distanceToNextFragment = DISTANCE_NOT_SET;
	private SAMRecord mateRecord;

	private String lastAddedReadSeqName;

	public PairedTemplateAssembler() {
		this(10000, 10000);
	}

	public PairedTemplateAssembler(int alignmentHorizon, int recordHorizon) {
		this.alignmentHorizon = alignmentHorizon;
		this.recordHorizon = recordHorizon;
	}

	public boolean isEmpty() {
		return recordQueue.isEmpty();
	}

	public void addSAMRecord(SAMRecord samRecord) {
		lastAddedReadSeqName = samRecord.getReferenceName();
		counter++;

		Record record = new Record();
		record.index = counter;
		record.samRecord = samRecord;
		recordQueue.add(record);

		if (!samRecord.getMateReferenceName().equals("=")
				&& !SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(samRecord.getMateReferenceName())
				&& !SAMRecord.NO_ALIGNMENT_REFERENCE_NAME.equals(samRecord.getReferenceName())
				&& !samRecord.getReferenceName().equals(samRecord.getMateReferenceName()))
			return;

		String spotName = samRecord.getReadName();

		if (recordsByNameMap.containsKey(spotName)) {
			Record foundRecord = recordsByNameMap.remove(spotName);

			foundRecord.nofRecordsToNextFragment = record.index - foundRecord.index;
			foundRecord.mateRecord = samRecord;
			record.nofRecordsToNextFragment = POINTEE_DISTANCE_NOT_SET;
			record.mateRecord = foundRecord.samRecord;
		} else
			recordsByNameMap.put(spotName, record);

	}

	public void addSAMRecordNoAssembly(SAMRecord samRecord) {
		lastAddedReadSeqName = samRecord.getReferenceName();
		counter++;

		Record record = new Record();
		record.needsAssembly = false;
		record.index = counter;
		record.samRecord = samRecord;
		recordQueue.add(record);
	}

	public SAMRecord nextSAMRecord() {
		Record record = recordQueue.peek();
		if (record == null)
			return null;
		if (!record.needsAssembly)
			return fetchNextSAMRecord();
		
		SAMRecord samRecord = record.samRecord;
		if (!samRecord.getReferenceName().equals(lastAddedReadSeqName))
			return fetchNextSAMRecord();

		if (record.nofRecordsToNextFragment == POINTEE_DISTANCE_NOT_SET)
			return fetchNextSAMRecord();

		if (!samRecord.getReadPairedFlag())
			return fetchNextSAMRecord();

		if (record.nofRecordsToNextFragment > -1)
			return fetchNextSAMRecord();

		if ((samRecord.getMateAlignmentStart() - samRecord.getAlignmentStart()) > alignmentHorizon)
			return fetchNextSAMRecord();

		if (recordQueue.size() > recordHorizon)
			return fetchNextSAMRecord();

		return null;
	}

	public SAMRecord fetchNextSAMRecord() {
		Record record = recordQueue.poll();
		if (record == null)
			return null;

		distanceToNextFragment = record.nofRecordsToNextFragment;
		mateRecord = record.mateRecord;
		Record remove = recordsByNameMap.remove(record.samRecord.getReadName());
		return record.samRecord;
	}

	public int distanceToNextFragment() {
		return distanceToNextFragment;
	}

	public SAMRecord getMateRecord() {
		return mateRecord;
	}

	public void clear() {
		counter = 0;
		distanceToNextFragment = -1;
		recordQueue.clear();
		recordsByNameMap.clear();
	}

	public void dumpHead() {
		System.out.println("PairedTemplateAssembler record queue size: " + recordQueue.size());
		System.out.println("PairedTemplateAssembler records by name map size: " + recordsByNameMap.size());
		Record head = recordQueue.peek();
		if (head == null)
			System.out.println("PairedTemplateAssembler head is null.");
		else {
			if (head.samRecord != null)
				System.out.println(head.samRecord.getSAMString());
			else
				System.out.println("PairedTemplateAssembler head: samRecord is null.");
			if (head.mateRecord != null)
				System.out.println(head.mateRecord.getSAMString());
			else
				System.out.println("PairedTemplateAssembler head: mateRecord is null.");
			System.out.println("distance=" + head.nofRecordsToNextFragment);

		}
	}
}

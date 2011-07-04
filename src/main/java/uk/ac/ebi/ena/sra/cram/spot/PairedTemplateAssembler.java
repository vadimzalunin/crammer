package uk.ac.ebi.ena.sra.cram.spot;

import java.util.Map;
import java.util.Queue;
import java.util.TreeMap;
import java.util.concurrent.LinkedBlockingQueue;

import net.sf.samtools.SAMRecord;

public class PairedTemplateAssembler {

	private Queue<Record> recordQueue = new LinkedBlockingQueue<PairedTemplateAssembler.Record>();
	private Map<String, Record> recordsByNameMap = new TreeMap<String, PairedTemplateAssembler.Record>();

	private static class Record {
		SAMRecord samRecord;
		int nofRecordsToNextFragment = -1;
		int index = -1;
	}

	private int counter = 0;
	private int alignmentHorizon;
	private int recordHorizon;

	private int distanceToNextFragment = -1;

	public PairedTemplateAssembler(int alignmentHorizon, int recordHorizon) {
		super();
		this.alignmentHorizon = alignmentHorizon;
		this.recordHorizon = recordHorizon;
	}

	public void addSAMRecord(SAMRecord samRecord) {
		counter++;

		Record record = new Record();
		record.index = counter;
		record.samRecord = samRecord;
		recordQueue.add(record);

		if (recordsByNameMap.containsKey(samRecord.getReadName())) {
			Record foundRecord = recordsByNameMap.get(samRecord.getReadName());

//			/*
//			 * The following can be simplified because there is no need to
//			 * compare the whole record for equality. We need to know if the
//			 * record is the same (which came from another reader) or is it the
//			 * mate pair with the same name.
//			 */
//			if (samRecord.equals(foundRecord.samRecord))
//				return;

			foundRecord.nofRecordsToNextFragment = record.index
					- foundRecord.index;
			recordsByNameMap.remove(samRecord.getReadName());
			record.nofRecordsToNextFragment = -2;
		} else
			recordsByNameMap.put(samRecord.getReadName(), record);

	}

	public SAMRecord nextSAMRecord() {
		Record record = recordQueue.peek();
		if (record == null)
			return null;
		SAMRecord samRecord = record.samRecord;

		if (!samRecord.getReadPairedFlag())
			return fetchNextSAMRecord();

		if (record.nofRecordsToNextFragment > -1)
			return fetchNextSAMRecord();

		if ((samRecord.getMateAlignmentStart() - samRecord.getAlignmentStart()) > alignmentHorizon)
			return fetchNextSAMRecord();

		if (samRecord.getMateAlignmentStart() < samRecord.getAlignmentStart())
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
		recordsByNameMap.remove(record.samRecord.getReadName());
		return record.samRecord;
	}

	public int distanceToNextFragment() {
		return distanceToNextFragment;
	}

	public void clear() {
		counter = 0;
		distanceToNextFragment = -1;
		recordQueue.clear();
		recordsByNameMap.clear();
	}
}

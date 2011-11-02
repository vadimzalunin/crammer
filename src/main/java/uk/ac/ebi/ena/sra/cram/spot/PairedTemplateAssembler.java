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
		SAMRecord mateRecord;
		int nofRecordsToNextFragment = -1;
		int index = -1;
	}

	private int counter = 0;
	private int alignmentHorizon;
	private int recordHorizon;

	private int distanceToNextFragment = -1;
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
		lastAddedReadSeqName = samRecord.getReferenceName() ;
		counter++;

		Record record = new Record();
		record.index = counter;
		record.samRecord = samRecord;
		recordQueue.add(record);

		String spotName = getSpotName(samRecord.getReadName());

		if (recordsByNameMap.containsKey(spotName)) {
			Record foundRecord = recordsByNameMap.remove(spotName);

			foundRecord.nofRecordsToNextFragment = record.index
					- foundRecord.index;
			foundRecord.mateRecord = samRecord;
			record.nofRecordsToNextFragment = -2;
			record.mateRecord = foundRecord.samRecord;
		} else
			recordsByNameMap.put(spotName, record);

	}

	private static final String getSpotName(String readName) {
		if (readName.endsWith(".1") || readName.endsWith(".2"))
			return readName.substring(0, readName.length() - 2);
		return readName;
	}

	public SAMRecord nextSAMRecord() {
		Record record = recordQueue.peek();
		if (record == null)
			return null;
		SAMRecord samRecord = record.samRecord;
		if (!samRecord.getReferenceName().equals(lastAddedReadSeqName))
			return fetchNextSAMRecord();

		if (record.nofRecordsToNextFragment == -2)
			return fetchNextSAMRecord();

		if (!samRecord.getReadPairedFlag())
			return fetchNextSAMRecord();

		if (record.nofRecordsToNextFragment > -1)
			return fetchNextSAMRecord();

		if ((samRecord.getMateAlignmentStart() - samRecord.getAlignmentStart()) > alignmentHorizon)
			return fetchNextSAMRecord();

		// if (samRecord.getMateAlignmentStart() <
		// samRecord.getAlignmentStart())
		// return fetchNextSAMRecord();

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
		Record remove = recordsByNameMap.remove(getSpotName(record.samRecord
				.getReadName()));
		// System.out
		// .printf("Fetched: distance=%d; read name=%s; remove.distance=%d; remove.read name=%s\n",
		// distanceToNextFragment, record.samRecord.getReadName(),
		// remove == null ? null : remove.nofRecordsToNextFragment,
		// remove == null ? null : remove.samRecord.getReadName());
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
}

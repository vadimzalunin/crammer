package uk.ac.ebi.ena.sra.cram.bam;

import java.io.File;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceRecord;

public class Test {

	public static void main(String[] args) throws InterruptedException {
		// File bamFile = new File("y:/Data/psyringae/psyringae.bam");
		File bamFile = new File("y:/Data/SangerExample/paired/5120_1.bam");
		SAMFileReader reader = new SAMFileReader(bamFile);
		List<SAMSequenceRecord> seqList = new LinkedList<SAMSequenceRecord>();
		for (SAMSequenceRecord sr : reader.getFileHeader()
				.getSequenceDictionary().getSequences()) {
			if (sr.getSequenceName().matches("^[\\d+]$")) {
				seqList.add(sr);
				System.out.println(sr.getSequenceName());
			}
		}

		SAMRecordIterator it = reader.query("1", 0, 100000, true);
		long counter_1 = 0;
		while (it.hasNext()) {
			it.next();
			counter_1++;
		}
		it.close();
		System.out.println("Count for 1: " + counter_1);

		ExecutorService es = Executors.newCachedThreadPool();

		BAMFileQueryQueues qq = new BAMFileQueryQueues(es, bamFile);
		qq.setDefaultQueueCapacity(100);

		List<BlockingQueue<SAMRecord>> list = new LinkedList<BlockingQueue<SAMRecord>>();

		for (SAMSequenceRecord sr : seqList) {
				list.addAll(qq.getQueuesForQuery(sr.getSequenceName(), 0,
						100000, true, 1));
		}

		long counter = 0;
		while (!list.isEmpty()) {
			Iterator<BlockingQueue<SAMRecord>> iterator = list.iterator();
			while (iterator.hasNext()) {
				if (iterator.next().take() == SAMRecordIteratorJob.STOP_SAMRECORD)
					iterator.remove();
				else
					counter++;
			}
		}
		System.out.println(counter);
		es.shutdownNow();
	}
}

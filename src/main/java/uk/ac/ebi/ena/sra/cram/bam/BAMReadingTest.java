package uk.ac.ebi.ena.sra.cram.bam;

import java.io.File;
import java.util.Queue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class BAMReadingTest {

	public static void main(String[] args) throws InterruptedException {
		File bamFile = new File(args[0]);
		String seqName = args.length > 1 ? args[1] : "1";
		int start = args.length > 2 ? Integer.valueOf(args[2]) : 0;
		int end = args.length > 3 ? Integer.valueOf(args[3]) : 0;
		boolean overlap = args.length > 4 ? Boolean.valueOf(args[4]) : false;

		System.out.printf("File=%s; seq=%s; start=%d; end=%d; overlap=%b\n",
				bamFile.getAbsolutePath(), seqName, start, end, overlap);

		long timeStart;
		long timeStop;
		ExecutorService es = Executors.newFixedThreadPool(1);

		for (int i = 0; i < 10; i++) {

			timeStart = System.currentTimeMillis();
			readWithQueue(es, bamFile, seqName, start, end, overlap);
			timeStop = System.currentTimeMillis();
			System.out.println("queue: " + (timeStop - timeStart));
			
			timeStart = System.currentTimeMillis();
			readSimple(bamFile, seqName, start, end, overlap);
			timeStop = System.currentTimeMillis();
			System.out.println("simple: " + (timeStop - timeStart));
		}
		
		es.shutdownNow();
	}

	private static void readSimple(File bamFile, String seqName, int start,
			int end, boolean overlap) {
		SAMFileReader reader = new SAMFileReader(bamFile);
		SAMRecordIterator iterator = reader.query(seqName, start, end, overlap);
		long alStart = 0L;
		while (iterator.hasNext()) {
			SAMRecord record = iterator.next();
			alStart += record.getAlignmentStart();
		}
		System.out.println(alStart);
	}

	private static void readWithQueue(ExecutorService es, File bamFile,
			String seqName, int start, int end, boolean overlap)
			throws InterruptedException {

		BAMFileQueryQueues qq = new BAMFileQueryQueues(es, bamFile);

		Queue<SAMRecord> queue = qq.getQueueForQuery(seqName,
				start, end, overlap);

		long alStart = 0L;
		SAMRecord record = null;
		while ((record = queue.poll()) != SAMRecordIteratorJob.STOP_SAMRECORD) {
			if (record == null) continue ;
			alStart += record.getAlignmentStart();
		}
		System.out.println(alStart);
		qq.stop(queue) ;
	}
}

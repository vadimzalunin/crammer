package uk.ac.ebi.ena.sra.cram.bam;

import java.util.Queue;
import java.util.concurrent.atomic.AtomicBoolean;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class SAMRecordIteratorJob implements Runnable {
	public static SAMRecord STOP_SAMRECORD = new SAMRecord(null);
	private static final long DEFAULT_SLEEP_TIME = 100;
	private SAMRecordIterator iterator;
	private Queue<SAMRecord> queue;
	private AtomicBoolean abort = new AtomicBoolean(false);
	private AtomicBoolean finished = new AtomicBoolean(false);
	private long sleepTime;
	private Throwable exception;

	public SAMRecordIteratorJob(SAMRecordIterator iterator,
			Queue<SAMRecord> queue, long sleepTime) {
		this.iterator = iterator;
		this.queue = queue;
		this.sleepTime = sleepTime;
	}

	public SAMRecordIteratorJob(SAMRecordIterator iterator,
			Queue<SAMRecord> queue) {
		this(iterator, queue, DEFAULT_SLEEP_TIME);
	}

	@Override
	public void run() {
		long counter = 0;
		long nullCounter = 0;

		try {
			while (!abort.get() && iterator.hasNext()) {
				SAMRecord record = iterator.next();
				counter++;
				while (!abort.get() && !queue.offer(record)) {
					nullCounter++;
					// try {
					// Thread.sleep(sleepTime);
					// } catch (InterruptedException e) {
					// return;
					// }
				}
			}
			queue.add(STOP_SAMRECORD);
		} catch (Throwable e) {
			e.printStackTrace();
			exception = e;
		} finally {
			System.err.println("SAMRecordIteratorJob exiting.");
			System.err.println("nullCounter=" + nullCounter);
			finished.set(true);
			try {
				iterator.close();
			} catch (Throwable t) {
				t.printStackTrace();
			}
		}
		System.out.println("SAMRecordIteratorJob: " + counter);
	}

	public void abort() {
		abort.set(true);
	}

	public boolean isFinished() {
		return finished.get();
	}

	public Throwable getException() {
		return exception;
	}
}

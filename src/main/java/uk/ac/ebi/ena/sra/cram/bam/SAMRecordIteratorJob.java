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
package uk.ac.ebi.ena.sra.cram.bam;

import java.util.Queue;
import java.util.concurrent.atomic.AtomicBoolean;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class SAMRecordIteratorJob implements Runnable {
	public final static SAMRecord STOP_SAMRECORD = new SAMRecord(null);
	private static final long DEFAULT_SLEEP_TIME = 100;
	private SAMRecordIterator iterator;
	private Queue<SAMRecord> queue;
	private AtomicBoolean abort = new AtomicBoolean(false);
	private AtomicBoolean finished = new AtomicBoolean(false);
	private long sleepTime;
	private Throwable exception;

	public SAMRecordIteratorJob(SAMRecordIterator iterator, Queue<SAMRecord> queue, long sleepTime) {
		this.iterator = iterator;
		this.queue = queue;
		this.sleepTime = sleepTime;
	}

	public SAMRecordIteratorJob(SAMRecordIterator iterator, Queue<SAMRecord> queue) {
		this(iterator, queue, DEFAULT_SLEEP_TIME);
	}

	@Override
	public void run() {
		long tries = 0;

		try {
			while (!abort.get() && iterator.hasNext()) {
				SAMRecord record = iterator.next();
				tries = 0;
				while (!abort.get() && !queue.offer(record)) {
					Thread.sleep(sleepTime);
					tries++ ;
				}
			}

			tries = 0;
			while (!abort.get() && !queue.offer(STOP_SAMRECORD)) 
				Thread.sleep(sleepTime);
			
		} catch (Throwable e) {
			e.printStackTrace();
			exception = e;
		} finally {
			finished.set(true);
//			System.out.println("SAMRecordIteratorJob: tries=" + tries);
		}
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

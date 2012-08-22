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

import java.io.File;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.LinkedBlockingQueue;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class BAMFileQueryQueues {
	private ExecutorService es;
	private File bamFile;
	private int defaultQueueCapacity = 100000;

	private Map<Queue<SAMRecord>, SAMRecordIteratorJob> queues = new HashMap<Queue<SAMRecord>, SAMRecordIteratorJob>();

	public void stop(Queue<SAMRecord> queue) {
		SAMRecordIteratorJob job = queues.get(queue);
		if (job != null) {
			job.abort();
			queues.remove(queue);
		}
	}

	public BAMFileQueryQueues(ExecutorService es, File bamFile) {
		super();
		this.es = es;
		this.bamFile = bamFile;
	}

	public Collection<BlockingQueue<SAMRecord>> getQueuesForQuery(String sequenceName, int start, int end,
			boolean overlaping, int nofQueues) {
		List<BlockingQueue<SAMRecord>> list = new LinkedList<BlockingQueue<SAMRecord>>();
		int step = (int) Math.ceil(((float) end - start) / nofQueues);
		for (int i = 0; i < nofQueues; i++) {
			BlockingQueue<SAMRecord> q = getBlockingQueueForQuery(sequenceName, start + (step * i), start
					+ (step * (i + 1)), overlaping);
			list.add(q);
		}

		return list;
	}

	public BlockingQueue<SAMRecord> getBlockingQueueForQuery(String sequenceName, int start, int end,
			boolean overlapping) {
		SAMFileReader reader = new SAMFileReader(bamFile);
		SAMRecordIterator iterator = reader.query(sequenceName, start, end, overlapping);
		return getBlockingQueue(iterator);
	}

	public Queue<SAMRecord> getQueueForQuery(String sequenceName, int start, int end, boolean overlapping) {
		SAMFileReader reader = new SAMFileReader(bamFile);
		SAMRecordIterator iterator = reader.query(sequenceName, start, end, overlapping);
		Queue<SAMRecord> queue = new ConcurrentLinkedQueue<SAMRecord>();
		startQuery(iterator, queue);
		return queue;
	}

	public BlockingQueue<SAMRecord> getBlockingQueue(SAMRecordIterator iterator) {
		BlockingQueue<SAMRecord> queue = new LinkedBlockingQueue<SAMRecord>(defaultQueueCapacity);

		startQuery(iterator, queue);

		return queue;
	}

	public void startQuery(SAMRecordIterator iterator, Queue<SAMRecord> queue) {
		SAMRecordIteratorJob queuer = new SAMRecordIteratorJob(iterator, queue);

		es.submit(queuer);
	}

	public int getDefaultQueueCapacity() {
		return defaultQueueCapacity;
	}

	public void setDefaultQueueCapacity(int defaultQueueCapacity) {
		this.defaultQueueCapacity = defaultQueueCapacity;
	}
}

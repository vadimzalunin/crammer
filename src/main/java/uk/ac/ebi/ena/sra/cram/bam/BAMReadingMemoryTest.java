package uk.ac.ebi.ena.sra.cram.bam;

import java.io.File;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;
import java.util.concurrent.atomic.AtomicBoolean;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceRecord;

public class BAMReadingMemoryTest {

	public static void main(String[] args) {
		File bamFile = new File(args[0]);

		final Map<String, SAMRecord> recordMap = Collections
				.synchronizedMap(new TreeMap<String, SAMRecord>());
		final AtomicBoolean stop = new AtomicBoolean(false);
		final long sleepMillis = 1000;
		new Thread(new Runnable() {

			@Override
			public void run() {
				while (!stop.get()) {
					System.out.printf(
							"records=%d k;\tfree memory=%d MB;\ttime=%d\n",
							recordMap.size() / 1000, Runtime.getRuntime()
									.freeMemory() / (1024 * 1024),
							System.currentTimeMillis());
					try {
						Thread.sleep(sleepMillis);
					} catch (InterruptedException e) {
						e.printStackTrace();
						return;
					}
				}
			}
		}).start();

		try {
			SAMFileReader
					.setDefaultValidationStringency(ValidationStringency.SILENT);
			SAMFileReader reader = new SAMFileReader(bamFile);
			List<String> seqNames = new LinkedList<String>();
			for (SAMSequenceRecord seq : reader.getFileHeader()
					.getSequenceDictionary().getSequences())
				seqNames.add(seq.getSequenceName());

			for (String seqName : seqNames) {
				SAMRecordIterator iterator = reader.query(seqName, 0, 0, true);
				try {
					while (iterator.hasNext()) {
						SAMRecord record = iterator.next();
						recordMap.put(record.getReadName(), record);
					}
				} finally {
					try {
						iterator.close();
					} catch (Exception e) {
						e.printStackTrace();
						break;
					}
				}
			}
		} finally {
			stop.set(true);
		}

	}

	private static void dump(long size) {
		System.out.println();
	}
}

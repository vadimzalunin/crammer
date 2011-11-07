package uk.ac.ebi.ena.sra.cram.impl;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertThat;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;

import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.SeekableFileStream;

import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;

import uk.ac.ebi.ena.sra.cram.CramIndexer;
import uk.ac.ebi.ena.sra.cram.CramTools;
import uk.ac.ebi.ena.sra.cram.format.CramFormatException;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.compression.CramCompressionException;
import uk.ac.ebi.ena.sra.cram.index.CramIndex;
import uk.ac.ebi.ena.sra.cram.index.RecordPointer;

public class TestCramIterators {
	private static long nofRecords = 100L;
	private static File bamFile = new File("data/set5/input.bam");
	private static File refFile = new File("data/set5/ref.fa");
	private static File cramFile;
	private static File indexFile;
	private static CramIndex index ;
	private static RecordPointer pointer;

	@BeforeClass
	public static void createCramFile() throws Exception {
		cramFile = File.createTempFile(bamFile.getName(), ".cram");
		cramFile.deleteOnExit();
		String command = String
				.format("-l error cram --input-bam-file %s --reference-fasta-file %s --output-cram-file %s --max-records 100",
						bamFile.getAbsolutePath(), refFile.getAbsolutePath(),
						cramFile.getAbsolutePath());
		CramTools.main(command.split("\\s+"));

		indexFile = File.createTempFile(cramFile.getName(), ".crai");
		indexFile.deleteOnExit();
		command = String
				.format("--input-cram-file %s --reference-fasta-file %s --index-file %s",
						cramFile.getAbsolutePath(), refFile.getAbsolutePath(),
						indexFile.getAbsolutePath());
		CramIndexer.main(command.split("\\s+"));
		
		index = CramIndex.fromFile(indexFile);
		pointer = index.findRecordPointerAt("20", 0);
	}

	public void createIndex() {
	}

	private void doTest(CloseableIterator<CramRecord> iterator, long size) {
		long counter = 0;
		while (iterator.hasNext()) {
			CramRecord record = iterator.next();
			counter++;
		}

		assertThat(counter, is(size));
	}

	@Test
	public void test1() throws FileNotFoundException, IOException {
		CloseableIterator<CramRecord> iterator = new CramIterator(
				new FileInputStream(cramFile),
				ReferenceSequenceFileFactory.getReferenceSequenceFile(refFile));

		doTest(iterator, nofRecords);
	}

	@Test
	public void test2() throws FileNotFoundException, IOException,
			CramFormatException, CramCompressionException {
		CloseableIterator<CramRecord> iterator = new CRAMPreemptiveIterator(
				new FileInputStream(cramFile),
				ReferenceSequenceFileFactory.getReferenceSequenceFile(refFile),
				null);

		doTest(iterator, nofRecords);
	}

	@Test
	public void test3() throws FileNotFoundException, IOException,
			CramFormatException, CramCompressionException {
		CloseableIterator<CramRecord> iterator = new CramRandomAccessIterator(
				new SeekableFileStream(cramFile),
				ReferenceSequenceFileFactory.getReferenceSequenceFile(refFile),
				pointer);

		doTest(iterator, nofRecords);
	}
	
	@Test
	public void test4() throws FileNotFoundException, IOException,
			CramFormatException, CramCompressionException {
		CloseableIterator<CramRecord> iterator = new CramPreemptiveRandomAccessIterator(
				new SeekableFileStream(cramFile),
				ReferenceSequenceFileFactory.getReferenceSequenceFile(refFile),
				pointer);

		doTest(iterator, nofRecords);
	}

}

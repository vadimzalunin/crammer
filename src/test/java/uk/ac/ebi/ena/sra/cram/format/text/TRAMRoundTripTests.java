package uk.ac.ebi.ena.sra.cram.format.text;

import static org.hamcrest.CoreMatchers.*;
import static org.junit.Assert.assertThat;

import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.ByteBuffer;
import java.util.Arrays;
import java.util.Collection;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;
import org.xerial.snappy.Snappy;

import uk.ac.ebi.ena.sra.cram.CramTools;
import uk.ac.ebi.ena.sra.cram.Utils;
import uk.ac.ebi.ena.sra.cram.format.CramHeader;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.CramRecordBlock;
import uk.ac.ebi.ena.sra.cram.impl.ByteArraySequenceBaseProvider;
import uk.ac.ebi.ena.sra.cram.impl.CramHeaderIO;
import uk.ac.ebi.ena.sra.cram.impl.RecordCodecFactory;
import uk.ac.ebi.ena.sra.cram.impl.RestoreBases;
import uk.ac.ebi.ena.sra.cram.impl.RestoreQualityScores;
import uk.ac.ebi.ena.sra.cram.impl.SequentialCramReader;

@RunWith(value = Parameterized.class)
public class TRAMRoundTripTests {
	private File folder;

	private static File findDataSetFolder(String datasetName) {
		File folder = new File(String.format("data/%s/", datasetName));
		if (folder.exists())
			return folder;

		return new File(String.format("../data/%s/", datasetName));
	}

	public static void main(String[] args) throws Exception {
		String sourceDatasetName = "set1";
		String destDatasetName = "set4";

		File sourceFolder = findDataSetFolder(sourceDatasetName);
		if (!sourceFolder.exists() || !sourceFolder.isDirectory())
			throw new RuntimeException("Dataset folder not found: "
					+ sourceDatasetName);

		File destFolder = findDataSetFolder(destDatasetName);
		if (!destFolder.exists() || !destFolder.isDirectory())
			throw new RuntimeException("Dataset folder not found: "
					+ destDatasetName);

		File inputBamPath = new File(sourceFolder, "input.bam");
		File refPath = new File(sourceFolder, "ref.fa");
		File readQualityMaskFile = new File(destFolder, "readmask.pos");

		ReferenceSequenceFile referenceSequenceFile = ReferenceSequenceFileFactory
				.getReferenceSequenceFile(refPath);

		// generate CRAM from BAM:
		File cramFile = new File(destFolder, "input.cram");
		// File cramFile = File.createTempFile(destDatasetName, ".cram");
		// cramFile.deleteOnExit();

		String cmd1 = String
				.format("-l ERROR cram --print-cram-records --max-records 10000 --input-bam-file %s --reference-fasta-file %s --output-cram-file %s --read-quality-mask-file %s --capture-unmapped-quality-scores --capture-substitution-quality-scores --capture-masked-quality-scores",
						inputBamPath, refPath, cramFile.getAbsolutePath(),
						readQualityMaskFile.getAbsolutePath());

		CramTools.main(cmd1.split(" "));

		File textCramFile = new File(destFolder, "input.tram");
		FileWriter textCramWriter = new FileWriter(textCramFile);

		DataInputStream cramDIS = new DataInputStream(new FileInputStream(
				cramFile));
		CramHeader cramHeader = CramHeaderIO.read(cramDIS);
		ReferenceSequence nextSequence = referenceSequenceFile
				.getSequence(cramHeader.getReferenceSequences().iterator()
						.next().getName());
		byte[] refBases = referenceSequenceFile.getSubsequenceAt(
				nextSequence.getName(), 1, nextSequence.length()).getBases();
		Utils.capitaliseAndCheckBases(refBases, true);
		ByteArraySequenceBaseProvider provider = new ByteArraySequenceBaseProvider(
				refBases);
		SequentialCramReader reader = new SequentialCramReader(cramDIS,
				provider);

		System.out.println();
		CramRecordFormat format = new CramRecordFormat();
		CramRecordBlock block = null;
		while ((block = reader.readBlock()) != null) {

			for (int i = 0; i < block.getRecordCount(); i++) {
				CramRecord record = reader.readRecord();
				
				String s = format.writeRecord(record);
				CramRecord derivedRecord = format.fromString(s);
				if (!record.equals(derivedRecord)) {
					System.err.println(s);
					System.err.println(record.toString());
					System.err.println(derivedRecord.toString());
					throw new RuntimeException(s);
				}
				textCramWriter.write(s);
				textCramWriter.write('\n');
			}
		}

		textCramWriter.close();

		// System.out.printf("BAM file size=%d; snapped=%d\n",
		// inputBamPath.length(), getSnappedSize(inputBamPath));
		// System.out.printf("CRAM file size=%d; snapped=%d\n",
		// cramFile.length(),
		// getSnappedSize(cramFile));
		// System.out.printf("TRAM file size=%d; snapped=%d\n",
		// textCramFile.length(), getSnappedSize(textCramFile));
	}

	// private static byte[] getContent(File file) throws IOException {
	// if (file.length() >= Integer.MAX_VALUE)
	// throw new RuntimeException("File size is too big: " + file.length());
	//
	// BufferedInputStream bis = new BufferedInputStream(new FileInputStream(
	// file));
	//
	// ByteBuffer byteBuffer = ByteBuffer.allocate((int) file.length());
	// byte[] readBuf = new byte[1024 * 4];
	// int len = 0;
	// while ((len = bis.read(readBuf)) != -1) {
	// byteBuffer.put(readBuf, 0, len);
	// }
	// bis.close();
	//
	// byteBuffer.flip();
	// byte[] data = new byte[byteBuffer.limit()];
	// byteBuffer.get(data);
	//
	// return data;
	// }
	//
	// private static long getSnappedSize(File file) throws IOException {
	// return Snappy.compress(getContent(file)).length;
	// }

	public TRAMRoundTripTests(String datasetName) {
		folder = findDataSetFolder(datasetName);
	}

	@Parameters
	public static Collection<Object[]> data() {
		Object[][] data = new Object[][] { { "set4" }, };
		return Arrays.asList(data);
	}

	@Test
	public void test1() throws Exception {
		File tramFile = new File(folder, "input.tram");
		File refFile = new File(folder, "ref.fa");
		File bamFile = new File(folder, "input.bam");

		SAMFileReader
				.setDefaultValidationStringency(ValidationStringency.SILENT);
		SAMFileReader samFileReader = new SAMFileReader(bamFile);
		SAMRecordIterator iterator = samFileReader.iterator();

		ReferenceSequenceFile referenceSequenceFile = ReferenceSequenceFileFactory
				.getReferenceSequenceFile(refFile);
		byte[] refBases = referenceSequenceFile.getSequence("1").getBases();
		ByteArraySequenceBaseProvider provider = new ByteArraySequenceBaseProvider(
				refBases);
		RestoreBases restoreBases = new RestoreBases(provider, "1");
		RestoreQualityScores restoreQualityScores = new RestoreQualityScores();

		CramRecordFormat format = new CramRecordFormat();

		BufferedReader reader = new BufferedReader(new FileReader(tramFile));
		String line = null;
		int readCounter = 1;
		while ((line = reader.readLine()) != null) {
			CramRecord record = format.fromString(line);

			String derivedString = format.writeRecord(record);

			assertThat(derivedString, notNullValue());
			assertThat(derivedString, equalTo(line));
			
			SAMRecord samRecord = iterator.next();
			byte[] restoredReadBases;
			if (record.isReadMapped())
				restoredReadBases = restoreBases.restoreReadBases(record);
			else
				restoredReadBases = record.getReadBases();


			assertThat("Mismatch for record: " + samRecord.getReadName()
					+ " for read number " + readCounter, new String(
					restoredReadBases),
					equalTo(new String(samRecord.getReadBases())));

			byte[] restoredScore = restoreQualityScores
					.restoreQualityScores(record);
			byte[] bamScores = samRecord.getBaseQualities();
			for (int i = 0; i < restoredScore.length; i++) {
				if (restoredScore[i] != 32)
					assertThat(restoredScore[i], is((byte) (bamScores[i] + 33)));
			}

			readCounter++;
		}
	}
}

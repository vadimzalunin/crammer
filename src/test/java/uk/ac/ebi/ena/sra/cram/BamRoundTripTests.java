package uk.ac.ebi.ena.sra.cram;

import java.io.File;
import java.util.Arrays;
import java.util.Collection;

import net.sf.picard.sam.BuildBamIndex;
import net.sf.samtools.SAMFileReader;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

@RunWith(value = Parameterized.class)
public class BamRoundTripTests {
	private String datasetName;

	private String inputBamPath;
	private String refPath;

	public BamRoundTripTests(String datasetName) {
		this.datasetName = datasetName;

		// let's find the dataset. It could be in 'data/setX' or in
		// '../data/setX':
		String prefix = String.format("data/%s/", datasetName);
		inputBamPath = String.format("%s/input.bam", prefix);
		if (!new File(inputBamPath).exists()) {
			prefix = String.format("../data/%s/", datasetName);
			inputBamPath = String.format("%s/input.bam", prefix);
		}
		refPath = String.format("%s/ref.fa", prefix);
	}

	@Parameters
	public static Collection<Object[]> data() {
		Object[][] data = new Object[][] { { "set1" }, { "set2" }, { "set3" },

		};
		return Arrays.asList(data);
	}

	@Test
	public void test1() throws Exception {
		File cramFile = File.createTempFile(datasetName, ".cram");
		cramFile.deleteOnExit();

		String cmd1 = String
				.format("-l INFO cram --input-bam-file %s --reference-fasta-file %s --output-cram-file %s",
						inputBamPath, refPath, cramFile.getAbsolutePath());
		System.out.println(cmd1);

		CramTools.main(cmd1.split(" "));

		File bamFile = File.createTempFile(datasetName, ".bam");
		bamFile.deleteOnExit();

		String cmd2 = String
				.format("bam --input-cram-file %s --reference-fasta-file %s --output-bam-file %s",
						cramFile.getAbsolutePath(), refPath, bamFile);
		System.out.println(cmd2);

		CramTools.main(cmd2.split(" "));

		File indexFile = new File(bamFile.getParentFile(), bamFile.getName()
				+ ".bai");
		indexFile.deleteOnExit();

		SAMFileReader reader = new SAMFileReader(bamFile);
		BuildBamIndex.createIndex(reader, indexFile);
		reader.close();

		String cmd3 = String
				.format("-l INFO cram --input-bam-file %s --reference-fasta-file %s --output-cram-file %s",
						bamFile, refPath, cramFile.getAbsolutePath());
		System.out.println(cmd3);

		CramTools.main(cmd3.split(" "));

	}
}

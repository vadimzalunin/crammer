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
package uk.ac.ebi.ena.sra.cram;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertThat;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import net.sf.picard.sam.BuildBamIndex;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

@RunWith(value = Parameterized.class)
public class BamRoundTripTests {
	private String datasetName;

	private String inputBamPath;
	private String refPath;
	private String model;

	public BamRoundTripTests(String datasetName, String ref, String model) {
		this.datasetName = datasetName;

		if (ref == null) {
			// let's find the dataset. It could be in 'data/setX' or in
			// '../data/setX':
			String prefix = String.format("data/%s/", datasetName);
			inputBamPath = String.format("%s/input.bam", prefix);
			if (!new File(inputBamPath).exists()) {
				prefix = String.format("../data/%s/", datasetName);
				inputBamPath = String.format("%s/input.bam", prefix);
			}
			refPath = String.format("%s/ref.fa", prefix);
		} else
			refPath = ref;

		this.model = model;
	}

	@Parameters
	public static Collection<Object[]> data() {
		String[] datasets = new String[] { "set1", "set2", "set3", "set4" };

		List<String> modelList = new ArrayList<String>();
		modelList.add("--capture-substitution-quality-scores --capture-insertion-quality-scores");
		modelList.add("--capture-all-quality-scores --include-unmapped-reads --preserve-read-names --capture-all-tags");

		List<Object[]> data = new ArrayList<Object[]>();

		for (String dataset : datasets)
			for (String model : modelList)
				data.add(new Object[] { dataset, null, model });

		return data;
	}

	@Test
	public void test1() throws Exception {
		File cramFileGeneration1 = File.createTempFile(datasetName, ".cram");
		cramFileGeneration1.deleteOnExit();

//		String model = "--capture-all-quality-scores --include-unmapped-reads --preserve-read-names --capture-all-tags";
		// String model =
		// "--capture-substitution-quality-scores --capture-insertion-quality-scores"
		// ;

		String cmd1 = String.format(
				"-l INFO cram --input-bam-file %s --reference-fasta-file %s --output-cram-file %s %s", inputBamPath,
				refPath, cramFileGeneration1.getAbsolutePath(), model);
		// System.out.println(cmd1);

		CramTools.main(cmd1.split(" "));

		File bamFile = File.createTempFile(datasetName, ".bam");
		bamFile.deleteOnExit();
		// System.out.println(bamFile.getAbsolutePath());

		String cmd2 = String.format("-l INFO bam --input-cram-file %s --reference-fasta-file %s --output-bam-file %s",
				cramFileGeneration1.getAbsolutePath(), refPath, bamFile);
		// System.out.println(cmd2);

		CramTools.main(cmd2.split(" "));

		File indexFile = new File(bamFile.getParentFile(), bamFile.getName() + ".bai");
		indexFile.deleteOnExit();

		SAMFileReader reader = new SAMFileReader(bamFile);
		reader.setValidationStringency(ValidationStringency.SILENT);
		BuildBamIndex.createIndex(reader, indexFile);
		reader.close();

		File cramFileGeneration2 = File.createTempFile(datasetName, ".cram");
		cramFileGeneration2.deleteOnExit();

		String cmd3 = String.format(
				"-l INFO cram --input-bam-file %s --reference-fasta-file %s --output-cram-file %s %s", bamFile,
				refPath, cramFileGeneration2.getAbsolutePath(), model);
		// System.out.println(cmd3);

		CramTools.main(cmd3.split(" "));

		File bamFileGeneration2 = File.createTempFile(datasetName, ".bam");
		bamFileGeneration2.deleteOnExit();
		// System.out.println(bamFileGeneration2.getAbsolutePath());

		String cmd4 = String.format("-l ERROR bam --input-cram-file %s --reference-fasta-file %s --output-bam-file %s",
				cramFileGeneration2.getAbsolutePath(), refPath, bamFileGeneration2);
		// System.out.println(cmd4);

		CramTools.main(cmd4.split(" "));

		assertThat(isContentSame(bamFile, bamFileGeneration2), is(true));
	}

	private static boolean isContentSame(File... files) throws NoSuchAlgorithmException, IOException {
		if (files == null)
			throw new NullPointerException();
		if (files.length < 2)
			throw new IllegalAccessError("Expecting at least 2 files.");

		byte[] hash = null;
		for (File file : files) {
			if (hash == null)
				hash = getHash(file);
			else if (!Arrays.equals(hash, getHash(file)))
				return false;
		}
		return true;
	}

	private static byte[] getHash(File file) throws NoSuchAlgorithmException, IOException {
		if (!file.isFile())
			throw new IllegalArgumentException("Not a file: " + file.getAbsolutePath());

		MessageDigest digest = MessageDigest.getInstance("SHA-1");
		digest.reset();

		FileInputStream fis = new FileInputStream(file);
		byte[] buf = new byte[1024];
		int len = 0;
		while ((len = fis.read(buf)) != -1) {
			digest.update(buf, 0, len);
		}

		byte[] hash = digest.digest();
		// System.out.println(Arrays.toString(hash));
		return hash;
	}

}

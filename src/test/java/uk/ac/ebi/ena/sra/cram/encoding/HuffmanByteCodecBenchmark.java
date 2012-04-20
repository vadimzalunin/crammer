package uk.ac.ebi.ena.sra.cram.encoding;

import static org.hamcrest.CoreMatchers.equalTo;
import static org.junit.Assert.assertThat;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Scanner;

import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;

import uk.ac.ebi.ena.sra.compression.huffman.HuffmanCode;
import uk.ac.ebi.ena.sra.compression.huffman.HuffmanTree;
import uk.ac.ebi.ena.sra.cram.Utils;
import uk.ac.ebi.ena.sra.cram.format.ByteFrequencies;
import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitInputStream;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitOutputStream;

public class HuffmanByteCodecBenchmark {
	private File dir = new File("data/set8/");
	private File qsFile = new File(dir, "qs");
	private HuffmanByteCodec2 codec;
	private List<String> lines;

	@Before
	public void before() throws IOException {
		ByteFrequencies bf = read(new File(dir, "freqs"));

		HuffmanTree<Byte> tree = HuffmanCode.buildTree(bf.getFrequencies(), Utils.autobox(bf.getValues()));
		codec = new HuffmanByteCodec2(tree);
		

		Scanner scanner = new Scanner(qsFile);

		lines = new ArrayList<String>(1000);
		while (scanner.hasNextLine()) {
			String line = scanner.nextLine();
			lines.add(line);
		}
	}

	@Test
	@Ignore
	public void test1() throws FileNotFoundException {
		ByteFrequencies f = new ByteFrequencies();
		Scanner scanner = new Scanner(qsFile);

		while (scanner.hasNextLine()) {
			String line = scanner.nextLine();

			for (byte b : line.getBytes())
				f.add(b);
		}

		PrintStream ps = new PrintStream(new File(dir, "freqs"));
		ps.println(Arrays.toString(f.getValues()));
		ps.println(Arrays.toString(f.getFrequencies()));
	}

	@Test
	public void test2() throws IOException {
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);
		
		for (String line:lines) {
			for (byte b : line.getBytes())
				codec.write(bos, b);
		}

		bos.close();
		byte[] buf = baos.toByteArray() ;

		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		BitInputStream bis = new DefaultBitInputStream(bais);

		for (int i = 0; i < lines.size(); i++) {
			String line = lines.get(i);
			byte[] bytes = new byte[line.length()];
			for (int b = 0; b < bytes.length; b++) {
				bytes[b] = codec.read(bis);
			}

			assertThat(bytes, equalTo(line.getBytes()));
		}
	}

	private ByteFrequencies read(File file) throws FileNotFoundException {
		Scanner scanner = new Scanner(file);
		String valueLine = scanner.nextLine();
		String freqLine = scanner.nextLine();
		scanner.close();

		valueLine = valueLine.substring(1, valueLine.length() - 1);
		freqLine = freqLine.substring(1, freqLine.length() - 1);

		String[] valueLineSplit = valueLine.split(", ");
		String[] freqLineSplit = freqLine.split(", ");

		ByteFrequencies bf = new ByteFrequencies();
		for (int i = 0; i < valueLineSplit.length; i++) {
			byte byteValue = Byte.valueOf(valueLineSplit[i]);
			int intFreq = Integer.valueOf(freqLineSplit[i]);
			bf.add(byteValue, intFreq);
		}

		return bf;
	}

}

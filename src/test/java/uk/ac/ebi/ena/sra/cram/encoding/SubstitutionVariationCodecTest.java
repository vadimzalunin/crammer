package uk.ac.ebi.ena.sra.cram.encoding;

import static org.hamcrest.CoreMatchers.equalTo;
import static org.junit.Assert.assertThat;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;

import org.junit.Test;

import uk.ac.ebi.ena.sra.compression.huffman.HuffmanCode;
import uk.ac.ebi.ena.sra.compression.huffman.HuffmanTree;
import uk.ac.ebi.ena.sra.cram.format.SubstitutionVariation;
import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitInputStream;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitOutputStream;

public class SubstitutionVariationCodecTest {

	@Test
	public void test1() throws IOException {
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);

		SubstitutionVariationCodec codec = new SubstitutionVariationCodec();
		codec.baseChangeCodec = new BaseChangeCodec();
		HuffmanTree<Byte> qualityScoreTree = HuffmanCode.buildTree(new int[] {
				100, 90, 80 }, new Byte[] { 33, 34, 35 });

		SubstitutionVariation v = new SubstitutionVariation();
		v.setBaseChange(new BaseChange((byte)'A', (byte)'C')) ;
		v.setPosition(-1);
		codec.write(bos, v);

		bos.flush();

		ByteArrayInputStream bais = new ByteArrayInputStream(baos.toByteArray());
		BitInputStream bis = new DefaultBitInputStream(bais);
		SubstitutionVariation v2 = codec.read(bis);

		assertThat(v2, equalTo(v));

	}
}

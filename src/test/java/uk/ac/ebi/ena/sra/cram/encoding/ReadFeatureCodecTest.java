package uk.ac.ebi.ena.sra.cram.encoding;

import static org.hamcrest.CoreMatchers.*;
import static org.junit.Assert.assertThat;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import org.junit.Test;

import uk.ac.ebi.ena.sra.compression.huffman.HuffmanCode;
import uk.ac.ebi.ena.sra.compression.huffman.HuffmanTree;
import uk.ac.ebi.ena.sra.cram.Utils;
import uk.ac.ebi.ena.sra.cram.format.ReadFeature;
import uk.ac.ebi.ena.sra.cram.format.SubstitutionVariation;
import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitInputStream;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitOutputStream;

public class ReadFeatureCodecTest {
	@Test
	public void test1() throws IOException {
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		BitOutputStream bos = new DefaultBitOutputStream(baos);

		SubstitutionVariationCodec codec = new SubstitutionVariationCodec();
		codec.baseChangeCodec = new BaseChangeCodec();
		HuffmanTree<Byte> qualityScoreTree = HuffmanCode.buildTree(new int[] {
				100, 90, 80 }, new Byte[] { 33, 34, 35 });
		HuffmanCodec<Byte> qualityScoreCodec = new HuffmanCodec<Byte>(
				qualityScoreTree);
		codec.qualityScoreCodec = qualityScoreCodec;

		ReadFeatureCodec featuresCodec = new ReadFeatureCodec();
		featuresCodec.substitutionCodec = codec;
		featuresCodec.featureOperationCodec = new ReadFeatureOperatorCodec(
				new int[] { 100, 90, 10, 10, 1 }, Utils.autobox("S$IDN"
						.getBytes()));
		featuresCodec.inReadPosCodec = new GolombRiceCodec(1);

		SubstitutionVariation v = new SubstitutionVariation();
		v.setBaseChange(new BaseChange((byte) 'C', (byte) 'A'));
		v.setQualityScore((byte) '!');
		v.setPosition(20);
		SubstitutionVariation v2 = new SubstitutionVariation();
		v2.setBaseChange(new BaseChange((byte) 'C', (byte) 'A'));
		v2.setQualityScore((byte) '"');
		v2.setPosition(23);

		List<ReadFeature> features = new ArrayList<ReadFeature>();
		features.add(v);
		features.add(v2);
		featuresCodec.write(bos, features);
		bos.write(true, 10) ;
		bos.write(false, 10) ;

		bos.flush();
		System.out.println(Utils.toBitString(baos.toByteArray()));

		ByteArrayInputStream bais = new ByteArrayInputStream(baos.toByteArray());
		BitInputStream bis = new DefaultBitInputStream(bais);
		Collection<ReadFeature> features2 = featuresCodec.read(bis);

		assertThat(features2, notNullValue());
		assertThat(features2.isEmpty(), is(false));
		assertThat(features2.size(), is(2));
		assertThat(features2.iterator().next(),
				instanceOf(SubstitutionVariation.class));

		Iterator<ReadFeature> iterator = features2.iterator();
		SubstitutionVariation v3 = (SubstitutionVariation) iterator.next();

		assertThat(v3, equalTo(v));

		SubstitutionVariation v4 = (SubstitutionVariation) iterator.next();
		
		for (int i=0; i<10; i++)
			 if (!bis.readBit()) 
				 throw new RuntimeException("No magick") ;
			for (int i=0; i<10; i++)
				if (bis.readBit()) 
					throw new RuntimeException("No magick") ;

		assertThat(v4, equalTo(v2));

	}

}

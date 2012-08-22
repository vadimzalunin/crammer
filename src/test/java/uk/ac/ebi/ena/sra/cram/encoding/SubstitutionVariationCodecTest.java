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
		HuffmanTree<Byte> qualityScoreTree = HuffmanCode
				.buildTree(new int[] { 100, 90, 80 }, new Byte[] { 33, 34, 35 });

		SubstitutionVariation v = new SubstitutionVariation();
		v.setBaseChange(new BaseChange((byte) 'A', (byte) 'C'));
		v.setPosition(-1);
		codec.write(bos, v);

		bos.flush();

		ByteArrayInputStream bais = new ByteArrayInputStream(baos.toByteArray());
		BitInputStream bis = new DefaultBitInputStream(bais);
		SubstitutionVariation v2 = codec.read(bis);

		assertThat(v2, equalTo(v));

	}
}

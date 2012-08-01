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

import java.io.IOException;

import uk.ac.ebi.ena.sra.cram.format.InsertionVariation;
import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;
import uk.ac.ebi.ena.sra.cram.io.NullBitOutputStream;

public class InsertionVariationCodec implements BitCodec<InsertionVariation> {
	public BitCodec<byte[]> insertBasesCodec;
	private long byteCounter = 0 ;
	private long totalLen = 0 ;

	@Override
	public InsertionVariation read(BitInputStream bis) throws IOException {
		// position is not read here because we need to keep track of previous
		// values read from the codec. See ReadFeatureCodec.
		long position = -1L;
		byte[] insertion = insertBasesCodec.read(bis);

		InsertionVariation v = new InsertionVariation();
		v.setPosition((int) position);
		v.setSequence(insertion);
		return v;
	}

	@Override
	public long write(BitOutputStream bos, InsertionVariation v)
			throws IOException {
		long len = 0L;

		len += insertBasesCodec.write(bos, v.getSequence());
//		System.out.println(new String(v.getSequence()) + " encoded in " +  len + " bits.");
		
		byteCounter += v.getSequence().length ;
		totalLen += len ;

		return len;
	}

	@Override
	public long numberOfBits(InsertionVariation v) {
		try {
			return write(NullBitOutputStream.INSTANCE, v);
		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
	
	@Override
	public String toString() {
		return String.format("Insertion codec: %d bases total, %d bits, %.2f bits per base.", byteCounter, totalLen, (float)totalLen/byteCounter) ;
	}

}

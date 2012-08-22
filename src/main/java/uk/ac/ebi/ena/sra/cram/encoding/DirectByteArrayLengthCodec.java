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

import uk.ac.ebi.ena.sra.cram.Utils;
import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;

public class DirectByteArrayLengthCodec implements BitCodec<byte[]> {

	@Override
	public byte[] read(BitInputStream bis) throws IOException {
		int len = bis.readByte();
		len <<= 8;
		len |= bis.readByte();
		len <<= 8;
		len |= bis.readByte();
		len <<= 8;
		len |= bis.readByte();

		byte[] bytes = new byte[len];
		bis.readAlignedBytes(bytes);

		return bytes;
	}

	@Override
	public long write(BitOutputStream bos, byte[] object) throws IOException {
		bos.write(Utils.toBytes(object.length));
		bos.write(object);
		return 8 * (4 + object.length);
	}

	@Override
	public long numberOfBits(byte[] object) {
		return 8 * (4 + object.length);
	}

}

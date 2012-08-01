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

import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;

public class BetaCodec implements BitCodec<Long> {
	private long offset = 0L;
	private int readNofBits;

	@Override
	public final Long read(BitInputStream bis) throws IOException {
		return bis.readLongBits(readNofBits) - offset;
	}

	@Override
	public final long write(BitOutputStream bos, Long value) throws IOException {
		if (value + offset < 0)
			throw new IllegalArgumentException("Value is less then offset: "
					+ value);

		int nofBits = (int) numberOfBits(value);
		long newValue = value + offset;
		bos.write(newValue, nofBits);
		return nofBits;
	}

	@Override
	public final long numberOfBits(Long value) {
		if (value > (1L << readNofBits))
			throw new IllegalArgumentException(
					"Value written is bigger then allowed: value=" + value
							+ ", max nof bits=" + readNofBits);
		
		return readNofBits;
		// long newValue = value + offset;
		// return (int) Math.floor(Math.log(newValue) / Math.log(2));
	}

	public long getOffset() {
		return offset;
	}

	public void setOffset(long offset) {
		this.offset = offset;
	}

	public int getReadNofBits() {
		return readNofBits;
	}

	public void setReadNofBits(int readNofBits) {
		this.readNofBits = readNofBits;
	}

}

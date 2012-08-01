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

public class GolombCodec implements BitCodec<Long> {
	private long m;
	private boolean quotientBit = true;
	private Long offset = 0L;
	private boolean deltaCodec = false;
	private long previousValue = 0L;

	public GolombCodec(int m) {
		this(m, true, 0L);
	}

	public GolombCodec(long m, boolean quotientBit, Long offset) {
		if (m < 2)
			throw new IllegalArgumentException("M parameter must be at least 2.");
		this.m = m;
		this.quotientBit = quotientBit;
		this.offset = offset;
	}

	@Override
	public final Long read(final BitInputStream bis) throws IOException {
		long quotient = 0L;
		while (bis.readBit() == quotientBit)
			quotient++;

		long numbits = quotient + 1;
		long ceiling = (long) (Math.log(m) / Math.log(2)+1);
		long reminder = bis.readBits((int) (ceiling - 1));
		numbits += ceiling - 1;
		if (reminder >= Math.pow(2, ceiling) - m) {
			reminder <<= 1;
			reminder |= bis.readBits(1);
			reminder -= Math.pow(2, ceiling) - m;
		}

		return (quotient * m + reminder) - offset;
	}

	@Override
	public final long write(final BitOutputStream bos, final Long value)
			throws IOException {
		long newValue = value + offset;
		long quotient = (long) (newValue / m);
		long reminder = newValue % m;
		long ceiling = (long) (Math.log(m) / Math.log(2)+1);

		long len = quotient + 1;
		bos.write(quotientBit, quotient);
		bos.write(!quotientBit);

		if (reminder < Math.pow(2, ceiling) - m) {
			bos.write(reminder, (int) ceiling - 1);
			len += ceiling - 1;
		} else {
			bos.write((int) (reminder + Math.pow(2, ceiling) - m),
					(int) ceiling);
			len += ceiling;
		}
		return len;
	}

	@Override
	public final long numberOfBits(Long value) {
		long newValue = value + offset;
		long quotient = (long) (newValue / m);
		long reminder = newValue % m;
		long ceiling = (long) (Math.log(m) / Math.log(2)+1);
		long l = quotient + 1;

		if (reminder < Math.pow(2, ceiling) - m)
			l += ceiling - 1;
		else
			l += ceiling;

		return l;
	}

	public long getM() {
		return m;
	}

	public boolean isQuotientBit() {
		return quotientBit;
	}

	public Long getOffset() {
		return offset;
	}

	public void setM(long m) {
		this.m = m;
	}

	public void setQuotientBit(boolean quotientBit) {
		this.quotientBit = quotientBit;
	}

	public void setOffset(Long offset) {
		this.offset = offset;
	}

	// @Override
	// public void reset() {
	// // TODO Auto-generated method stub
	//
	// }

}

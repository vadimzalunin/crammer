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
package net.sf.cram.encoding;

import java.io.IOException;

import net.sf.cram.io.BitInputStream;
import net.sf.cram.io.BitOutputStream;


public class GolombLongCodec implements BitCodec<Long> {
	private int m;
	private boolean quotientBit = true;
	private long offset = 0L;

	public GolombLongCodec(int m) {
		this(0, m, true);
	}

	public GolombLongCodec(long offset, int m) {
		this(offset, m, true);
	}

	public GolombLongCodec(long offset, int m, boolean quotientBit) {
		if (m < 2)
			throw new IllegalArgumentException(
					"M parameter must be at least 2.");
		this.m = m;
		this.quotientBit = quotientBit;
		this.offset = offset;
	}

	@Override
	public final Long read(final BitInputStream bis) throws IOException {
		long quotient = 0L;
		while (bis.readBit() == quotientBit)
			quotient++;

		long ceiling = (long) (Math.log(m) / Math.log(2) + 1);
		long reminder = bis.readBits((int) (ceiling - 1));
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
		long ceiling = (long) (Math.log(m) / Math.log(2) + 1);

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
		long ceiling = (long) (Math.log(m) / Math.log(2) + 1);
		long l = quotient + 1;

		if (reminder < Math.pow(2, ceiling) - m)
			l += ceiling - 1;
		else
			l += ceiling;

		return l;
	}

	public int getM() {
		return m;
	}

	public boolean isQuotientBit() {
		return quotientBit;
	}

	public Long getOffset() {
		return offset;
	}

	public void setM(int m) {
		this.m = m;
	}

	public void setQuotientBit(boolean quotientBit) {
		this.quotientBit = quotientBit;
	}

	public void setOffset(Long offset) {
		this.offset = offset;
	}

	@Override
	public Long read(BitInputStream bis, int len) throws IOException {
		throw new RuntimeException("Multi-value read method not defined.");
	}
}

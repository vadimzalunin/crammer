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

public class SubexpCodec implements BitCodec<Long> {
	private long offset = 0L;
	private long k = 2;
	private boolean unaryBit = true;

	public SubexpCodec(long offset, long k, boolean unaryBit) {
		super();
		this.offset = offset;
		this.k = k;
		this.unaryBit = unaryBit;
	}

	public SubexpCodec(long k) {
		this.k = k;
	}

	@Override
	public final Long read(BitInputStream bis) throws IOException {
		long u = 0L;
		while (bis.readBit() == unaryBit)
			u++;

		long nubits = u + 1;
		long b = 0;
		long n = 0;
		if (u == 0) {
			b = k;
			n = bis.readBits((int) b);
		} else {
			b = u + k - 1;
			n = (1 << b) | bis.readBits((int) b);
		}

		nubits += b;
		return n - offset;
	}

	@Override
	public final long write(BitOutputStream bos, Long value) throws IOException {
		if (value + offset < 0)
			throw new IllegalArgumentException("Value is less then offset: "
					+ value);

		long newValue = value + offset;
		long b = 0;
		long u = 0;
		if (newValue < (1L << k)) {
			b = k;
			u = 0;
		} else {
			b = (long) (Math.log(newValue) / Math.log(2));
			u = b - k + 1;
		}

		bos.write(unaryBit, u);
		bos.write(!unaryBit);

		bos.write(newValue, (int) (b));
		return u + 1 + b;
	}

	@Override
	public final long numberOfBits(Long value) {
		long newValue = value + offset;
		long b = 0;
		long u = 0;
		if (newValue < (1L << k)) {
			b = k;
			u = 0;
		} else {
			b = (long) Math.floor(Math.log(newValue) / Math.log(2));
			u = b - k + 1;
		}
		return u + 1 + b;
	}

	public long getOffset() {
		return offset;
	}

	public void setOffset(long offset) {
		this.offset = offset;
	}

	public long getK() {
		return k;
	}

	public void setK(long k) {
		this.k = k;
	}

	public boolean isUnaryBit() {
		return unaryBit;
	}

	public void setUnaryBit(boolean unaryBit) {
		this.unaryBit = unaryBit;
	}

}

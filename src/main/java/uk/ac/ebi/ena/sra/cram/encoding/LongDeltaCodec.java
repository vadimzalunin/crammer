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

public class LongDeltaCodec implements BitCodec<Long> {
	private BitCodec<Long> delegate;
	private long previousValue = 0L;

	public LongDeltaCodec(BitCodec<Long> delegate, long offset) {
		this.delegate = delegate;
		this.previousValue = offset;
	}

	@Override
	public Long read(BitInputStream bis) throws IOException {
		long value = previousValue + delegate.read(bis);
		previousValue = value;
		return value;
	}

	@Override
	public long write(BitOutputStream bos, Long object) throws IOException {
		long value = delegate.write(bos, object - previousValue);
		previousValue = value;
		return value;
	}

	@Override
	public long numberOfBits(Long object) {
		return delegate.numberOfBits(object - previousValue);
	}

	public void reset() {
		previousValue = 0L;
	}

}

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

public class BaseSequenceCodec implements BitCodec<byte[]> {
	public enum BaseCodecType {
		FLAT, RAISED, STEEP;
	}

	private BitCodec<byte[]> delegate;

	public BaseSequenceCodec(BaseCodecType type, byte[] order) {
		switch (type) {
		case FLAT:
			delegate = new BaseSequenceFlatCodec(order);
			break;
		case RAISED:
			delegate = new BaseSequenceRaisedCodec(order);
			break;
		case STEEP:
			delegate = new BaseSequenceSteepCodec(order);
			break;

		default:
			throw new IllegalArgumentException("Unknown base codec type: "
					+ type);
		}
	}

	@Override
	public byte[] read(BitInputStream bis) throws IOException {
		return delegate.read(bis);
	}

	@Override
	public long write(BitOutputStream bis, byte[] object) throws IOException {
		return delegate.write(bis, object);
	}

	@Override
	public long numberOfBits(byte[] object) {
		return delegate.numberOfBits(object);
	}

}

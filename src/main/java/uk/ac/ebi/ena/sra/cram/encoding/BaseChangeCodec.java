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

public class BaseChangeCodec implements BitCodec<BaseChange> {

	@Override
	public BaseChange read(BitInputStream bis) throws IOException {
		return new BaseChange(bis.readBits(2));
	}

	@Override
	public long write(BitOutputStream bis, BaseChange baseChange) throws IOException {
		bis.write(baseChange.getChange(), 2);
		return 2;
	}

	@Override
	public long numberOfBits(BaseChange baseChange) {
		return 2;
	}

}

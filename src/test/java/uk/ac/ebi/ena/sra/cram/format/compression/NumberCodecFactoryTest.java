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
package uk.ac.ebi.ena.sra.cram.format.compression;

import static org.hamcrest.CoreMatchers.*;
import static org.junit.Assert.assertThat;

import org.junit.Test;

public class NumberCodecFactoryTest {

	@Test
	public void test_all_encodings_have_codec() throws CramCompressionException {
		for (EncodingAlgorithm encoding : EncodingAlgorithm.values()) {
			NumberCodecStub stub = NumberCodecFactory.createStub(encoding);
			assertThat(stub, is(notNullValue()));
			assertThat(stub.getEncoding(), equalTo(encoding));
		}
	}

	@Test
	public void test_string_representation_round_trip()
			throws CramCompressionException {
		for (EncodingAlgorithm encoding : EncodingAlgorithm.values()) {
			NumberCodecStub stub = NumberCodecFactory.createStub(encoding);
			String stringRepresentation = stub.getStringRepresentation();
			stub.initFromString(stringRepresentation);
			assertThat(
					"String representaion round trip failed for " + encoding,
					stub.getStringRepresentation(),
					equalTo(stringRepresentation));
		}
	}
}

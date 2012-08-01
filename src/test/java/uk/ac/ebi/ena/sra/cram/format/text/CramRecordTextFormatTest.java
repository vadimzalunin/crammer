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
package uk.ac.ebi.ena.sra.cram.format.text;

import static org.hamcrest.CoreMatchers.equalTo;
import static org.junit.Assert.assertThat;

import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import uk.ac.ebi.ena.sra.cram.format.CramRecord;

@RunWith(value = Parameterized.class)
public class CramRecordTextFormatTest {
	private String recordSpec;

	public CramRecordTextFormatTest(String recordSpec) {
		this.recordSpec = recordSpec;
	}

	@Parameters
	public static Collection<Object[]> data() {
		Object[][] data = new Object[][] {
				{ "1	223022	54	*	POS	M26D1M21y M5g 	*	*" },
				{ "*	*	54	*	POS	*	AAAAAAAAAAAAAAAA	!!!!!!!!!!!!!!!!" },

		};
		return Arrays.asList(data);
	}

	@Test
	public void test1() throws IOException {
		CramRecordFormat format = new CramRecordFormat();

		CramRecord record = format.fromString(recordSpec);

		String derivedSpec = format.writeRecord(record);
		
		System.out.println(record.toString());
		System.out.println(derivedSpec);

		assertThat(derivedSpec, equalTo(recordSpec));
		
		CramRecord derivedRecord = format.fromString(derivedSpec) ;
		assertThat(derivedRecord, equalTo(record)) ;
	}

}

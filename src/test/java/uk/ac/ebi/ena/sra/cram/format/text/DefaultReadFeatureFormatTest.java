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

import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import uk.ac.ebi.ena.sra.cram.format.ReadFeature;

@RunWith(value = Parameterized.class)
public class DefaultReadFeatureFormatTest {
	private String spec;
	private int readLen;
	private DefaultReadFeaturesFormat format = new DefaultReadFeaturesFormat();

	public DefaultReadFeatureFormatTest(String spec, int readLen) {
		super();
		this.spec = spec;
		this.readLen = readLen;
	}

	@Parameters
	public static Collection<Object[]> data() {
		Object[][] data = new Object[][] { { "", 0 }, { "M5", 5 },
				{ "M2wxyzT$M2IAC.M5D2A!C!M13A!C!", 31 },

		};
		return Arrays.asList(data);
	}

	@Test
	public void test1() {
		List<ReadFeature> features = format.asReadFeatureList(spec);

		StringBuilder sb = new StringBuilder();
		format.addFeaturesToStringBuilder(features, readLen, sb);

		assertThat(sb.toString(), equalTo(spec));
	}
}

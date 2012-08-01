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
package uk.ac.ebi.ena.sra.cram.impl;

import static org.hamcrest.CoreMatchers.equalTo;
import static org.junit.Assert.assertThat;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import net.sf.samtools.Cigar;

import org.junit.Ignore;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import uk.ac.ebi.ena.sra.cram.format.ReadFeature;
import uk.ac.ebi.ena.sra.cram.format.text.DefaultReadFeaturesFormat;

@RunWith(value = Parameterized.class)
public class ReadFeatures2CigarTest {
	private String readFeaturesString;
	private String cigarString;
	private int readLength;

	public ReadFeatures2CigarTest(String readFeaturesString,
			String cigarString, int readLength) {
		super();
		this.readFeaturesString = readFeaturesString;
		this.cigarString = cigarString;
		this.readLength = readLength;
	}

	@Parameters
	public static Collection<Object[]> data() {
		Object[][] data = new Object[][] { 
				{ "a c g t M70a ", "4S70M1I1M", 76 },
				{ "M39t c c a c t c c a g c c a a c!", "39M15S", 54 },
				{ "M36", "36M", 36 },
				{ "M2y M19z M25y M5", "54M", 54 },
				{ "M27x M6x M3ITCCACTCCAGCCAA.M1c ", "38M14I1M1S", 54 },
				{ "z M2t M1x M48", "3M1I50M", 54 },
				{ "M26D1M21y M4t M1", "26M1D26M1I1M", 54 },
				{ "M85a M2g ", "85M1I2M1I1M", 90 },
				{ "a c g t a c g t a c g t a c g t M23a c g t a c g t a c g t a c g ", "16S23M15S", 54 },
				{ "M18D2D703M32", "18M705D32M", 50 },
				{ "M18D703D2M32", "18M705D32M", 50 },
				
//				9M2D398N1D41M vs 9M401D
//				4S3I43M vs 7S43M
//				5S1I6M3177N38M vs 6S6M3177D38M
//				22M2899N22M1I5S vs 22M2899D22M6S
//				22M2899N22M1I5S vs 22M2899D22M6S
//				31M1D2420N1D19M vs 31M2422D
//				22M1D118N1D28M vs 22M120D

		};
		return Arrays.asList(data);
	}

	@Test
	public void test() {
		DefaultReadFeaturesFormat drfFormat = new DefaultReadFeaturesFormat();
		List<ReadFeature> features = drfFormat
				.asReadFeatureList(readFeaturesString);

		ReadFeatures2Cigar rf2Cigar = new ReadFeatures2Cigar();
		Cigar cigar = rf2Cigar.getCigar2(features, readLength);

		assertThat(cigar.toString(), equalTo(cigarString));
	}
}

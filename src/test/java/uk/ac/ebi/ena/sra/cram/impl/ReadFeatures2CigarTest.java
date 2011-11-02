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
//				{ "a c g t M70a ", "4I70M1I1M", 76 },
//				{ "M39t c c a c t c c a g c c a a c!", "39M15I", 54 },
//				{ "M36", "36M", 36 },
//				{ "M2y M19z M25y M5", "54M", 54 },
//				{ "M27x M6x M3ITCCACTCCAGCCAA.M1c ", "38M14I1M1I", 54 },
//				{ "z M2t M1x M48", "3M1I50M", 54 },
//				{ "M26D1M21y M4t M1", "26M1D26M1I1M", 54 },
//				{ "M85a M2g ", "85M1I2M1I1M", 90 },
				{ "a c g t a c g t a c g t a c g t M23a c g t a c g t a c g t a c g ", "16S23M15S", 54 },

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

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

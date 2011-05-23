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

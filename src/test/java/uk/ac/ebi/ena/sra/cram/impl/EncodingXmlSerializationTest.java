package uk.ac.ebi.ena.sra.cram.impl;

import static org.hamcrest.CoreMatchers.*;
import static org.junit.Assert.assertThat;

import java.beans.XMLDecoder;
import java.beans.XMLEncoder;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;

import org.junit.Test;

import uk.ac.ebi.ena.sra.cram.format.Encoding;
import uk.ac.ebi.ena.sra.cram.format.compression.EncodingAlgorithm;

public class EncodingXmlSerializationTest {

	@Test
	public void test1() {
		Encoding encoding = new Encoding(EncodingAlgorithm.GOLOMB, "1,0,1");
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		XMLEncoder xmlEncoder = new XMLEncoder(baos);
		xmlEncoder.writeObject(encoding);
		xmlEncoder.close();

		System.out.println(new String (baos.toByteArray()));
		ByteArrayInputStream bais = new ByteArrayInputStream(baos.toByteArray());
		XMLDecoder xmlDecoder = new XMLDecoder(bais);
		Object object = xmlDecoder.readObject();

		assertThat(object, notNullValue());
		assertThat(object, instanceOf(Encoding.class));
		
		Encoding encoding2 = (Encoding) object ;
		
		assertThat(encoding2, equalTo(encoding)) ;

	}

}

package uk.ac.ebi.ena.sra.cram.impl;

import static org.hamcrest.CoreMatchers.equalTo;
import static org.junit.Assert.assertThat;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import uk.ac.ebi.ena.sra.cram.encoding.BitCodec;
import uk.ac.ebi.ena.sra.cram.format.CramCompression;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.CramRecordBlock;
import uk.ac.ebi.ena.sra.cram.format.compression.CramCompressionException;
import uk.ac.ebi.ena.sra.cram.format.text.CramRecordFormat;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitInputStream;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitOutputStream;
import uk.ac.ebi.ena.sra.cram.stats.CramStats;

@RunWith(value = Parameterized.class)
public class CramRecordCodecRoundTripTests {
	private String recordSpec;

	public CramRecordCodecRoundTripTests(String recordSpec) {
		this.recordSpec = recordSpec;
	}

	@Parameters
	public static Collection<Object[]> data() {
		Object[][] data = new Object[][] {
				{ "123	*	*	POS	*	*	*" },
				{ "124	36	*	POS	M2w!x!y!z!T$M2IAC.M5D2A!C!M13A!C!M2	*	*" },
				{ "*	34	*	NEG	*	GTGCGGATGCTCTCCTCCAGTTTGGGCTCGTGGTGTGTGTCCAGCAGGGACTGG	BBBBBBB=BBBBBBBBBBBBBBBB?BBBBB?BB?BBBBB?BBBBB?BBB??BBB" },

		};
		return Arrays.asList(data);
	}

	@Test
	public void test1() throws IOException, CramCompressionException {
		CramRecordFormat format = new CramRecordFormat();
		CramRecord record = format.fromString(recordSpec);

		CramStats stats = new CramStats();
		Logger.getLogger(CramStats.class).setLevel(Level.ERROR);
		stats.addRecord(record);

		RecordCodecFactory factory = new RecordCodecFactory();
		CramRecordBlock block = new CramRecordBlock();
		block.setUnmappedReadQualityScoresIncluded(true) ;
		CramCompression compression = new CramCompression();
		block.setCompression(compression);

		stats.adjustBlock(block);
		BitCodec<CramRecord> codec = factory.createRecordCodec(block, null);

		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		DefaultBitOutputStream bos = new DefaultBitOutputStream(baos);

		codec.write(bos, record);
		bos.close();

		ByteArrayInputStream bais = new ByteArrayInputStream(baos.toByteArray());
		DefaultBitInputStream bis = new DefaultBitInputStream(bais);
		CramRecord derivedRecord = codec.read(bis);

		assertThat(derivedRecord, equalTo(record));

	}
}
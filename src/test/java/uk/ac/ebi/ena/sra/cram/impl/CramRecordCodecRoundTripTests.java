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
import uk.ac.ebi.ena.sra.cram.format.CramHeader;
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
		Object[][] data = new Object[][] { { "1	223022	54	*	POS	M26D1M21y M5g 	*	*" },
				{ "*	*	16	*	POS	*	AAAAAAAAAAAAAAAA	!!!!!!!!!!!!!!!!" }, };
		return Arrays.asList(data);
	}

	@Test
	public void test1() throws IOException, CramCompressionException {
		CramRecordFormat format = new CramRecordFormat();
		CramRecord record = format.fromString(recordSpec);

		CramStats stats = new CramStats(new CramHeader(), null);
		Logger.getLogger(CramStats.class).setLevel(Level.ERROR);
		stats.addRecord(record);

		RecordCodecFactory factory = new RecordCodecFactory();
		CramRecordBlock block = new CramRecordBlock();
		block.setUnmappedReadQualityScoresIncluded(true);
		CramCompression compression = CramCompression.createDefaultCramCompression();
		block.setCompression(compression);

		stats.adjustBlock(block);
		BitCodec<CramRecord> codec = factory.createRecordCodec(null, block, null);

		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		DefaultBitOutputStream bos = new DefaultBitOutputStream(baos);

		codec.write(bos, record);
		bos.close();

		ByteArrayInputStream bais = new ByteArrayInputStream(baos.toByteArray());
		DefaultBitInputStream bis = new DefaultBitInputStream(bais);
		codec = factory.createRecordCodec(null, block, null);
		CramRecord derivedRecord = codec.read(bis);
		derivedRecord.setFlags(derivedRecord.getFlags());

		assertThat(derivedRecord, equalTo(record));

	}
}

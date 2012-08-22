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
package uk.ac.ebi.ena.sra.cram.encoding;

import static org.hamcrest.CoreMatchers.*;
import static org.junit.Assert.assertThat;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;

import org.junit.Test;
import org.junit.runner.RunWith;
import org.junit.runners.Parameterized;
import org.junit.runners.Parameterized.Parameters;

import uk.ac.ebi.ena.sra.cram.Utils;
import uk.ac.ebi.ena.sra.cram.encoding.BaseSequenceCodec.BaseCodecType;
import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitInputStream;
import uk.ac.ebi.ena.sra.cram.io.DefaultBitOutputStream;

@RunWith(value = Parameterized.class)
public class BaseSequenceCodecTest {
	private boolean dumpState = false;

	private BitCodec<byte[]> codec;
	private String sequence;
	private int expectedLength;

	public BaseSequenceCodecTest(BaseCodecType codecType, String order, String sequence, int expectedLength) {
		this.codec = new BaseSequenceCodec(codecType, order.getBytes());
		this.sequence = sequence;
		this.expectedLength = expectedLength;
	}

	@Parameters
	public static Collection<Object[]> data() {
		Object[][] data = new Object[][] { { BaseCodecType.FLAT, "ACGTNS", "", 3 },
				{ BaseCodecType.FLAT, "ACGTNS", "A", 6 }, { BaseCodecType.FLAT, "ACGTNS", "AA", 9 },
				{ BaseCodecType.FLAT, "ACGTNS", "ACGTN", 5 * 3 + 3 },
				{ BaseCodecType.FLAT, "ACGTNS", "ACGACTNNCNGACTACNC", "ACGACTNNCNGACTACNC".length() * 3 + 3 },
				{ BaseCodecType.RAISED, "ACGTNS", "", 4 }, { BaseCodecType.RAISED, "ACGTNS", "A", 6 },
				{ BaseCodecType.RAISED, "ACGTNS", "C", 6 }, { BaseCodecType.RAISED, "ACGTNS", "G", 6 },
				{ BaseCodecType.RAISED, "ACGTNS", "T", 7 }, { BaseCodecType.RAISED, "ACGTNS", "N", 8 },
				{ BaseCodecType.RAISED, "ACGTNS", "ACGACTNNCNGACTACNC", 50 }, { BaseCodecType.STEEP, "ACGTNS", "", 5 },
				{ BaseCodecType.STEEP, "ACGTNS", "A", 6 }, { BaseCodecType.STEEP, "ACGTNS", "C", 7 },
				{ BaseCodecType.STEEP, "ACGTNS", "G", 8 }, { BaseCodecType.STEEP, "ACGTNS", "T", 9 },
				{ BaseCodecType.STEEP, "ACGTNS", "N", 10 },
				{ BaseCodecType.STEEP, "ACGTNS", "ACGACTNNCNGACTACNC", 55 }, };
		return Arrays.asList(data);
	}

	@Test
	public void test() throws IOException {
		ByteArrayOutputStream baos = new ByteArrayOutputStream(1024);
		BitOutputStream bos = new DefaultBitOutputStream(baos);

		byte[] bytes = sequence.getBytes();

		long len = codec.write(bos, bytes);
		if (expectedLength > -1)
			assertThat((int) len, is(expectedLength));

		bos.flush();

		byte[] buf = baos.toByteArray();

		if (dumpState)
			System.out.printf("%s\t%s\t%d\t%s\n", sequence, Utils.toBitString(buf).substring(0, (int) len), len,
					Utils.toBitString(buf));

		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		BitInputStream bis = new DefaultBitInputStream(bais);

		byte[] restoredSequence = codec.read(bis);

		assertThat("Expecting " + sequence + " but got " + new String(restoredSequence), restoredSequence,
				equalTo(bytes));
	}
}

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
package uk.ac.ebi.ena.sra.cram.io;

import static org.hamcrest.CoreMatchers.is;
import static org.junit.Assert.assertThat;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.IOException;

import org.junit.Test;

public class BitStreamIntegrationTest {

	@Test
	public void test_3_bits_of_byte() throws IOException {
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		DefaultBitOutputStream bos = new DefaultBitOutputStream(baos);

		byte value = 1 << 1;
		bos.write(value, 3);
		bos.flush();

		ByteArrayInputStream bais = new ByteArrayInputStream(baos.toByteArray());
		DefaultBitInputStream bis = new DefaultBitInputStream(bais);
		int readBits = bis.readBits(3);

		assertThat((byte) readBits, is(value));
	}

	@Test
	public void test_all_possible_bytes() throws IOException {
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		DefaultBitOutputStream bos = new DefaultBitOutputStream(baos);

		for (byte value = Byte.MIN_VALUE; value < Byte.MAX_VALUE; value++) {
			baos.reset();
			bos.write(value, 8);
			bos.flush();

			ByteArrayInputStream bais = new ByteArrayInputStream(
					baos.toByteArray());
			DefaultBitInputStream bis = new DefaultBitInputStream(bais);
			int readBits = bis.readBits(8);

			assertThat((byte) readBits, is(value));
		}
	}

	@Test
	public void test_bits_1_to_8() throws IOException {
		for (int bit = 1; bit < 8; bit++) {
			ByteArrayOutputStream baos = new ByteArrayOutputStream();
			DefaultBitOutputStream bos = new DefaultBitOutputStream(baos);
			byte value = 1;
			bos.write(value, bit);
			bos.flush();

			ByteArrayInputStream bais = new ByteArrayInputStream(
					baos.toByteArray());
			DefaultBitInputStream bis = new DefaultBitInputStream(bais);
			int readBits = bis.readBits(bit);

			assertThat((byte) readBits, is(value));
		}
	}

	@Test
	public void testWrite_read_int_by_bit() throws IOException {
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		DefaultBitOutputStream bos = new DefaultBitOutputStream(baos);
		int value = 1;
		bos.write(value, 32);
		bos.flush();
		byte[] buf = baos.toByteArray();

		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		DefaultBitInputStream bis = new DefaultBitInputStream(bais);
		for (int i = 0; i < 31; i++) {
			int readBits = bis.readBits(1);
			assertThat(readBits, is(0));
		}

		int readBits = bis.readBits(1);
		assertThat(readBits, is(1));
	}

	@Test
	public void testWrite_read_32_bits() throws IOException {
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		DefaultBitOutputStream bos = new DefaultBitOutputStream(baos);
		int value = 1234567890;
		bos.write(value, 32);
		bos.flush();
		byte[] buf = baos.toByteArray();

		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		DefaultBitInputStream bis = new DefaultBitInputStream(bais);

		int readBits = bis.readBits(32);
		assertThat(readBits, is(value));
	}

	@Test
	public void testWrite_MaxInteger() throws IOException {
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		DefaultBitOutputStream bos = new DefaultBitOutputStream(baos);
		int value = Integer.MAX_VALUE;
		bos.write(value, 32);
		bos.flush();
		byte[] buf = baos.toByteArray();

		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		DefaultBitInputStream bis = new DefaultBitInputStream(bais);

		int readBits = bis.readBits(32);
		assertThat(readBits, is(value));
	}

	@Test
	public void testWrite_MinInteger() throws IOException {
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		DefaultBitOutputStream bos = new DefaultBitOutputStream(baos);
		int value = Integer.MIN_VALUE;
		bos.write(value, 32);
		bos.flush();
		byte[] buf = baos.toByteArray();

		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		DefaultBitInputStream bis = new DefaultBitInputStream(bais);

		int readBits = bis.readBits(32);
		assertThat(readBits, is(value));
	}

	@Test (timeout=1000)
	public void testWriteReadBenchmark() throws IOException {
		int maxValues = 1000000;
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		DefaultBitOutputStream bos = new DefaultBitOutputStream(baos);

		for (int value = 0; value < maxValues; value++)
			bos.write(value, 32);

		bos.flush();

		byte[] buf = baos.toByteArray();
		assertThat(buf.length, is(maxValues*4)) ;
		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		DefaultBitInputStream bis = new DefaultBitInputStream(bais);

		for (int value = 0; value < maxValues; value++) {

			int readBits = bis.readBits(32);
			assertThat(readBits, is(value));
		}
	}
	
	@Test// (timeout=1000)
	public void testWriteReadLongBenchmark() throws IOException {
		int maxValues = 1000000;
		ByteArrayOutputStream baos = new ByteArrayOutputStream();
		DefaultBitOutputStream bos = new DefaultBitOutputStream(baos);

		for (int value = 0; value < maxValues; value++)
			bos.write(value, 32);

		bos.flush();
		bos.close() ;

		byte[] buf = baos.toByteArray();
		assertThat(buf.length, is(maxValues*4)) ;
		ByteArrayInputStream bais = new ByteArrayInputStream(buf);
		DefaultBitInputStream bis = new DefaultBitInputStream(bais);

		for (long value = 0; value < maxValues; value++) {

			long readBits = bis.readLongBits(32);
			assertThat(readBits, is(value));
		}
		bis.close() ;
	}
	
	@Test
	public void test_Align() throws IOException {
		for (int bits=0; bits<64; bits++) {
			for (int markerLen = 0; markerLen<64; markerLen++) {
				long marker = ~(-1 << (markerLen / 2));
				
				ByteArrayOutputStream baos = new ByteArrayOutputStream();
				DefaultBitOutputStream bos = new DefaultBitOutputStream(baos);
//				DebuggingBitOuputStream bos = new DebuggingBitOuputStream(System.err, '\n', new DefaultBitOutputStream(baos)) ; 
				
				bos.write(0L, bits) ;
				bos.alignToByte() ;
				bos.write(marker, markerLen) ;
				bos.write(0L, bits) ;
				bos.write(marker, markerLen) ;
				
				bos.close() ;
//				bos.flush() ;
				
				ByteArrayInputStream bais = new ByteArrayInputStream(baos.toByteArray());
				DefaultBitInputStream bis = new DefaultBitInputStream(bais);
				
				assertThat(bis.readBits(bits), is(0)) ;
				bis.alignToByte() ;
				assertThat(bis.ensureMarker(marker, markerLen), is(true)) ;
				assertThat(bis.readBits(bits), is(0)) ;
				assertThat(bis.ensureMarker(marker, markerLen), is(true)) ;
				
				bis.close() ;
				
			}
		}
	}
}

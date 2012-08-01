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
package uk.ac.ebi.ena.sra.cram.mask;

import static org.hamcrest.CoreMatchers.equalTo;
import static org.junit.Assert.assertThat;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.StringReader;
import java.util.Arrays;

import org.junit.Test;

public class SingleLineMaskReaderTest {

	private static final String arrayToString(int[][] data) {
		StringBuilder sb = new StringBuilder();

		for (int[] array : data) {
			boolean first = true;
			for (int value : array) {
				if (!first) {
					sb.append(" ");
				} else
					first = false;
				sb.append(value);
			}
			sb.append("\n");
		}
		return sb.toString();
	}

	@Test
	public void testIntegerListMaskFactory() throws ReadMaskFormatException, IOException {
		int[][] data = new int[][] { {}, { 1 }, { 1, 3 },
				{ 1, 100, 101, 200, 600 }, {}, { 13 } };

		String content = arrayToString(data);
		ReadMaskReader reader = new SingleLineMaskReader(new BufferedReader(
				new StringReader(content)), new IntegerListMaskFactory());

		int[][] newData = new int[data.length][];
		int i = 0;
		PositionMask mask = null ;
		while ((mask = reader.readNextMask()) != null) {
			newData[i++] = mask.getMaskedPositions();
		}

		assertThat(newData, equalTo(data));
	}

	private static final String arrayToFasta(int[][] data) {
		StringBuilder sb = new StringBuilder();

		for (int[] array : data) {
			if (array.length > 0) {
				byte[] bytes = new byte[array[array.length-1]] ;
				Arrays.fill(bytes, FastaByteArrayMaskFactory.DEFAULT_NON_MASK_BYTE) ;
				for (int pos:array) {
					bytes[pos-1] = FastaByteArrayMaskFactory.DEFAULT_MASK_BYTE ;
				}
				sb.append(new String (bytes)) ;
			}
			sb.append("\n");
		}
		return sb.toString();
	}

	@Test
	public void testFastaByteArrayMaskFactory() throws ReadMaskFormatException,
			IOException {
		int[][] data = new int[][] { {}, { 1 }, { 1, 3 },
				{ 1, 100, 101, 200, 600 }, {}, { 13 } };

		String content = arrayToFasta(data);
		ReadMaskReader reader = new SingleLineMaskReader(new BufferedReader(
				new StringReader(content)), new FastaByteArrayMaskFactory());

		int[][] newData = new int[data.length][];
		int i = 0;
		PositionMask mask = null ;
		while ((mask = reader.readNextMask()) != null) {
			newData[i++] = mask.getMaskedPositions();
		}

		assertThat(newData, equalTo(data));
	}
}

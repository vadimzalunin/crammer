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
package uk.ac.ebi.ena.sra.cram.format;

import java.util.Arrays;

public class DiByteFrequencies {
	private int[][] frequencies;

	public DiByteFrequencies() {
		frequencies = new int[256][256];
	}

	public DiByteFrequencies(byte[][] values, int[] freqs) {
		frequencies = new int[256][256];
		for (int i = 0; i < values.length; i++)
			frequencies[0xFF & values[i][0]][0xFF & values[i][1]] = freqs[i];
	}

	public void add(byte value1, byte value2) {
		frequencies[0xFF & value1][0xFF & value2]++;
	}

	public void add(byte value1, byte value2, int howMany) {
		frequencies[0xFF & value1][0xFF & value2] += howMany;
	}

	public int getFrequency(byte value1, byte value2) {
		return frequencies[0xFF & value1][0xFF & value2];
	}

	public byte[][] getValues() {
		int size = 0;
		for (int i = 0; i < frequencies.length; i++) {
			for (int j = 0; j < frequencies[0].length; j++)
				if (frequencies[i][j] > 0)
					size++;
		}

		int valueIndex = 0;
		byte[][] collapsedValueArray = new byte[size][2];
		for (int i = 0; i < frequencies.length; i++) {
			for (int j = 0; j < frequencies[0].length; j++)
				if (frequencies[i][j] > 0) {
					collapsedValueArray[valueIndex][0] = (byte) i;
					collapsedValueArray[valueIndex++][1] = (byte) j;
				}
		}

		return collapsedValueArray;
	}

	public int[] getFrequencies() {
		int size = 0;
		for (int i = 0; i < frequencies.length; i++) {
			for (int j = 0; j < frequencies[0].length; j++)
				if (frequencies[i][j] > 0)
					size++;
		}

		int frequencyIndex = 0;
		int[] collapsedFrequencyArray = new int[size];
		for (int i = 0; i < frequencies.length; i++) {
			for (int j = 0; j < frequencies[0].length; j++)
				if (frequencies[i][j] > 0)
					collapsedFrequencyArray[frequencyIndex++] = frequencies[i][j];
		}

		return collapsedFrequencyArray;
	}

	@Override
	public String toString() {
		StringBuffer sb = new StringBuffer();
		sb.append("[");
		byte[][] values = getValues();
		for (int i = 0; i < values.length; i++) {
			sb.append((char) values[i][0]).append((char) values[i][1]);
			if (i < values.length - 1)
				sb.append(", ");
		}
		sb.append("]\n");
		sb.append(Arrays.toString(getFrequencies()));
		return sb.toString();
	}
}

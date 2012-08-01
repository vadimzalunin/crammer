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
package net.sf.picard.sam;

public class AlignmentSliceQuery {
	public String sequence;
	public int start;
	public int end;

	public AlignmentSliceQuery(String spec) {
		String[] chunks = spec.split(":");

		sequence = chunks[0];

		if (chunks.length > 1) {
			chunks = chunks[1].split("-");
			start = Integer.valueOf(chunks[0]);
			if (chunks.length == 2)
				end = Integer.valueOf(chunks[1]);
		}

	}
}

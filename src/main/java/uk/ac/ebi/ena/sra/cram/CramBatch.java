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
package uk.ac.ebi.ena.sra.cram;

import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

public class CramBatch {

	public static void main(String[] args) throws Exception {
		List<String> sources = new ArrayList<String>();
		if (args.length == 0)
			sources.add("cram.options");
		else
			sources.add(args[0]);

		for (String sourceFileName : sources) {
			File optionsFile = new File(sourceFileName);
			InputStream is;
			if (!optionsFile.exists())
				is = System.in;
			else
				is = new FileInputStream(optionsFile);
			Scanner scanner = new Scanner(is);
			String[] options = scanner.nextLine().split(" ");
			CramTools.main(options);
		}
	}
}

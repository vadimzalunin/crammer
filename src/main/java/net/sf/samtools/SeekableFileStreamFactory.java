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
package net.sf.samtools;

import java.io.File;
import java.io.FileNotFoundException;

import net.sf.samtools.util.SeekableFileStream;
import net.sf.samtools.util.SeekableStream;

public class SeekableFileStreamFactory implements SeekableInputStreamFactory {
	private File file;

	public SeekableFileStreamFactory(File file) {
		super();
		this.file = file;
	}

	@Override
	public SeekableStream createSeekableStream() throws FileNotFoundException {
		return new SeekableFileStream(file);
	}

}

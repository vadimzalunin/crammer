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
package uk.ac.ebi.ena.sra.compression;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.util.zip.GZIPOutputStream;

import org.apache.tools.bzip2.CBZip2OutputStream;

public class ArrayCompressionReport {
	private String name;
	// private long lzmaSize;
	private long gzipSize;
	private long bzip2Size;
	private long arraySize;

	public ArrayCompressionReport(String name) {
		this.name = name;
	}

	public void run(byte[] array) throws IOException {
		arraySize = array.length;
		//
		// ByteArrayOutputStream lzmaBAOS = new ByteArrayOutputStream();
		// SevenZip.Compression.LZMA.Encoder encoder = new
		// SevenZip.Compression.LZMA.Encoder();
		// encoder.Code(new ByteArrayInputStream(array), lzmaBAOS, -1, -1,
		// null);
		// lzmaBAOS.flush();
		// lzmaSize = lzmaBAOS.size();
		// lzmaBAOS.close();
		// lzmaBAOS = null;

		ByteArrayOutputStream gzipBAOS = new ByteArrayOutputStream();
		GZIPOutputStream gzipOS = new GZIPOutputStream(gzipBAOS);
		gzipOS.write(array);
		gzipOS.close();
		gzipOS = null;
		gzipBAOS.flush();
		gzipSize = gzipBAOS.size();
		gzipBAOS.close();
		gzipBAOS = null;

		ByteArrayOutputStream bzip2BAOS = new ByteArrayOutputStream();
		CBZip2OutputStream bz2OS = new CBZip2OutputStream(bzip2BAOS);
		bz2OS.write(array);
		bz2OS.close();
		bz2OS = null;
		bzip2BAOS.flush();
		bzip2Size = bzip2BAOS.size();
		bzip2BAOS.close();
		bzip2BAOS = null;
	}

	//
	// public long getLzmaSize() {
	// return lzmaSize;
	// }

	public long getGzipSize() {
		return gzipSize;
	}

	public long getBzip2Size() {
		return bzip2Size;
	}

	public long getArraySize() {
		return arraySize;
	}

	@Override
	public String toString() {
		StringBuffer sb = new StringBuffer(name == null ? "ArrayCompressionReport" : name);
		sb.append(": [");
		sb.append("array: ").append(arraySize).append("; ");
		sb.append("bzip2: ").append(bzip2Size).append("; ");
		sb.append("gzip: ").append(gzipSize).append("] ");
		// sb.append("lzma: ").append(lzmaSize).append("] ");
		return sb.toString();
	}
}

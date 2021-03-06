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
package net.sf.cram.io;

import java.io.IOException;

public interface BitOutputStream {

	public void write(int b, int nbits) throws IOException;

	public void write(long b, int nbits) throws IOException;

	public void write(byte b, int nbits) throws IOException;

	public void write(boolean bit) throws IOException;

	public void write(boolean bit, long repeat) throws IOException;

	public void flush() throws IOException;

	public void close() throws IOException;

	public int alignToByte() throws IOException;

	public void write(byte[] data) throws IOException;

	public void write(byte b) throws IOException;
	
}

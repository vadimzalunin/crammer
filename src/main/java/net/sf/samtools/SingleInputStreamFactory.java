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

import java.io.FileNotFoundException;
import java.io.InputStream;

public class SingleInputStreamFactory implements InputStreamFactory {
	private InputStream delegate;

	public SingleInputStreamFactory(InputStream delegate) {
		this.delegate = delegate;
	}

	@Override
	public InputStream createInputStream() throws FileNotFoundException,
			ExhaustedFactoryException {
		if (delegate == null)
			throw new ExhaustedFactoryException();
		InputStream is = delegate;
		delegate = null;
		return is;
	}

}

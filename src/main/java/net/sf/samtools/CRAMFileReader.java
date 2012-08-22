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

import java.io.IOException;
import java.io.InputStream;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.util.CloseableIterator;
import net.sf.samtools.util.SeekableStream;
import uk.ac.ebi.ena.sra.cram.CramException;
import uk.ac.ebi.ena.sra.cram.Utils;
import uk.ac.ebi.ena.sra.cram.format.CramFormatException;
import uk.ac.ebi.ena.sra.cram.format.CramHeader;
import uk.ac.ebi.ena.sra.cram.format.compression.CramCompressionException;
import uk.ac.ebi.ena.sra.cram.impl.CRAMPreemptiveIterator;
import uk.ac.ebi.ena.sra.cram.impl.CramIterator;
import uk.ac.ebi.ena.sra.cram.impl.CramPreemptiveRandomAccessIterator;
import uk.ac.ebi.ena.sra.cram.index.CramIndex;
import uk.ac.ebi.ena.sra.cram.index.RecordPointer;

public class CRAMFileReader extends SAMFileReader.ReaderImplementation {

	private SeekableInputStreamFactory sisFactory = null;
	private InputStreamFactory isFactory = null;

	private ReferenceSequenceFile referenceSequenceFile;

	private CramIndex cramIndex;
	private SAMFileHeader samFileHeader;

	public CRAMFileReader(SeekableInputStreamFactory sisFactory, InputStreamFactory isFactory,
			ReferenceSequenceFile referenceSequenceFile, CramIndex cramIndex) {
		this.sisFactory = sisFactory;
		this.isFactory = isFactory;
		this.referenceSequenceFile = referenceSequenceFile;
		this.cramIndex = cramIndex;
	}

	@Override
	SAMFileHeader getFileHeader() {
		if (samFileHeader == null)
			try {
				readSAMFileHeader();
			} catch (IOException e) {
				throw new RuntimeException(e);
			} catch (ExhaustedFactoryException e) {
				throw new RuntimeException(e);
			}
		return samFileHeader;
	}

	void readSAMFileHeader() throws IOException, ExhaustedFactoryException {
		CramIterator it = null;
		try {
			it = new CramIterator(isFactory.createInputStream(), referenceSequenceFile);
		} catch (CramFormatException e) {
			System.err.println("Cram format exception: " + e.getMessage());
			throw new RuntimeException(e);
		}
		CramHeader cramHeader = it.getCramHeader();
		samFileHeader = Utils.cramHeader2SamHeader(cramHeader);
	}

	CloseableIterator<SAMRecord> getIterator() {
		InputStream is = null;
		try {
			is = isFactory.createInputStream();
		} catch (IOException e) {
			throw new RuntimeException(e);
		} catch (ExhaustedFactoryException e) {
			throw new RuntimeException(e);
		}
		CRAMPreemptiveIterator iterator = null;
		try {
			iterator = new CRAMPreemptiveIterator(is, referenceSequenceFile, null);
		} catch (CramFormatException e) {
			throw new RuntimeException(e);
		} catch (CramCompressionException e) {
			throw new RuntimeException(e);
		} catch (IOException e) {
			throw new RuntimeException(e);
		} catch (CramException e) {
			throw new RuntimeException(e);
		}

		return new CachingSAMRecordIterator(iterator, iterator.getCramHeader());
	}

	@Override
	void close() {

	}

	@Override
	void setValidationStringency(ValidationStringency validationStringency) {
		// TODO Auto-generated method stub

	}

	@Override
	ValidationStringency getValidationStringency() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	void enableFileSource(SAMFileReader reader, boolean enabled) {
		// throw new
		// RuntimeException("Method not defined for CRAM file reader.");
	}

	@Override
	void enableIndexCaching(boolean enabled) {
		// throw new
		// RuntimeException("Method not defined for CRAM file reader.");
	}

	@Override
	void enableIndexMemoryMapping(boolean enabled) {
		// throw new
		// RuntimeException("Method not defined for CRAM file reader.");
	}

	@Override
	void enableCrcChecking(boolean enabled) {
		// throw new
		// RuntimeException("Method not defined for CRAM file reader.");
	}

	@Override
	boolean hasIndex() {
		return cramIndex != null;
	}

	@Override
	BAMIndex getIndex() {
		throw new RuntimeException("Not a BAM file.");
	}

	@Override
	CloseableIterator<SAMRecord> getIterator(SAMFileSpan fileSpan) {
		throw new RuntimeException("Method not defined for CRAM file reader.");
	}

	@Override
	SAMFileSpan getFilePointerSpanningReads() {
		throw new RuntimeException("Method not defined for CRAM file reader.");
	}

	@Override
	CloseableIterator<SAMRecord> query(String sequence, int start, int end, boolean contained) {
		if (cramIndex == null || sisFactory == null) {
			InputStream is = null;
			try {
				is = isFactory.createInputStream();
			} catch (IOException e) {
				throw new RuntimeException(e);
			} catch (ExhaustedFactoryException e) {
				throw new RuntimeException(e);
			}
			CRAMPreemptiveIterator iterator = null;
			try {
				iterator = new CRAMPreemptiveIterator(is, referenceSequenceFile, null);
			} catch (CramFormatException e) {
				throw new RuntimeException(e);
			} catch (CramCompressionException e) {
				throw new RuntimeException(e);
			} catch (IOException e) {
				throw new RuntimeException(e);
			} catch (CramException e) {
				throw new RuntimeException(e);
			}

			return new PreemptiveSAMRecordIterator(new CachingSAMRecordIterator(iterator, iterator.getCramHeader()),
					start, end, contained);
		}
		if (sisFactory != null) {
			SeekableStream is = null;
			try {
				is = sisFactory.createSeekableStream();
			} catch (IOException e) {
				throw new RuntimeException(e);
			} catch (ExhaustedFactoryException e) {
				throw new RuntimeException(e);
			}
			CramPreemptiveRandomAccessIterator iterator = null;
			try {
				RecordPointer pointer = cramIndex.findRecordPointerAt(sequence, start);
				iterator = new CramPreemptiveRandomAccessIterator(is, referenceSequenceFile, pointer);
			} catch (CramFormatException e) {
				throw new RuntimeException(e);
			} catch (CramCompressionException e) {
				throw new RuntimeException(e);
			} catch (IOException e) {
				throw new RuntimeException(e);
			} catch (CramException e) {
				throw new RuntimeException(e);
			}

			return new PreemptiveSAMRecordIterator(new CachingSAMRecordIterator(iterator, iterator.getCramHeader()),
					start, end, contained);
		}

		// need skipable input stream:
		throw new RuntimeException("CRAM index found but the stream is not random access.");
	}

	@Override
	CloseableIterator<SAMRecord> queryAlignmentStart(String sequence, int start) {
		return query(sequence, start, -1, true);
	}

	@Override
	public CloseableIterator<SAMRecord> queryUnmapped() {
		InputStream is = null;
		try {
			is = isFactory.createInputStream();
		} catch (IOException e) {
			throw new RuntimeException(e);
		} catch (ExhaustedFactoryException e) {
			throw new RuntimeException(e);
		}
		CRAMPreemptiveIterator iterator = null;
		try {
			iterator = new CRAMPreemptiveIterator(is, referenceSequenceFile, null);
		} catch (CramFormatException e) {
			throw new RuntimeException(e);
		} catch (CramCompressionException e) {
			throw new RuntimeException(e);
		} catch (IOException e) {
			throw new RuntimeException(e);
		} catch (CramException e) {
			throw new RuntimeException(e);
		}

		return new PreemptiveUnmappedSAMRecordIterator(new CachingSAMRecordIterator(iterator, iterator.getCramHeader()));
	}

}

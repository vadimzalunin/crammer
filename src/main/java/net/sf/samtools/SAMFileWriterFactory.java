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
/*
 * The MIT License
 *
 * Copyright (c) 2009 The Broad Institute
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package net.sf.samtools;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;

import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.samtools.util.BlockCompressedOutputStream;
import net.sf.samtools.util.IOUtil;
import net.sf.samtools.util.Md5CalculatingOutputStream;
import net.sf.samtools.util.RuntimeIOException;

/**
 * Create a SAMFileWriter for writing SAM or BAM.
 */
public class SAMFileWriterFactory {

	private static boolean DefaultCreateIndexWhileWriting = false;
	private boolean createIndex = DefaultCreateIndexWhileWriting;
	private static boolean defaultCreateMd5File = false;
	private boolean createMd5File = defaultCreateMd5File;

	private Integer maxRecordsInRam;

	/**
	 * Sets the default for whether to create md5Files for BAM files this
	 * factory.
	 */
	public static void setDefaultCreateMd5File(final boolean createMd5File) {
		defaultCreateMd5File = createMd5File;
	}

	/**
	 * Sets whether to create md5Files for BAMs from this factory.
	 */
	public SAMFileWriterFactory setCreateMd5File(final boolean createMd5File) {
		this.createMd5File = createMd5File;
		return this;
	}

	/**
	 * Sets the default for subsequent SAMFileWriterFactories that do not
	 * specify whether to create an index. If a BAM (not SAM) file is created,
	 * the setting is true, and the file header specifies coordinate order, then
	 * a BAM index file will be written along with the BAM file.
	 * 
	 * @param setting
	 *            whether to attempt to create a BAM index while creating the
	 *            BAM file
	 */
	public static void setDefaultCreateIndexWhileWriting(final boolean setting) {
		DefaultCreateIndexWhileWriting = setting;
	}

	/**
	 * Convenience method allowing
	 * newSAMFileWriterFactory().setCreateIndex(true); Equivalent to
	 * SAMFileWriterFactory.setDefaultCreateIndexWhileWriting(true);
	 * newSAMFileWriterFactory(); If a BAM (not SAM) file is created, the
	 * setting is true, and the file header specifies coordinate order, then a
	 * BAM index file will be written along with the BAM file.
	 * 
	 * @param setting
	 *            whether to attempt to create a BAM index while creating the
	 *            BAM file.
	 * @return this factory object
	 */
	public SAMFileWriterFactory setCreateIndex(final boolean setting) {
		this.createIndex = setting;
		return this;
	}

	/**
	 * Before creating a writer that is not presorted, this method may be called
	 * in order to override the default number of SAMRecords stored in RAM
	 * before spilling to disk (c.f. SAMFileWriterImpl.MAX_RECORDS_IN_RAM). When
	 * writing very large sorted SAM files, you may need call this method in
	 * order to avoid running out of file handles. The RAM available to the JVM
	 * may need to be increased in order to hold the specified number of records
	 * in RAM. This value affects the number of records stored in subsequent
	 * calls to one of the make...() methods.
	 * 
	 * @param maxRecordsInRam
	 *            Number of records to store in RAM before spilling to temporary
	 *            file when creating a sorted SAM or BAM file.
	 */
	public SAMFileWriterFactory setMaxRecordsInRam(int maxRecordsInRam) {
		this.maxRecordsInRam = maxRecordsInRam;
		return this;
	}

	/**
	 * Create a BAMFileWriter that is ready to receive SAMRecords. Uses default
	 * compression level.
	 * 
	 * @param header
	 *            entire header. Sort order is determined by the sortOrder
	 *            property of this arg.
	 * @param presorted
	 *            if true, SAMRecords must be added to the SAMFileWriter in
	 *            order that agrees with header.sortOrder.
	 * @param outputFile
	 *            where to write the output.
	 */
	public SAMFileWriter makeBAMWriter(final SAMFileHeader header, final boolean presorted, final File outputFile) {
		return makeBAMWriter(header, presorted, outputFile, BlockCompressedOutputStream.getDefaultCompressionLevel());
	}

	/**
	 * 
	 * Create a BAMFileWriter that is ready to receive SAMRecords.
	 * 
	 * @param header
	 *            entire header. Sort order is determined by the sortOrder
	 *            property of this arg.
	 * @param presorted
	 *            if true, SAMRecords must be added to the SAMFileWriter in
	 *            order that agrees with header.sortOrder.
	 * @param outputFile
	 *            where to write the output.
	 * @param compressionLevel
	 *            Override default compression level with the given value,
	 *            between 0 (fastest) and 9 (smallest).
	 */
	public SAMFileWriter makeBAMWriter(final SAMFileHeader header, final boolean presorted, final File outputFile,
			final int compressionLevel) {
		try {
			boolean createMd5File = this.createMd5File && IOUtil.isRegularPath(outputFile);
			if (this.createMd5File && !createMd5File) {
				System.err.println("Cannot create MD5 file for BAM because output file is not a regular file: "
						+ outputFile.getAbsolutePath());
			}
			final BAMFileWriter ret = createMd5File ? new BAMFileWriter(new Md5CalculatingOutputStream(
					new FileOutputStream(outputFile, false), new File(outputFile.getAbsolutePath() + ".md5")),
					outputFile, compressionLevel) : new BAMFileWriter(outputFile, compressionLevel);
			boolean createIndex = this.createIndex && IOUtil.isRegularPath(outputFile);
			if (this.createIndex && !createIndex) {
				System.err.println("Cannot create index for BAM because output file is not a regular file: "
						+ outputFile.getAbsolutePath());
			}
			initializeBAMWriter(ret, header, presorted, createIndex);
			return ret;
		} catch (IOException ioe) {
			throw new RuntimeIOException("Error opening file: " + outputFile.getAbsolutePath());
		}
	}

	private void initializeBAMWriter(BAMFileWriter writer, SAMFileHeader header, boolean presorted, boolean createIndex) {
		writer.setSortOrder(header.getSortOrder(), presorted);
		writer.setHeader(header);
		if (createIndex && writer.getSortOrder().equals(SAMFileHeader.SortOrder.coordinate)) {
			writer.enableBamIndexConstruction();
		}
		if (maxRecordsInRam != null) {
			writer.setMaxRecordsInRam(maxRecordsInRam);
		}
	}

	/**
	 * Create a SAMTextWriter that is ready to receive SAMRecords.
	 * 
	 * @param header
	 *            entire header. Sort order is determined by the sortOrder
	 *            property of this arg.
	 * @param presorted
	 *            if true, SAMRecords must be added to the SAMFileWriter in
	 *            order that agrees with header.sortOrder.
	 * @param outputFile
	 *            where to write the output.
	 */
	public SAMFileWriter makeSAMWriter(final SAMFileHeader header, final boolean presorted, final File outputFile) {
		try {
			final SAMTextWriter ret = this.createMd5File ? new SAMTextWriter(new Md5CalculatingOutputStream(
					new FileOutputStream(outputFile, false), new File(outputFile.getAbsolutePath() + ".md5")))
					: new SAMTextWriter(outputFile);
			ret.setSortOrder(header.getSortOrder(), presorted);
			if (maxRecordsInRam != null) {
				ret.setMaxRecordsInRam(maxRecordsInRam);
			}
			ret.setHeader(header);
			return ret;
		} catch (IOException ioe) {
			throw new RuntimeIOException("Error opening file: " + outputFile.getAbsolutePath());
		}
	}

	/**
	 * Create a SAMTextWriter for writing to a stream that is ready to receive
	 * SAMRecords. This method does not support the creation of an MD5 file
	 * 
	 * @param header
	 *            entire header. Sort order is determined by the sortOrder
	 *            property of this arg.
	 * @param presorted
	 *            if true, SAMRecords must be added to the SAMFileWriter in
	 *            order that agrees with header.sortOrder.
	 * @param stream
	 *            the stream to write records to.
	 */
	public SAMFileWriter makeSAMWriter(final SAMFileHeader header, final boolean presorted, final OutputStream stream) {
		final SAMTextWriter ret = new SAMTextWriter(stream);
		ret.setSortOrder(header.getSortOrder(), presorted);
		if (maxRecordsInRam != null) {
			ret.setMaxRecordsInRam(maxRecordsInRam);
		}
		ret.setHeader(header);
		return ret;
	}

	public SAMFileWriter makeCRAMWriter(ReferenceSequenceFile rsf, final SAMFileHeader header, final File outputFile) {
		try {
			CRAMFileWriter writer = new CRAMFileWriter(rsf, new BufferedOutputStream(new FileOutputStream(outputFile)),
					header);
			boolean createIndex = this.createIndex && IOUtil.isRegularPath(outputFile);
			if (this.createIndex && !createIndex) {
				System.err.println("Cannot create index for BAM because output file is not a regular file: "
						+ outputFile.getAbsolutePath());
			}
			return writer;
		} catch (IOException ioe) {
			throw new RuntimeIOException("Error opening file: " + outputFile.getAbsolutePath());
		}
	}

	/**
	 * Create either a SAM or a BAM writer based on examination of the
	 * outputFile extension.
	 * 
	 * @param header
	 *            entire header. Sort order is determined by the sortOrder
	 *            property of this arg.
	 * @param presorted
	 *            presorted if true, SAMRecords must be added to the
	 *            SAMFileWriter in order that agrees with header.sortOrder.
	 * @param outputFile
	 *            where to write the output. Must end with .sam or .bam.
	 * @return SAM or BAM writer based on file extension of outputFile.
	 */
	public SAMFileWriter makeSAMOrBAMWriter(final SAMFileHeader header, final boolean presorted, final File outputFile) {
		final String filename = outputFile.getName();
		if (filename.endsWith(".cram")) {
			ReferenceSequenceFile rsf = ReferenceDiscovery.findReferenceSequenceFileOrFail(outputFile);
			if (rsf == null)
				throw new RuntimeException("Reference fasta file is required for CRAM compression.");
			return makeCRAMWriter(rsf, header, outputFile);
		}
		if (filename.endsWith(".bam")) {
			return makeBAMWriter(header, presorted, outputFile);
		}
		if (filename.endsWith(".sam")) {
			return makeSAMWriter(header, presorted, outputFile);
		}
		return makeBAMWriter(header, presorted, outputFile);
	}

}

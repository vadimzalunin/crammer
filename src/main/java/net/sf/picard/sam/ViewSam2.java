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

package net.sf.picard.sam;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.PositionalArguments;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
import net.sf.samtools.Defaults;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

/**
 * A modified version of ViewSam which can output alignment slices. Most likely
 * this class will disappear soon.
 * 
 * @author vadim@ebi.ac.uk
 * 
 *         Very simple command that just reads a SAM or BAM file and writes out
 *         the header and each records to standard out.
 * 
 * @author tfennell@broad.mit.edu
 */
public class ViewSam2 extends CommandLineProgram {
	public static enum AlignmentStatus {
		Aligned, Unaligned, All
	}

	public static enum PfStatus {
		PF, NonPF, All
	}

	@Usage
	public final String USAGE = getStandardUsagePreamble() + "Prints a SAM or BAM file to the screen.";
	@Option(shortName = StandardOptionDefinitions.INPUT_SHORT_NAME, doc = "The SAM or BAM file to view.")
	public File INPUT;

	@Option(doc = "Print out all reads, just the aligned reads or just the unaligned reads.")
	public AlignmentStatus ALIGNMENT_STATUS = AlignmentStatus.All;

	@Option(doc = "Print out all reads, just the PF reads or just the non-PF reads.")
	public PfStatus PF_STATUS = PfStatus.All;

	@Option(doc = "Whether to exclude records overlapping alignment slice borders. ", optional = true)
	public boolean contained = false;

	@PositionalArguments()
	public List<String> queries;

	public static void main(final String[] args) {
		// hack to avoid 'java.io.IOException: Not enough storage is available
		// to process this command':
		// Defaults.BUFFER_SIZE = 1024 * 32;
		System.setProperty("samjdk.buffer_size", String.valueOf(1024 * 32));
		new ViewSam2().instanceMain(args);
	}

	@Override
	protected int doWork() {
		if (queries == null || queries.isEmpty()) {
			queries = new ArrayList<String>(1);
			queries.add(null);
		} else
			for (String q : queries)
				System.out.println(q);

		IoUtil.assertFileIsReadable(INPUT);
		File index = new File(INPUT + ".bai");
		if (!index.exists())
			index = new File(INPUT + ".crai");
		if (!index.exists())
			index = null;

		final SAMFileReader in = new SAMFileReader(INPUT, index);
		final SAMFileHeader header = in.getFileHeader();

		final SAMFileWriter out = new SAMFileWriterFactory().makeSAMWriter(header, true, System.out);

		for (String query : queries) {
			SAMRecordIterator iterator = null;
			if (query == null)
				iterator = in.iterator();
			else {
				AlignmentSliceQuery asq = null;
				try {
					asq = new AlignmentSliceQuery(query);
				} catch (Exception e) {
					throw new RuntimeException("Malformed alignment query: " + query);
				}
				iterator = in.query(asq.sequence, asq.start, asq.end, contained);
			}
			while (iterator.hasNext()) {
				SAMRecord rec = iterator.next();
				if (System.out.checkError()) {
					return 0;
				}

				if (this.ALIGNMENT_STATUS == AlignmentStatus.Aligned && rec.getReadUnmappedFlag())
					continue;
				if (this.ALIGNMENT_STATUS == AlignmentStatus.Unaligned && !rec.getReadUnmappedFlag())
					continue;

				if (this.PF_STATUS == PfStatus.PF && rec.getReadFailsVendorQualityCheckFlag())
					continue;
				if (this.PF_STATUS == PfStatus.NonPF && !rec.getReadFailsVendorQualityCheckFlag())
					continue;

				out.addAlignment(rec);
			}
		}
		out.close();

		return 0;
	}

}

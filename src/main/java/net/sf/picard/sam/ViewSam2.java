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
import java.util.List;

import net.sf.picard.cmdline.CommandLineProgram;
import net.sf.picard.cmdline.Option;
import net.sf.picard.cmdline.PositionalArguments;
import net.sf.picard.cmdline.StandardOptionDefinitions;
import net.sf.picard.cmdline.Usage;
import net.sf.picard.io.IoUtil;
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
		new ViewSam2().instanceMain(args);
	}

	@Override
	protected int doWork() {
		if (queries != null)
			for (String q : queries)
				System.out.println(q);

		IoUtil.assertFileIsReadable(INPUT);
		final SAMFileReader in = new SAMFileReader(INPUT);
		final SAMFileHeader header = in.getFileHeader();

		final SAMFileWriter out = new SAMFileWriterFactory().makeSAMWriter(header, true, System.out);

		for (String query : queries) {
			SAMRecordIterator iterator = null;
			if (query == null)
				iterator = in.iterator();
			else {
				AlignmentSliceQuery asq = null;
				try {
					asq = createAlignmentSliceQuery(query);
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

	private static AlignmentSliceQuery createAlignmentSliceQuery(String spec) {
		String[] chunks = spec.split(":");

		AlignmentSliceQuery q = new AlignmentSliceQuery();
		q.sequence = chunks[0];

		if (chunks.length > 1) {
			chunks = chunks[1].split("-");
			q.start = Integer.valueOf(chunks[0]);
			if (chunks.length == 2)
				q.end = Integer.valueOf(chunks[1]);
		}

		return q;
	}

	private static class AlignmentSliceQuery {
		private String sequence;
		private int start;
		private int end;
	}
}

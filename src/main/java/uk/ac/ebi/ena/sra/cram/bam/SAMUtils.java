package uk.ac.ebi.ena.sra.cram.bam;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.StringWriter;
import java.io.Writer;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.RuntimeIOException;
import uk.ac.ebi.ena.sra.cram.CramException;

public class SAMUtils {
	private static final String FIELD_SEPARATOR = "\t";

	/**
	 * Copied from net.sf.samtools.SAMTextWriter Attributes are commented
	 * because they are package-visible.
	 * 
	 * @param out
	 * @param alignment
	 */
	public static final void writeAlignment(Writer out,
			final SAMRecord alignment) {
		try {
			out.write(alignment.getReadName());
			out.write(FIELD_SEPARATOR);
			out.write(Integer.toString(alignment.getFlags()));
			out.write(FIELD_SEPARATOR);
			out.write(alignment.getReferenceName());
			out.write(FIELD_SEPARATOR);
			out.write(Integer.toString(alignment.getAlignmentStart()));
			out.write(FIELD_SEPARATOR);
			out.write(Integer.toString(alignment.getMappingQuality()));
			out.write(FIELD_SEPARATOR);
			out.write(alignment.getCigarString());
			out.write(FIELD_SEPARATOR);

			// == is OK here because these strings are interned
			if (alignment.getReferenceName() == alignment
					.getMateReferenceName()
					&& SAMRecord.NO_ALIGNMENT_REFERENCE_NAME != alignment
							.getReferenceName()) {
				out.write("=");
			} else {
				out.write(alignment.getMateReferenceName());
			}
			out.write(FIELD_SEPARATOR);
			out.write(Integer.toString(alignment.getMateAlignmentStart()));
			out.write(FIELD_SEPARATOR);
			out.write(Integer.toString(alignment.getInferredInsertSize()));
			out.write(FIELD_SEPARATOR);
			out.write(alignment.getReadString());
			out.write(FIELD_SEPARATOR);
			out.write(alignment.getBaseQualityString());
			// SAMBinaryTagAndValue attribute = alignment.getBinaryAttributes();
			// while (attribute != null) {
			// out.write(FIELD_SEPARATOR);
			// out.write(tagCodec.encode(tagUtil.makeStringTag(attribute.tag),
			// attribute.value));
			// attribute = attribute.getNext();
			// }
			out.write("\n");

		} catch (IOException e) {
			throw new RuntimeIOException(e);
		}
	}

	public static final String toString(SAMRecord record) {
		StringWriter writer = new StringWriter();
		writeAlignment(writer, record);
		return writer.toString();
	}

	public static IndexedFastaSequenceFile createIndexedFastaSequenceFile(
			File file) throws CramException, FileNotFoundException {
		if (IndexedFastaSequenceFile.canCreateIndexedFastaReader(file)) {
			IndexedFastaSequenceFile ifsFile = new IndexedFastaSequenceFile(
					file);

			return ifsFile;
		} else
			throw new CramException(
					"Reference fasta file is not indexed or index file not found. Try executing 'samtools faidx "
							+ file.getAbsolutePath() + "'");
	}
}

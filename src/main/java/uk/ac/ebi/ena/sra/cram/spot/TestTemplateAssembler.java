package uk.ac.ebi.ena.sra.cram.spot;

import java.io.File;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

class TestTemplateAssembler {

	public static void main(String[] args) {
		SAMFileReader
				.setDefaultValidationStringency(ValidationStringency.SILENT);
		SAMFileReader reader = new SAMFileReader(new File(
				"y:/Data/SangerExample/paired/5120_1.bam"));

		TemplateAssembler ta = new TemplateAssembler(1000, 1000);

		int maxRecords = 1000;
		SAMRecordIterator iterator = reader.iterator();
		while (iterator.hasNext()) {
			SAMRecord samRecord = iterator.next();
			// if
			// (!"IL24_5120:1:29:15064:14308#0".equals(samRecord.getReadName()))
			// continue;

			ta.addSAMRecord(samRecord);
			SAMSpot spot = ta.getNextAssembledTemplate();
			if (spot != null) {
				System.out.printf("%s", spot.getName());
				for (SAMRecordHolder holder : spot) {
					System.out.printf("\t%s (%d)", holder == null ? "?" : holder
							.getSamRecord().getAlignmentStart(), holder.getIndex());
				}
				System.out.println();
			}
			if (maxRecords-- <= 0)
				break;
		}

		SAMSpot spot = null;
		while ((spot = ta.getNextTemplate()) != null) {
			System.out.printf("%s", spot.getName());
			for (SAMRecordHolder holder : spot) {
				System.out.printf("\t%s (%d)", holder == null ? "?" : holder
						.getSamRecord().getAlignmentStart(), holder == null ? -1 : holder.getIndex());
			}
			System.out.println();
		}
	}
}

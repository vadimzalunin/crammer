package uk.ac.ebi.ena.sra.cram;

import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.reference.ReferenceSequenceFileFactory;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceRecord;

import org.apache.commons.math.random.EmpiricalDistribution;
import org.apache.commons.math.random.EmpiricalDistributionImpl;
import org.apache.commons.math.stat.Frequency;
import org.apache.commons.math.stat.descriptive.SummaryStatistics;

import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.format.SubstitutionVariation;
import uk.ac.ebi.ena.sra.cram.impl.CramRecordStaticFactory;

import com.beust.jcommander.JCommander;
import com.beust.jcommander.Parameter;
import com.beust.jcommander.ParameterException;
import com.beust.jcommander.Parameters;
import com.beust.jcommander.converters.FileConverter;

public class BamStats {

	public static void usage(JCommander jc) {
		StringBuilder sb = new StringBuilder();
		sb.append("\n");
		jc.usage(sb);

		System.out.println(sb.toString());
	}

	public static void main(String[] args) throws IOException {
		Params params = new Params();
		JCommander jc = new JCommander(params);
		try {
			jc.parse(args);
		} catch (ParameterException e) {
			e.printStackTrace();
			return;
		}

		if (args.length < 1 || params.help) {
			usage(jc);
			return;
		}

		ReferenceSequenceFile referenceSequenceFile = ReferenceSequenceFileFactory
				.getReferenceSequenceFile(new File(params.reference));

		SAMFileReader
				.setDefaultValidationStringency(ValidationStringency.SILENT);

		SAMFileReader reader = null;
		if (params.bamFile != null) {
			reader = new SAMFileReader(params.bamFile);
		} else {
			usage(jc);
			return;
			// URL url = new URL(params.bamLocation);
			//
			// if (params.localIndexFile == null) {
			// URL indexURL = new URL(params.bamLocation + ".bai");
			// InputStream indexFileStream = indexURL.openStream();
			// File indexFile = File.createTempFile("bai", null);
			// indexFile.deleteOnExit();
			// FileOutputStream indexFOS = new FileOutputStream(indexFile);
			// IOUtils.copy(indexFileStream, indexFOS);
			// indexFOS.close();
			// reader = new SAMFileReader(url, indexFile, false);
			// } else
			// reader = new SAMFileReader(url, params.localIndexFile, false);

		}

		if (!reader.hasIndex())
			System.out.println("Index: not found");

		SAMFileHeader header = reader.getFileHeader();
		System.out.println("Sort order: " + header.getSortOrder().name());

		System.out.println("Sequences: "
				+ header.getSequenceDictionary().getSequences().size());

		List<String> seqNames = new ArrayList<String>();
		for (SAMSequenceRecord seq : header.getSequenceDictionary()
				.getSequences()) {
			seqNames.add(seq.getSequenceName());
			if (params.maxSequences > 0
					&& seqNames.size() >= params.maxSequences)
				break;
		}

		long recordCounter = 0L;
		long mappedCounter = 0L;
		// long perfectMatchCounter = 0L;
		long insCounter = 0L;
		long delCounter = 0L;
		long subsCounter = 0L;
		long otherCigarOpsCounter = 0L;
		long negStrandCounter = 0L;

		CramRecord prevCramRecord = null;
		SAMRecord prevSamRecord = null;

		Frequency alignmentStartDistanceFrequency = new Frequency();

		Frequency insLenFreq = new Frequency();
		Frequency delLenFreq = new Frequency();

		Frequency insPosFreq = new Frequency();
		Frequency delPosFreq = new Frequency();
		Frequency subsPosFreq = new Frequency();

		Frequency subsBaseFreq = new Frequency();
		Frequency subsRefBaseFreq = new Frequency();

		Frequency cigarOpsFreq = new Frequency();

		Frequency mappedQualFreq = new Frequency();
		Frequency unmappedQualFreq = new Frequency();

		Frequency inferredInsertFreq = new Frequency();

		EmpiricalDistributionImpl inferredInsertDistr = new EmpiricalDistributionImpl(
				10);

		for (String seqName : seqNames) {
			int pos = 0;

			SAMRecordIterator iterator = reader.queryOverlapping(seqName, 0, 0);
			if (iterator.hasNext())
				System.out.println("Reading reads for sequence: " + seqName);

			ReferenceSequence sequence = referenceSequenceFile
					.getSequence(seqName);
			byte[] refSequence = referenceSequenceFile.getSubsequenceAt(
					seqName, 0, sequence.length()).getBases();
			while (iterator.hasNext()) {
				SAMRecord record = iterator.next();

				if (recordCounter++ > params.maxRecordsPerSequence
						&& params.maxRecordsPerSequence > 0)
					break;

				if (record.getReadNegativeStrandFlag())
					negStrandCounter++;

				if (record.getReadUnmappedFlag()) {
					for (byte q : record.getBaseQualities())
						unmappedQualFreq.addValue(q);
					continue;
				}

				for (byte q : record.getBaseQualities())
					mappedQualFreq.addValue(q);

				mappedCounter++;
				alignmentStartDistanceFrequency.addValue(record
						.getAlignmentStart() - pos);
				pos = record.getAlignmentStart();

				Cigar cigar = record.getCigar();

				int posInRead = 1;
				for (CigarElement ce : cigar.getCigarElements()) {
					cigarOpsFreq.addValue(ce.getOperator().ordinal());
					switch (ce.getOperator()) {
					case D:
						delCounter++;
						delLenFreq.addValue(ce.getLength());
						delPosFreq.addValue(posInRead);
						break;
					case I:
						insCounter++;
						insLenFreq.addValue(ce.getLength());
						insPosFreq.addValue(posInRead);
						break;
					case M:

						break;
					default:
						otherCigarOpsCounter++;
						break;
					}
					posInRead += ce.getLength();
				}

				CramRecord cramRecord = CramRecordStaticFactory.newRecord3(record,
						refSequence);

				if (cramRecord.getSubstitutionVariations() != null) {
					for (SubstitutionVariation sv : cramRecord
							.getSubstitutionVariations()) {
						subsBaseFreq.addValue(sv.getBase());
						subsCounter++;
						subsPosFreq.addValue(sv.getPosition());
						subsRefBaseFreq.addValue(sv.getRefernceBase());
					}
				}

				if (prevSamRecord == null)
					inferredInsertFreq.addValue(0);
				else
					inferredInsertFreq.addValue(record.getInferredInsertSize()
							- prevSamRecord.getInferredInsertSize());

				prevCramRecord = cramRecord;
				prevSamRecord = record;

			}
			iterator.close();
			System.out.println("Record counter: " + recordCounter);
		}

		System.out.println("recordCounter: " + recordCounter);
		System.out.println("mappedCounter: " + mappedCounter);
		System.out.println("insCounter: " + insCounter);
		System.out.println("delCounter: " + delCounter);
		System.out.println("subsCounter: " + subsCounter);
		System.out.println("otherCigarOpsCounter: " + otherCigarOpsCounter);
		System.out.println("negStrandCounter: " + negStrandCounter);

		System.out.println("alignmentStartDistanceFrequency: ");
		System.out.println(alignmentStartDistanceFrequency);

		System.out.println("insLenFreq: ");
		System.out.println(insLenFreq);

		System.out.println("delLenFreq: ");
		System.out.println(delLenFreq);

		System.out.println("insPosFreq: ");
		System.out.println(insPosFreq);

		System.out.println("delPosFreq: ");
		System.out.println(delPosFreq);

		System.out.println("subsPosFreq: ");
		System.out.println(subsPosFreq);

		System.out.println("cigarOpsFreq: ");
		System.out.println(cigarOpsFreq);

		System.out.println("Cigar legend:");
		for (CigarOperator op : CigarOperator.values())
			System.out.println(op.ordinal() + "=" + op.name());

		System.out.println("mappedQualFreq: ");
		System.out.println(mappedQualFreq);

		System.out.println("unmappedQualFreq: ");
		System.out.println(unmappedQualFreq);

		System.out.println("inferredInsertFreq: ");
		System.out.println(frequencyToString(inferredInsertFreq, 0,
				params.maxStatsLines));

		System.out.println("subsRefBaseFreq: ");
		System.out.println(subsRefBaseFreq);

		System.out.println("subsBaseFreq: ");
		System.out.println(subsBaseFreq.toString());

		List<Double> dlist = new ArrayList<Double>();
		Iterator<Comparable<?>> iter = inferredInsertFreq.valuesIterator();
		while (iter.hasNext()) {
			long value = (Long) iter.next();
			if (value > -1)
				dlist.add(new Double(value));
		}
		double[] darray = new double[dlist.size()];
		for (int i = 0; i < darray.length; i++)
			darray[i] = dlist.get(i);

		inferredInsertDistr.load(darray);
		System.out.println("InferredInsertSize distribution: (min/max/count)");
		System.out.println(distributionToString(inferredInsertDistr));
	}

	private static String frequencyToString(Frequency frequency, long min,
			long max) {
		NumberFormat nf = NumberFormat.getPercentInstance();
		StringBuffer outBuffer = new StringBuffer();
		outBuffer.append("Value \t Freq. \t Pct. \t Cum Pct. \n");
		Iterator<Comparable<?>> iter = frequency.valuesIterator();
		while (iter.hasNext()) {
			Comparable<?> value = iter.next();
			if (min > (Long) value || max < ((Long) value))
				continue;
			outBuffer.append(value);
			outBuffer.append('\t');
			outBuffer.append(frequency.getCount(value));
			outBuffer.append('\t');
			outBuffer.append(nf.format(frequency.getPct(value)));
			outBuffer.append('\t');
			outBuffer.append(nf.format(frequency.getCumPct(value)));
			outBuffer.append('\n');
		}
		return outBuffer.toString();
	}

	private static String distributionToString(EmpiricalDistribution distr) {
		StringBuffer outBuffer = new StringBuffer();
		outBuffer.append("Min \t Max \t N \n");
		for (SummaryStatistics stat : distr.getBinStats())
			if (stat.getN() > 0)
				outBuffer.append(String.format("%.2f\t%.2f\t%d\n",
						stat.getMin(), stat.getMax(), stat.getN()));
		return outBuffer.toString();
	}

	@Parameters(commandDescription = "BAM file statistics")
	static class Params {
		@Parameter(names = { "--bam-file" }, converter = FileConverter.class)
		File bamFile;

		@Parameter(names = { "--index-file" }, converter = FileConverter.class)
		File localIndexFile;

		// @Parameter(names = { "--bam-url" })
		// String bamLocation;

		@Parameter(names = { "--max-sequences" })
		int maxSequences = 0;

		@Parameter(names = { "--max-records-per-sequence" })
		long maxRecordsPerSequence = 0;

		@Parameter(names = { "--reference" })
		String reference;

		@Parameter(names = { "--help", "-h" })
		boolean help = false;

		@Parameter(names = { "--max-stats-lines" })
		int maxStatsLines = 1000;
	}
}

package uk.ac.ebi.ena.sra.cram.bam;

import static org.hamcrest.CoreMatchers.equalTo;
import static org.junit.Assert.assertThat;

import java.io.ByteArrayInputStream;
import java.io.IOException;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.StringUtil;

import org.junit.Test;

import uk.ac.ebi.ena.sra.cram.SequenceBaseProvider;
import uk.ac.ebi.ena.sra.cram.bam.Sam2CramRecordFactory.TREAT_TYPE;
import uk.ac.ebi.ena.sra.cram.format.CramRecord;
import uk.ac.ebi.ena.sra.cram.impl.ByteArraySequenceBaseProvider;
import uk.ac.ebi.ena.sra.cram.impl.ReadFeatures2Cigar;
import uk.ac.ebi.ena.sra.cram.impl.RestoreBases;

public class Sam2CramRecordTest {

	/**
	 * Creates a SAMRecord out of its string representation. Replaces spaces
	 * with a tab.
	 * 
	 * @param samRecordString
	 *            SAM-style string
	 * @param referenceBases
	 *            exact bases to which the record should map
	 * @return SAMRecord object
	 */
	private static SAMRecord createSAMRecord(String samRecordString, byte[] referenceBases) {
		String[] words = samRecordString.split("\\s+");
		if (words.length < 11)
			throw new IllegalArgumentException("Not enough words in the SAM record: " + samRecordString);
		String seqName = words[2];
		String correctedSAMRecordString = StringUtil.join("\t", words);

		ByteArrayInputStream bais = new ByteArrayInputStream(correctedSAMRecordString.getBytes());
		SAMFileReader reader = new SAMFileReader(bais);
		SAMSequenceRecord samSequenceRecord = new SAMSequenceRecord(seqName, referenceBases.length);
		reader.getFileHeader().addSequence(samSequenceRecord);
		SAMRecord samRecord = reader.iterator().next();
		return samRecord;
	}

	@Test
	public void test_Reversed_read() throws IOException {
		String samRecordString = "ERR005143.135209	16	gi|66043271|ref|NC_007005.1|	1	37	4M1I30M2D1M	*	0	0	CCCGCCCACGGCTGCTCCAGCTGCTGCTGTAGCGAC	JXR>>WQNPhOhUdh`hhhhhchhhhhhhhhhhhhh	XT:A:U	NM:i:4	X0:i:1	X1:i:0	XM:i:1	XO:i:1	XG:i:1	MD:Z:2T31^CG1";
		String refString = "CCTGCCACGGCTGCTCCAGCTGCTGCTGTAGCGACGC";
		byte[] referenceBases = refString.getBytes();

		SAMRecord samRecord = createSAMRecord(samRecordString, referenceBases);

		Sam2CramRecordFactory factory = new Sam2CramRecordFactory(refString.getBytes());
		CramRecord cramRecord = factory.createCramRecord(samRecord);
		SequenceBaseProvider provider = new ByteArraySequenceBaseProvider(refString.getBytes());
		byte[] restoredBases = new RestoreBases(provider, samRecord.getReferenceName()).restoreReadBases(cramRecord);

		assertThat(new String(restoredBases), equalTo(new String(samRecord.getReadBases())));
	}

	@Test
	public void test1() throws IOException {
		String samRecordString = "ERR005143.532371	0	gi|66043271|ref|NC_007005.1|	1	25	36M	*	0	0	GTGTCAGTGGAACTTTGGCCGCCGTGCGTGGAGCTT	hhhhhhhhhhMbhehhTaYBNJ@]LN\\INDO>WDNf	XT:A:U	NM:i:2	X0:i:1	X1:i:0	XM:i:2	XO:i:0	XG:i:0	MD:Z:19A2A13";
		String refString = "GTGTCAGTGGAACTTTGGCAGCAGTGCGTGGAGCTT";
		byte[] referenceBases = refString.getBytes();

		SAMRecord samRecord = createSAMRecord(samRecordString, referenceBases);

		Sam2CramRecordFactory factory = new Sam2CramRecordFactory(refString.getBytes());
		CramRecord cramRecord = factory.createCramRecord(samRecord);
		SequenceBaseProvider provider = new ByteArraySequenceBaseProvider(refString.getBytes());
		byte[] restoredBases = new RestoreBases(provider, samRecord.getReferenceName()).restoreReadBases(cramRecord);

		assertThat(new String(restoredBases), equalTo(new String(samRecord.getReadBases())));
	}

	@Test
	public void test2() throws IOException {
		String samRecordString = "ERR005143.1517506	16	gi|66043271|ref|NC_007005.1|	1	37	31M1D5M	*	0	0	TCGACATCAAGCTCGATGACTCGAAACTGATGTCGA	Rhghhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh	XT:A:U	NM:i:2	X0:i:1	X1:i:0	XM:i:1	XO:i:1	XG:i:1	MD:Z:31^C3A1";
		String refString = "TCGACATCAAGCTCGATGACTCGAAACTGATCGTCAA";
		byte[] referenceBases = refString.getBytes();

		SAMRecord samRecord = createSAMRecord(samRecordString, referenceBases);

		Sam2CramRecordFactory factory = new Sam2CramRecordFactory(refString.getBytes());

		CramRecord cramRecord = factory.createCramRecord(samRecord);
		SequenceBaseProvider provider = new ByteArraySequenceBaseProvider(refString.getBytes());
		byte[] restoredBases = new RestoreBases(provider, samRecord.getReferenceName()).restoreReadBases(cramRecord);

		assertThat(new String(restoredBases), equalTo(new String(samRecord.getReadBases())));
	}

	@Test
	public void test3() throws IOException {
		// System.out.println(new String(Utils.getBasesFromReferenceFile(
		// "y:/Data/psyringae/psyringae.fa",
		// "gi|66043271|ref|NC_007005.1|", 486970-1, 40)));

		String samRecordString = "ERR005143.1967430       0       gi|66043271|ref|NC_007005.1|    2  36      1S2M2D33M       *       0       0       GTCCAAACAGCGCCAGCAGGTTGGGATCGCGCTGTT    hhhhhhhh[hhhhfS[hLhZ[PhS\\IOSKKTHPKP\\    XT:A:U  NM:i:3     X0:i:1  X1:i:0  XM:i:1  XO:i:1  XG:i:1  MD:Z:2^GG1C31";
		String refString = "CTCGGCCAACAGCGCCAGCAGGTTGGGATCGCGCTGTTTGG";
		String ss = "GTC  CAAACAGCGCCAGCAGGTTGGGATCGCGCTGTT";
		byte[] referenceBases = refString.getBytes();

		SAMRecord samRecord = createSAMRecord(samRecordString, referenceBases);

		Sam2CramRecordFactory factory = new Sam2CramRecordFactory(refString.getBytes());
		factory.setTreatSoftClipsAs(TREAT_TYPE.INSERTION);
		CramRecord cramRecord = factory.createCramRecord(samRecord);
		SequenceBaseProvider provider = new ByteArraySequenceBaseProvider(refString.getBytes());
		byte[] restoredBases = new RestoreBases(provider, samRecord.getReferenceName()).restoreReadBases(cramRecord);

		assertThat(new String(restoredBases), equalTo(new String(samRecord.getReadBases())));
	}

	@Test
	public void test4() throws IOException {
		// System.out.println(new String(Utils.getBasesFromReferenceFile(
		// "y:/Data/psyringae/psyringae.fa",
		// "gi|66043271|ref|NC_007005.1|", 675923, 40)));

		String samRecordString = "ERR005143.2606307       16      gi|66043271|ref|NC_007005.1|    1  37      32M2D3M1S       *       0       0       GAGCCCTGGCATGTCCTGGGTGAAGAGGGCTCCAGC    hhhh`_hhhXhhhhhhfhhhhhhhhhhhhhhhhhhh    XT:A:U  NM:i:3     X0:i:1  X1:i:0  XM:i:1  XO:i:1  XG:i:1  MD:Z:32^GG1C1";
		String refString = "GAGCCCTGGCATGTCCTGGGTGAAGAGGGCTCGGCCGGCGG";
		byte[] referenceBases = refString.getBytes();

		SAMRecord samRecord = createSAMRecord(samRecordString, referenceBases);

		Sam2CramRecordFactory factory = new Sam2CramRecordFactory(refString.getBytes());
		factory.setTreatSoftClipsAs(TREAT_TYPE.INSERTION);
		CramRecord cramRecord = factory.createCramRecord(samRecord);
		SequenceBaseProvider provider = new ByteArraySequenceBaseProvider(refString.getBytes());
		byte[] restoredBases = new RestoreBases(provider, samRecord.getReferenceName()).restoreReadBases(cramRecord);

		assertThat(new String(restoredBases), equalTo(new String(samRecord.getReadBases())));
	}

	@Test
	public void test5() throws IOException {
		// System.out.println(new String(Utils.getBasesFromReferenceFile(
		// "y:/Data/psyringae/psyringae.fa",
		// "gi|66043271|ref|NC_007005.1|", 343, 40)));

		String samRecordString = "ERR005143.135209        16      gi|66043271|ref|NC_007005.1|    1     37      4M1I30M2D1M     *       0       0       CCCGCCCACGGCTGCTCCAGCTGCTGCTGTAGCGAC    JXR>>WQNPhOhUdh`hhhhhchhhhhhhhhhhhhh    XT:A:U  NM:i:4     X0:i:1  X1:i:0  XM:i:1  XO:i:1  XG:i:1  MD:Z:2T31^CG1";
		String refString = "CCTGCCACGGCTGCTCCAGCTGCTGCTGTAGCGACGCCTGC";
		byte[] referenceBases = refString.getBytes();

		SAMRecord samRecord = createSAMRecord(samRecordString, referenceBases);

		Sam2CramRecordFactory factory = new Sam2CramRecordFactory(refString.getBytes());
		CramRecord cramRecord = factory.createCramRecord(samRecord);
		SequenceBaseProvider provider = new ByteArraySequenceBaseProvider(refString.getBytes());
		byte[] restoredBases = new RestoreBases(provider, samRecord.getReferenceName()).restoreReadBases(cramRecord);

		assertThat(new String(restoredBases), equalTo(new String(samRecord.getReadBases())));
	}

	@Test
	public void test6() throws IOException {
		// System.out.println(new String(Utils.getBasesFromReferenceFile(
		// "y:/Data/psyringae/psyringae.fa",
		// "gi|66043271|ref|NC_007005.1|", 343, 40)));

		String samRecordString = "HS3_6007:1:2103:12124:168409#3	99	MAL1	1	29	3S37M1I7M1D23M4S	=	1	595	CTGAACCCTGAACCCTAAACCCTAAACCCTAAACCCTAAACCCCTAAACCTAAACCCTAAACCCTGAACCTTCTA	HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHFHHHGHHGGHHHHFHHGHHEFFHHHEEHHFHHGHFHEHFEGE6	MD:Z:6A6G6G9A13^C23	XG:i:2	AM:i:29	NM:i:6	SM:i:29	XM:i:4	XO:i:2	QT:Z:HHHHHHHF	RT:Z:TTAGGCAT	XT:A:M";
		String refString = "CTGAACCCTAAACCCTGAACCCTGAACCCTAAAACCTAAACCCCTAACACCTAAACCCTAAACCCTGAACCTTCTA";
		byte[] referenceBases = refString.getBytes();

		SAMRecord samRecord = createSAMRecord(samRecordString, referenceBases);

		Sam2CramRecordFactory factory = new Sam2CramRecordFactory(refString.getBytes());
		factory.setTreatSoftClipsAs(TREAT_TYPE.INSERTION);
		CramRecord cramRecord = factory.createCramRecord(samRecord);
		SequenceBaseProvider provider = new ByteArraySequenceBaseProvider(refString.getBytes());
		byte[] restoredBases = new RestoreBases(provider, samRecord.getReferenceName()).restoreReadBases(cramRecord);

		assertThat(new String(restoredBases), equalTo(new String(samRecord.getReadBases())));

		ReadFeatures2Cigar rf2Cigar = new ReadFeatures2Cigar();

		assertThat(rf2Cigar.getCigar2(cramRecord.getReadFeatures(), (int) cramRecord.getReadLength()),
				equalTo(samRecord.getCigar()));
	}
}

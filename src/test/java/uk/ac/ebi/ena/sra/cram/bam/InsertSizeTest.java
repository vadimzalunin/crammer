package uk.ac.ebi.ena.sra.cram.bam;


import java.io.BufferedOutputStream;
import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.PrintStream;

import net.sf.picard.sam.SamPairUtil;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMRecordUtil;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.CoordMath;
import net.sf.samtools.util.SequenceUtil;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Ignore;
import org.junit.Test;

import uk.ac.ebi.ena.sra.cram.Utils;

public class InsertSizeTest {

	@BeforeClass
	public static void setUpBeforeClass() throws Exception {
	}

	@AfterClass
	public static void tearDownAfterClass() throws Exception {
	}

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	@Ignore
	public void test1 () {
		ByteArrayOutputStream baos = new ByteArrayOutputStream(1024*1024) ;
		PrintStream ps = new PrintStream(baos) ;
		ps.println("@HD	VN:1.0	SO:coordinate")	;
		ps.println("@SQ	SN:20	LN:135374737")	;
		ps.println("SOLEXA4_1:1:88:906:645	99	20	61545	0	50M	=	61689	194	ATCACAATTAAAAGAACTAGAAAAGCAAGAGCAAACACATTCAAAAGCTA	HHHHIHHHHHHHHHHHHHHHGGHHHGHHHHHHEHHHHHHHHGIGHHHEHH	XT:A:R	NM:i:0	SM:i:0	AM:i:0	X0:i:2705	XM:i:0	XO:i:0	XG:i:0	MD:Z:50")	;
		ps.println("SOLEXA4_1:1:88:906:645	147	20	61689	0	50M	=	61545	-194	TTTAAAAGGACAAACAAAATTGATAGACCGCTAGCAAGACTAATAAAGAA	>A5>A545-1+-551BFFBF8</18B:HCHEHHHDGHHHHHGHBHGHCHH	XT:A:R	NM:i:2	SM:i:0	AM:i:0	X0:i:16	X1:i:1471	XM:i:2	XO:i:0	XG:i:0	MD:Z:3G7C38")	;
		ps.println("6	99	20	61545	0	50M	=	61689	193	ATCACAATTAAAAGAACTAGAAAAGCAAGAGCAAACACATTCAAAAGCTA	HHHHIHHHHHHHHHHHHHHHGGHHHGHHHHHHEHHHHHHHHGIGHHHEHH	MQ:i:0")	;
		ps.println("6	147	20	61689	0	50M	=	61545	-193	TTTAAAAGGACAAACAAAATTGATAGACCGCTAGCAAGACTAATAAAGAA	>A5>A545-1+-551BFFBF8</18B:HCHEHHHDGHHHHHGHBHGHCHH	MQ:i:0")	;
		
		ps.close() ;
		ByteArrayInputStream bais = new ByteArrayInputStream(baos.toByteArray()) ;
		
		SAMFileReader reader = new SAMFileReader(bais) ;
		SAMRecordIterator iterator = reader.iterator() ;
		
		SAMRecord record1 = iterator.next() ;
		System.out.println(record1.getSAMString());
		SAMRecord record2 = iterator.next() ;
		System.out.println(record2.getSAMString());
		System.out.println(SamPairUtil.computeInsertSize(record1, record2)) ;
		
		SAMRecord record3 = iterator.next() ;
		System.out.println(record3.getSAMString());
		SAMRecord record4 = iterator.next() ;
		System.out.println(record4.getSAMString());
		
		System.out.println(SamPairUtil.computeInsertSize(record3, record4)) ;
	}
	
	@Test
	@Ignore
	public void test2 () {
		ByteArrayOutputStream baos = new ByteArrayOutputStream(1024*1024) ;
		PrintStream ps = new PrintStream(baos) ;
		ps.println("@HD	VN:1.0	SO:coordinate")	;
		ps.println("@SQ	SN:20	LN:135374737")	;
		ps.println("SOLEXA7_30:3:14:1078:1926	83	20	23902	60	50M	=	23852	-100	CTTATCTCTTAATATTGAAACTCACTAGAATTTAATTCTAGTCCTCTTTT	@>AA9BABA?BC9+5>8@B?CB?3@AAA8CBB@@84:BBB;*;B7BBB?9	XT:A:U	NM:i:0	SM:i:37	AM:i:37	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:50")	;
		ps.println("SOLEXA7_30:3:14:1078:1926	163	20	23852	60	50M	=	23902	100	ATCTGACTTATCTTATTGGTCCTTTTAAGTCGTTTCCTCTTATCTGATCT	ACBCC:>A=ACBCC@@BAB=BBBB9@=;@=<A<><*7@9A=7==>;-<7@	XT:A:U	NM:i:0	SM:i:37	AM:i:37	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:50")	;
		ps.println("1	83	20	23902	60	50M	=	23852	-98	CTTATCTCTTAATATTGAAACTCACTAGAATTTAATTCTAGTCCTCTTTT	@>AA9BABA?BC9+5>8@B?CB?3@AAA8CBB@@84:BBB;*;B7BBB?9	MQ:i:60")	;
		ps.println("1	163	20	23852	60	50M	=	23902	100	ATCTGACTTATCTTATTGGTCCTTTTAAGTCGTTTCCTCTTATCTGATCT	ACBCC:>A=ACBCC@@BAB=BBBB9@=;@=<A<><*7@9A=7==>;-<7@	MQ:i:60")	;
		
		ps.close() ;
		ByteArrayInputStream bais = new ByteArrayInputStream(baos.toByteArray()) ;
		
		SAMFileReader reader = new SAMFileReader(bais) ;
		SAMRecordIterator iterator = reader.iterator() ;
		
		SAMRecord record1 = iterator.next() ;
		SAMRecord record2 = iterator.next() ;
		System.out.println(Utils.computeInsertSize(record1, record2)) ;
		System.out.println(Utils.computeInsertSize(record2, record1)) ;
		
		SAMRecord record3 = iterator.next() ;
		SAMRecord record4 = iterator.next() ;

		Utils.setLooseMateInfo(record3, record4, reader.getFileHeader()) ;
		System.out.println(record3.getSAMString());
		System.out.println(record4.getSAMString());
		
		Utils.setLooseMateInfo(record4, record3, reader.getFileHeader()) ;
		System.out.println(record3.getSAMString());
		System.out.println(record4.getSAMString());
	}
	@Test
	public void test3 () {
		ByteArrayOutputStream baos = new ByteArrayOutputStream(1024*1024) ;
		PrintStream ps = new PrintStream(baos) ;
		ps.println("@HD	VN:1.0	SO:coordinate")	;
		ps.println("@SQ	SN:20	LN:135374737")	;
		ps.println("SOLEXA8_49:3:81:832:1907	163	20	199704	60	50M	=	199864	204	GAAAGAAAGAGGGAGGAAGGGCAGAGGGAGCAGGGAGACTGTAGATCAGG	BAA@B@BBB>BBB?BB=@BB@@?A<AA@6??6???6?76<>55;5751<6	XT:A:U	NM:i:0	SM:i:37	AM:i:37	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:50")	;
		ps.println("SOLEXA8_49:3:81:832:1907	83	20	199864	60	45M5S	=	199704	-204	TCTTCCACAGGAATGTTGAGGATGACNTCCATGTCTGGGGTGCACTTGGG	%%%%%:3455043+286<88<54464&3;??@A?<@B@A@B@@8?;9@BB	XT:A:U	NM:i:1	SM:i:37	AM:i:37	X0:i:1	X1:i:0	XM:i:0	XO:i:0	XG:i:0	MD:Z:26A18	ZJ:Z:G_199909_204608")	;
		
		ps.close() ;
		ByteArrayInputStream bais = new ByteArrayInputStream(baos.toByteArray()) ;
		
		SAMFileReader reader = new SAMFileReader(bais) ;
		SAMRecordIterator iterator = reader.iterator() ;
		
		SAMRecord record1 = iterator.next() ;
		SAMRecord record2 = iterator.next() ;
		Utils.setLooseMateInfo(record1, record2, reader.getFileHeader()) ;
		System.out.println(SamPairUtil.computeInsertSize(record2, record1));
		System.out.println(record1.getSAMString());
		System.out.println(record2.getSAMString());
		
		Utils.setLooseMateInfo(record2, record1, reader.getFileHeader()) ;
		System.out.println(record1.getSAMString());
		System.out.println(record2.getSAMString());
	}
}

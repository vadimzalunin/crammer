package uk.ac.ebi.ena.sra.cram;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;

import uk.ac.ebi.ena.sra.cram.encoding.BaseChangeCodecTest;
import uk.ac.ebi.ena.sra.cram.encoding.BaseSequenceCodecTest;
import uk.ac.ebi.ena.sra.cram.encoding.GammaCodecTest;
import uk.ac.ebi.ena.sra.cram.encoding.GolombCodecTest;
import uk.ac.ebi.ena.sra.cram.encoding.GolombRiceCodecTest;
import uk.ac.ebi.ena.sra.cram.encoding.HuffmanCodecTest;
import uk.ac.ebi.ena.sra.cram.encoding.UnaryCodecTest;
import uk.ac.ebi.ena.sra.cram.format.text.TRAMRoundTripTests;
import uk.ac.ebi.ena.sra.cram.impl.CramRecordCodecRoundTripTests;
import uk.ac.ebi.ena.sra.cram.io.BitInputStreamTest;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStreamTest;
import uk.ac.ebi.ena.sra.cram.io.BitStreamIntegrationTest;
import uk.ac.ebi.ena.sra.cram.mask.FastaByteArrayMaskFactoryTest;
import uk.ac.ebi.ena.sra.cram.mask.IntegerListMaskFactoryTest;
import uk.ac.ebi.ena.sra.cram.mask.SingleLineMaskReaderTest;

@RunWith(Suite.class)
@Suite.SuiteClasses({ BaseChangeCodecTest.class, BaseSequenceCodecTest.class,
		UnaryCodecTest.class, GammaCodecTest.class, GolombRiceCodecTest.class,
		GolombCodecTest.class, HuffmanCodecTest.class,
		BitInputStreamTest.class, BitOutputStreamTest.class,
		BitStreamIntegrationTest.class, BamRoundTripTests.class,
		FastaByteArrayMaskFactoryTest.class, IntegerListMaskFactoryTest.class,
		SingleLineMaskReaderTest.class, /* text formating needs re-design: CramRecordTextFormatTest.class,
		DefaultReadFeatureFormatTest.class,*/
		CramRecordCodecRoundTripTests.class, TRAMRoundTripTests.class })
public class TestAll {

}

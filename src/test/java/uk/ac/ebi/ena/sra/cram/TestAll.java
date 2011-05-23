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
import uk.ac.ebi.ena.sra.cram.impl.CramRecordWriteReadTests;
import uk.ac.ebi.ena.sra.cram.io.BitInputStreamTest;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStreamTest;
import uk.ac.ebi.ena.sra.cram.io.BitStreamIntegrationTest;

@RunWith(Suite.class)
@Suite.SuiteClasses({ BaseChangeCodecTest.class, BaseSequenceCodecTest.class,
		UnaryCodecTest.class, GammaCodecTest.class, GolombRiceCodecTest.class,
		GolombCodecTest.class, HuffmanCodecTest.class,
		CramRecordWriteReadTests.class, BitInputStreamTest.class,
		BitOutputStreamTest.class, BitStreamIntegrationTest.class })
public class TestAll {

}

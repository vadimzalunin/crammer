package uk.ac.ebi.ena.sra.cram;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;

import uk.ac.ebi.ena.sra.cram.encoding.BaseChangeCodecTest;
import uk.ac.ebi.ena.sra.cram.encoding.BaseSequenceCodecTest;
import uk.ac.ebi.ena.sra.cram.encoding.GolombRiceCodecTest;
import uk.ac.ebi.ena.sra.cram.encoding.UnaryCodecTest;
import uk.ac.ebi.ena.sra.cram.impl.CramRecordWriteReadTests;
import uk.ac.ebi.ena.sra.cram.io.BitInputStreamTest;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStreamTest;
import uk.ac.ebi.ena.sra.cram.io.BitStreamIntegrationTest;

@RunWith(Suite.class)
@Suite.SuiteClasses({ BaseChangeCodecTest.class, BaseSequenceCodecTest.class,
		GolombRiceCodecTest.class, UnaryCodecTest.class,
		CramRecordWriteReadTests.class, BitInputStreamTest.class,
		BitOutputStreamTest.class, BitStreamIntegrationTest.class})
public class TestAll {

}

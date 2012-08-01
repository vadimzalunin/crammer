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
package uk.ac.ebi.ena.sra.cram;

import org.junit.runner.RunWith;
import org.junit.runners.Suite;

import uk.ac.ebi.ena.sra.cram.bam.InsertSizeTest;
import uk.ac.ebi.ena.sra.cram.bam.Sam2CramRecordTest;
import uk.ac.ebi.ena.sra.cram.encoding.BaseChangeCodecTest;
import uk.ac.ebi.ena.sra.cram.encoding.BaseSequenceCodecTest;
import uk.ac.ebi.ena.sra.cram.encoding.GammaCodecTest;
import uk.ac.ebi.ena.sra.cram.encoding.GolombCodecTest;
import uk.ac.ebi.ena.sra.cram.encoding.GolombRiceCodecTest;
import uk.ac.ebi.ena.sra.cram.encoding.HuffmanByteCodec2Test;
import uk.ac.ebi.ena.sra.cram.encoding.HuffmanByteCodecTest;
import uk.ac.ebi.ena.sra.cram.encoding.HuffmanCodecTest;
import uk.ac.ebi.ena.sra.cram.encoding.ReadFeatureCodecTest;
import uk.ac.ebi.ena.sra.cram.encoding.SubexpCodecTest;
import uk.ac.ebi.ena.sra.cram.encoding.SubstitutionVariationCodecTest;
import uk.ac.ebi.ena.sra.cram.encoding.UnaryCodecTest;
import uk.ac.ebi.ena.sra.cram.format.compression.NumberCodecFactoryTest;
import uk.ac.ebi.ena.sra.cram.format.text.TRAMRoundTripTests;
import uk.ac.ebi.ena.sra.cram.impl.CramRecordCodecRoundTripTests;
import uk.ac.ebi.ena.sra.cram.impl.ReadFeatures2CigarTest;
import uk.ac.ebi.ena.sra.cram.impl.RestoreBasesTest;
import uk.ac.ebi.ena.sra.cram.impl.TestCramIterators;
import uk.ac.ebi.ena.sra.cram.io.BitInputStreamTest;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStreamTest;
import uk.ac.ebi.ena.sra.cram.io.BitStreamIntegrationTest;
import uk.ac.ebi.ena.sra.cram.io.LongBufferBitInputStreamTest;
import uk.ac.ebi.ena.sra.cram.io.LongBufferBitStreamIntegrationTest;
import uk.ac.ebi.ena.sra.cram.mask.FastaByteArrayMaskFactoryTest;
import uk.ac.ebi.ena.sra.cram.mask.IntegerListMaskFactoryTest;
import uk.ac.ebi.ena.sra.cram.mask.SingleLineMaskReaderTest;

@RunWith(Suite.class)
// @formatter:off
@Suite.SuiteClasses({
	BamRoundTripTests.class,
	BaseChangeCodecTest.class,
	BaseSequenceCodecTest.class, 
	BitInputStreamTest.class,
	BitOutputStreamTest.class, 
	BitStreamIntegrationTest.class, 
	CramRecordCodecRoundTripTests.class, 
	FastaByteArrayMaskFactoryTest.class, 
	GammaCodecTest.class,
	GolombCodecTest.class, 
	GolombRiceCodecTest.class, 
	HuffmanByteCodec2Test.class, 
	HuffmanByteCodecTest.class, 
	HuffmanCodecTest.class,
	InsertSizeTest.class, 
	IntegerListMaskFactoryTest.class,
	LongBufferBitInputStreamTest.class,
	LongBufferBitStreamIntegrationTest.class,
	NumberCodecFactoryTest.class,
	ReadFeatureCodecTest.class, 
	ReadFeatures2CigarTest.class, 
	RestoreBasesTest.class, 
	Sam2CramRecordTest.class,
	SingleLineMaskReaderTest.class,
	SubexpCodecTest.class, 
	SubstitutionVariationCodecTest.class, 
	TRAMRoundTripTests.class, 
	TestCramIterators.class,
	TestLongJumps.class, 
	UnaryCodecTest.class
})
// @formatter:on
public class TestAll {

}

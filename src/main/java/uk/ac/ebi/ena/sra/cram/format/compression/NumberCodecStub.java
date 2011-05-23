package uk.ac.ebi.ena.sra.cram.format.compression;

import uk.ac.ebi.ena.sra.cram.encoding.BitCodec;

public interface NumberCodecStub extends BitCodec<Long>{

	public EncodingAlgorithm getEncoding () ;
	public String getStringRepresentation () ;
	public void initFromString (String spec) throws CramCompressionException ;
	public Object[] getParameters () ;
	public void setParameters (Object [] params) ;
}

package uk.ac.ebi.ena.sra.cram.encoding;

import java.io.IOException;

import uk.ac.ebi.ena.sra.cram.io.BitInputStream;
import uk.ac.ebi.ena.sra.cram.io.BitOutputStream;

public interface ByteArrayBitCodec {
	
	public String getName () ;

	public byte[] read(BitInputStream bis, int len) throws IOException;

	public long write(BitOutputStream bos, byte[] object) throws IOException;

	public long numberOfBits(byte[] object);
	
	public Stats getStats () ;
	
	public static class Stats {
		public long nofBis;
		public long arraysRead, arraysWritten;
		public long bytesRead, bytesWritten;
		
		public void reset() {
			nofBis = 0 ;
			arraysRead = 0 ;
			arraysWritten = 0 ;
			bytesRead = 0 ;
			bytesWritten = 0 ;
		}
	}
}

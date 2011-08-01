package uk.ac.ebi.ena.sra.cram.mask;

import java.io.IOException;

public interface ReadMaskReader {

	public PositionMask readNextMask() throws IOException, ReadMaskFormatException;

}
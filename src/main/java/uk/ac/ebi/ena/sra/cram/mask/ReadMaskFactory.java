package uk.ac.ebi.ena.sra.cram.mask;

public interface ReadMaskFactory<T> {

	public PositionMask createMask(T data) throws ReadMaskFormatException;
}

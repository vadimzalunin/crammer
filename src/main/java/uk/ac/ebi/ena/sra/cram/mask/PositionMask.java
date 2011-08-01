package uk.ac.ebi.ena.sra.cram.mask;

public interface PositionMask {

	public boolean isMasked (int position) ;
	public int[] getMaskedPositions () ;
	public boolean isEmpty () ;
	public int getMaskedCount() ;
	public int getMinMaskedPosition () ;
	public int getMaxMaskedPosition () ;
	public byte[] toByteArrayUsing (byte mask, byte nonMask) ;
}

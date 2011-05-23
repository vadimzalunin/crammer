package uk.ac.ebi.ena.sra.cram.format;

public interface ReadFeature {
	public static final byte STOP_OPERATOR = '$' ;

	/**
	 * @return zero-based position in the read
	 */
	public int getPosition();

	public void setPosition(int position);

	public byte getOperator();
}

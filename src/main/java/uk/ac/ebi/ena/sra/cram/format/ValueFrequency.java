package uk.ac.ebi.ena.sra.cram.format;

@Deprecated
public class ValueFrequency<T> {

	private T value;
	private int frequency;

	public ValueFrequency(T value, int frequency) {
		if (value == null)
			throw new NullPointerException("Value cannot be null.");
		this.value = value;
		this.frequency = frequency;
	}

	public int getFrequency() {
		return frequency;
	}

	public void setFrequency(int frequency) {
		this.frequency = frequency;
	}

	public T getValue() {
		return value;
	}

	public void setValue(T value) {
		this.value = value;
	}

	@Override
	public boolean equals(Object obj) {
		if (!(obj instanceof ValueFrequency))
			return false;
		ValueFrequency foe = (ValueFrequency) obj;
		return value.equals(foe.value) && frequency == foe.frequency;
	}
}

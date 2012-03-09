package uk.ac.ebi.ena.sra.cram.format;

import java.util.Arrays;
import java.util.Iterator;
import java.util.Set;
import java.util.SortedSet;

import org.apache.commons.collections.bag.HashBag;

public class IntFrequencies {
	private HashBag bag = new HashBag();

	public IntFrequencies() {
	}

	public IntFrequencies(int[] values, int[] freqs) {
		for (int i = 0; i < values.length; i++)
			bag.add(values[i], freqs[i]);
	}

	public void add(int value) {
		bag.add(value);
	}

	public void add(int value, int howMany) {
		bag.add(value, howMany);
	}

	public void add(int[] bytes) {
		for (int b : bytes)
			bag.add(b);
	}

	public int getFrequency(int value) {
		return bag.getCount(value);
	}

	public int[] getValues() {
		Set<?> uniqueSet = bag.uniqueSet() ;
		int[] collapsedValueArray = new int[uniqueSet.size()];

		Iterator<?> iterator = uniqueSet.iterator();
		int valueIndex = 0;
		while (iterator.hasNext())
			collapsedValueArray[valueIndex++] = ((Integer) iterator.next()).intValue();

		Arrays.sort(collapsedValueArray) ;
		return collapsedValueArray;
	}

	public int[] getFrequencies() {
		int[] values = getValues () ;
		int[] collapsedFrequencyArray = new int[values.length] ;
		int frequencyIndex = 0;
		for (int value:values)
			collapsedFrequencyArray[frequencyIndex++] = bag.getCount(value);

		return collapsedFrequencyArray;
	}

	@Override
	public String toString() {
		return bag.toString();
	}
}

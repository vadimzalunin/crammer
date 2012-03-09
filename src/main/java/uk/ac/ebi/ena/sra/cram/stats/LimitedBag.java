package uk.ac.ebi.ena.sra.cram.stats;

import java.util.Collection;
import java.util.Iterator;
import java.util.Set;

import org.apache.commons.collections.Bag;
import org.apache.commons.collections.bag.HashBag;
import org.apache.commons.math.stat.Frequency;

public class LimitedBag implements Bag {
	public static final int MAX_DISTINCT_VALUES = 1000;

	private String key;
	private Bag bag;
	private int maxDistinctValues;

	public LimitedBag(String key, int maxDistinctValues, Bag bag) {
		if (key == null)
			throw new NullPointerException("Key is null.");
		if (bag == null)
			throw new NullPointerException("Bag is null.");
		if (maxDistinctValues < 1)
			throw new IllegalArgumentException("Max distinct values must be more than one.");

		this.key = key;
		this.maxDistinctValues = maxDistinctValues;
		this.bag = bag;
	}

	public LimitedBag(String key, int maxDistinctValues) {
		this(key, maxDistinctValues, new HashBag());
	}

	public LimitedBag(String key) {
		this(key, MAX_DISTINCT_VALUES, new HashBag());
	}
	
	public boolean add(Object value) {
		if (isFull())
			return true;
		return bag.add(value);
	}

	public boolean isFull() {
		return bag.size() >= maxDistinctValues;
	}

	public String getKey() {
		return key;
	}

	public Bag getValueBag() {
		return bag;
	}

	public int getCount(Object object) {
		return bag.getCount(object);
	}

	public boolean add(Object object, int nCopies) {
		if (isFull())
			return true;
		return add(object, nCopies);
	}

	public boolean remove(Object object) {
		return bag.remove(object);
	}

	public boolean remove(Object object, int nCopies) {
		return bag.remove(object, nCopies);
	}

	public Set uniqueSet() {
		return bag.uniqueSet();
	}

	public int size() {
		return bag.size();
	}

	public boolean containsAll(Collection coll) {
		return bag.containsAll(coll);
	}

	public boolean isEmpty() {
		return bag.isEmpty();
	}

	public boolean contains(Object o) {
		return bag.contains(o);
	}

	public boolean removeAll(Collection coll) {
		return bag.removeAll(coll);
	}

	public boolean retainAll(Collection coll) {
		return bag.retainAll(coll);
	}

	public Object[] toArray() {
		return bag.toArray();
	}

	public Object[] toArray(Object[] a) {
		return bag.toArray(a);
	}

	public Iterator iterator() {
		return bag.iterator();
	}

	public boolean addAll(Collection c) {
		return bag.addAll(c);
	}

	public void clear() {
		bag.clear();
	}
}

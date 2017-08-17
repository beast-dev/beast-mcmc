/*
 * FrequencySet.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.util;

import java.util.*;

/**
 * Stores a set of objects with frequencies
 *
 * @version $Id: FrequencySet.java,v 1.9 2005/06/27 21:18:40 rambaut Exp $
 *
 * @author Andrew Rambaut
 */
public class FrequencySet<T>
{
	//
	// Public stuff
	//

	public FrequencySet() {}

	/** get number of objects */
	public int size()
	{
		return size;
	}

	/** get object in frequency order */
	public T get(int i) {
		if (!sorted) {
			sortByFrequency();
		}

		return list.get(i).object;
	}


    protected int getFrequency(T b) {
        Bin bin = binMap.get(b);
        if( bin == null ) {
            return 0;
        }
        return bin.frequency;
    }

    /** get frequency of ith object */
	public int getFrequency(int i)
	{
		if (!sorted) {
			sortByFrequency();
		}

		return list.get(i).frequency;
	}

	/** get sum of all frequencies */
	public int getSumFrequency() {
		int sum = 0;
		for (int i = 0, n = size(); i < n; i++) {
			sum += getFrequency(i);
		}

		return sum;
	}

	/** adds an object to the set */
	public void add(T object)
	{
		add(object, 1);
	}

	/**
     * adds an object to the set with an initial frequency, or if object already
     * in frequency set then frequency is incremented by 1. 
     */
	public void add(T object, int frequency) {

		Bin bin = binMap.get(object);
		if (bin != null) {
			bin.frequency += 1;
		} else {
			bin = new Bin(object, frequency);
			binMap.put(object, bin);
			size += 1;
			sorted = false;
		}
	}

	public Set<T> getKeySet() {
		return binMap.keySet();
	}

	/** The frequencySets are equal if their inner sets are equal */
	@Override
    public boolean equals(Object obj) {
		return (obj instanceof FrequencySet) && binMap.keySet().equals(((FrequencySet)obj).binMap.keySet());
	}

	/** sort by descending frequency */
	private void sortByFrequency() {

		list.clear();
		for (Bin bin : binMap.values()) {
			list.add(bin);
		}

		Collections.sort(list, frequencyComparator);
		sorted = true;
	}

	//
	// Private stuff
	//

	private List<Bin> list = new ArrayList<Bin>();
	private Map<T, Bin> binMap = new HashMap<T, Bin>();
	private boolean sorted = false;
	private int size = 0;

	private class Bin {
		T object;
        int frequency;

		public Bin(T object, int frequency) {
			this.object = object;
			this.frequency = frequency;
		}

		@Override
        public boolean equals(Object obj) {
			return object.equals(((Bin)obj).object);
		}

		public int hashCode() {
			return object.hashCode();
		}
	}

	private Comparator<Bin> frequencyComparator = new Comparator<Bin>() {
		public int compare(Bin bin1, Bin bin2) {
			return bin2.frequency - bin1.frequency;
		}
	};
}

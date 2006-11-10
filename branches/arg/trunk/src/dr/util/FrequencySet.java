/*
 * FrequencySet.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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
public class FrequencySet
{
	//
	// Public stuff
	//

	public FrequencySet()
	{
	}

	/** get number of objects */
	public int size()
	{
		return size;
		//return list.size();

	}

	/** get object in frequency order */
	public Object get(int i)
	{

		//return null;

		if (!sorted) {
			sortByFrequency();
		}

		return ((Bin)list.get(i)).object;
	}

	/** get frequency of ith object */
	public int getFrequency(int i)
	{
		if (!sorted) {
			sortByFrequency();
		}

		return ((Bin)list.get(i)).frequency;
	}

	/** get sum of all frequencies */
	public int getSumFrequency()
	{
		int sum = 0;
		for (int i = 0, n = size(); i < n; i++) {
			sum += getFrequency(i);
		}

		return sum;
	}

	/** adds an object to the set */
	public void add(Object object)
	{
		add(object, 1);
	}

	/** adds an object to the set with an initial frequency */
	public void add(Object object, int frequency)
	{
/*		if (!set.contains(object)) {
			Bin bin = new Bin(object, frequency);

			list.add(bin);
			set.add(object);
		} else {
			for (int i = 0, n = list.size(); i < n; i++) {
				Bin bin = (Bin)list.get(i);
				if (bin.object.equals(object)) {
					bin.frequency++;
					break;
				}
			}
		}*/

		Bin bin = (Bin)hash.get(object.toString());
		if (bin != null) {
			bin.frequency += 1;
		} else {
			bin = new Bin(object, frequency);
			hash.put(object.toString(), bin);
			size += 1;
			sorted = false;
		}

		/*boolean found = false;

		for (int i = 0, n = list.size(); i < n; i++) {
			Bin bin = (Bin)list.get(i);
			if (bin.object.equals(object)) {
				bin.frequency++;
				found = true;
				break;
			}
		}

		if (!found) {
			Bin bin = new Bin(object, frequency);

			list.add(bin);
			set.add(object);
		}
		*/

	}

	/** The frequencySets are equal if their inner sets are equal */
	public boolean equals(Object obj) {
		return set.equals(((FrequencySet)obj).set);
	}

	/** sort by descending frequency */
	private void sortByFrequency() {

		list.clear();
		for (Iterator i = hash.values().iterator(); i.hasNext();) {
			list.add(i.next());
		}

		Collections.sort(list, frequencyComparator);
		sorted = true;
	}

	//
	// Private stuff
	//

	private ArrayList list = new ArrayList();
	private Hashtable hash = new Hashtable();
	private HashSet set = new HashSet();
	private boolean sorted = false;
	private int size = 0;

	private class Bin {
		Object object;
		int frequency;

		public Bin(Object object, int frequency) {
			this.object = object;
			this.frequency = frequency;
		}

		public boolean equals(Object obj) {
			return object.equals(((Bin)obj).object);
		}

		public int hashCode() {
			return object.hashCode();
		}
	}

	private Comparator frequencyComparator = new Comparator() {
		public int compare(Object o1, Object o2) {
			return ((Bin)o2).frequency - ((Bin)o1).frequency;
		}

	};
}

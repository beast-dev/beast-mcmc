/*
 * Permutator.java
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

import java.lang.reflect.Array;
import java.util.Iterator;
import java.util.NoSuchElementException;

public class Permutator implements Iterator {

	private final int size;
	private final Object [] elements;  // copy of original 0 .. size-1
	private final Object ar;           // array for output,  0 .. size-1
	private final int [] permutation;  // perm of nums 1..size, perm[0]=0

	private boolean next = true;

	// int[], double[] array won't work :-(
	public Permutator (Object [] e) {
		size = e.length;
		elements = new Object [size];    // not suitable for primitives
		System.arraycopy (e, 0, elements, 0, size);
		ar = Array.newInstance (e.getClass().getComponentType(), size);
		System.arraycopy (e, 0, ar, 0, size);
		permutation = new int [size+1];
		for (int i=0; i<size+1; i++) {
			permutation [i]=i;
		}
	}

	private void formNextPermutation () {
		for (int i=0; i<size; i++) {
		// i+1 because perm[0] always = 0
		// perm[]-1 because the numbers 1..size are being permuted
			Array.set (ar, i, elements[permutation[i+1]-1]);
		}
	}

	public boolean hasNext() {
		return next;
	}

	public void remove() throws UnsupportedOperationException {
		throw new UnsupportedOperationException();
	}

	private void swap (final int i, final int j) {
		final int x = permutation[i];
		permutation[i] = permutation [j];
		permutation[j] = x;
	}

	// does not throw NoSuchElement; it wraps around!
	public Object next() throws NoSuchElementException {

		formNextPermutation ();  // copy original elements

		int i = size-1;
		while (permutation[i]>permutation[i+1]) i--;

		if (i==0) {
			next = false;
			for (int j=0; j<size+1; j++) {
				permutation [j]=j;
			}
			return ar;
		}

		int j = size;

		while (permutation[i]>permutation[j]) j--;
		swap (i,j);
		int r = size;
		int s = i+1;
		while (r>s) { swap(r,s); r--; s++; }

		return ar;
	}

	public String toString () {
		final int n = Array.getLength(ar);
		final StringBuffer sb = new StringBuffer ("[");
		for (int j=0; j<n; j++) {
			sb.append (Array.get(ar,j).toString());
			if (j<n-1) sb.append (",");
		}
		sb.append("]");
		return new String (sb);
	}

	public static void main (String [] args) {
		for (Iterator i = new Permutator(args); i.hasNext(); ) {
			i.next();
			System.out.println (i);
		}
	}
}

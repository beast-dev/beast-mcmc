/*
 * SequenceList.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evolution.sequence;

import dr.evolution.util.TaxonList;

import java.util.List;

/**
 * Interface for a list of sequences.
 *
 *
 * @author Andrew Rambaut
 */
public interface SequenceList extends TaxonList {

	/**
	 * @return a count of the number of sequences in the list.
	 */
	public int getSequenceCount();

	/**
	 * @return the ith sequence in the list.
	 */
	public Sequence getSequence(int i);

	/**
	 * Sets an named attribute for a given sequence.
	 * @param index the index of the sequence whose attribute is being set.
	 * @param name the name of the attribute.
	 * @param value the new value of the attribute.
	 */
	public void setSequenceAttribute(int index, String name, Object value);

	/**
	 * @return an object representing the named attributed for the given sequence.
	 * @param index the index of the sequence whose attribute is being fetched.
	 * @param name the name of the attribute of interest.
	 */
	public Object getSequenceAttribute(int index, String name);

	/**
	 * @return an immutable iterable List of sequences
	 */
	public List<Sequence> getSequences();
}


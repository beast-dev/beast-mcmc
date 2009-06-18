/*
 * Sequences.java
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

package dr.evolution.sequence;

import dr.evolution.util.Taxon;
import dr.util.Attributable;
import dr.util.Identifiable;

import java.util.*;

/**
 * Class for storing sequences.
 *
 * @version $Id: Sequences.java,v 1.12 2005/05/24 20:25:56 rambaut Exp $
 *
 * @author Andrew Rambaut
 */
public class Sequences implements SequenceList, Attributable, Identifiable {

 	/**
	 * Constructor.
	 */
	public Sequences() {
	}

   // **************************************************************
    // SequenceList IMPLEMENTATION
    // **************************************************************

	/**
	 * @return a count of the number of sequences in the list.
	 */
	public int getSequenceCount() {
		return sequences.size();
	}

	/**
	 * @return the ith sequence in the list.
	 */
	public Sequence getSequence(int index) {
		return (Sequence)sequences.get(index);
	}

	/**
	 * Sets an named attribute for a given sequence.
	 * @param index the index of the sequence whose attribute is being set.
	 * @param name the name of the attribute.
	 * @param value the new value of the attribute.
	 */
	public void setSequenceAttribute(int index, String name, Object value) {
		Sequence sequence = getSequence(index);
		sequence.setAttribute(name, value);
	}

	/**
	 * @return an object representing the named attributed for the given sequence.
	 * @param index the index of the sequence whose attribute is being fetched.
	 * @param name the name of the attribute of interest.
	 */
	public Object getSequenceAttribute(int index, String name) {
		Sequence sequence = getSequence(index);
		return sequence.getAttribute(name);
	}

    // **************************************************************
    // TaxonList IMPLEMENTATION
    // **************************************************************

	/**
	 * @return a count of the number of taxa in the list.
	 */
	public int getTaxonCount() {
		return getSequenceCount();
	}

	/**
	 * @return the taxon of the ith sequence.
	 */
	public Taxon getTaxon(int taxonIndex) {
		return getSequence(taxonIndex).getTaxon();
	}

	/**
	 * @return the ID of the taxon of the ith sequence. If it doesn't have
	 * a taxon, returns the ID of the sequence itself.
	 */
	public String getTaxonId(int taxonIndex) {
		Taxon taxon = getTaxon(taxonIndex);
		if (taxon != null)
			return taxon.getId();
		else
			throw new IllegalArgumentException("Illegal taxon index:" + taxonIndex);
	}

	/**
	 * returns the index of the taxon with the given id.
	 */
	public int getTaxonIndex(String id) {
		for (int i = 0, n = getTaxonCount(); i < n; i++) {
			if (getTaxonId(i).equals(id)) return i;
		}
		return -1;
	}

	/**
	 * returns the index of the given taxon.
	 */
	public int getTaxonIndex(Taxon taxon) {
		for (int i = 0, n = getTaxonCount(); i < n; i++) {
			if (getTaxon(i) == taxon) return i;
		}
		return -1;
	}

    public List<Taxon> asList() {
        List<Taxon> taxa = new ArrayList<Taxon>();
        for (int i = 0, n = getTaxonCount(); i < n; i++) {
            taxa.add(getTaxon(i));
        }
        return taxa;
    }

    /**
	 * Sets an named attribute for the taxon of a given sequence. If the sequence
	 * doesn't have a taxon then the attribute is added to the sequence itself.
	 * @param taxonIndex the index of the taxon whose attribute is being set.
	 * @param name the name of the attribute.
	 * @param value the new value of the attribute.
	 */
	public void setTaxonAttribute(int taxonIndex, String name, Object value) {
		Taxon taxon = getTaxon(taxonIndex);
		if (taxon != null)
			taxon.setAttribute(name, value);
		else
			setSequenceAttribute(taxonIndex, name, value);
	}

	/**
	 * @return an object representing the named attributed for the given taxon.
	 * @param taxonIndex the index of the taxon whose attribute is being fetched.
	 * @param name the name of the attribute of interest.
	 */
	public Object getTaxonAttribute(int taxonIndex, String name) {
		Taxon taxon = getTaxon(taxonIndex);
		if (taxon != null)
			return taxon.getAttribute(name);
		else
			return getSequenceAttribute(taxonIndex, name);
	}

    public Iterator<Taxon> iterator() {
        return new Iterator<Taxon>() {
            private int index = -1;

            public boolean hasNext() {
                return index < getTaxonCount() - 1;
            }

            public Taxon next() {
                index ++;
                return getTaxon(index);
            }

            public void remove() { /* do nothing */ }
        };
    }

    // **************************************************************
    // Identifiable IMPLEMENTATION
    // **************************************************************

	protected String id = null;

	/**
	 * @return the id.
	 */
	public String getId() {
		return id;
	}

	/**
	 * Sets the id.
	 */
	public void setId(String id) {
		this.id = id;
	}

    // **************************************************************
    // Sequences METHODS
    // **************************************************************

	/**
	 * Add a sequence to the sequence list
	 */
	public void addSequence(Sequence sequence) {
		sequences.add(sequence);
	}

	/**
	 * Insert a sequence to the sequence list at position
	 */
	public void insertSequence(int position, Sequence sequence) {
		sequences.insertElementAt(sequence, position);
	}

	/**
	 * Remove a sequence from the sequence list
	 * @return the sequence.
	 */
	public Sequence removeSequence(int index) {
		Sequence sequence = getSequence(index);
		sequences.removeElementAt(index);

		return sequence;
	}

    // **************************************************************
    // Attributable IMPLEMENTATION
    // **************************************************************

	/**
	 * Sets an named attribute for this object.
	 * @param name the name of the attribute.
	 * @param value the new value of the attribute.
	 */
	public void setAttribute(String name, Object value) {
		if (attributes == null)
			attributes = new Attributable.AttributeHelper();
		attributes.setAttribute(name, value);
	}

	/**
	 * @return an object representing the named attributed for this object.
	 * @param name the name of the attribute of interest.
	 */
	public Object getAttribute(String name) {
		if (attributes == null)
			return null;
		else
			return attributes.getAttribute(name);
	}

	/**
	 * @return an iterator of the attributes that this object has.
	 */
	public Iterator getAttributeNames() {
		if (attributes == null)
			return null;
		else
			return attributes.getAttributeNames();
	}

	// **************************************************************
	// INSTANCE VARIABLES
	// **************************************************************

	private Vector sequences = new Vector();

	private Attributable.AttributeHelper attributes = null;
}

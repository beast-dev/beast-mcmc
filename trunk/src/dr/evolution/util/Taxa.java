/*
 * Taxa.java
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

package dr.evolution.util;

import dr.util.Identifiable;

import java.util.ArrayList;

/**
 * Class for a list of taxa.
 *
 * @version $Id: Taxa.java,v 1.29 2006/09/05 13:29:34 rambaut Exp $
 *
 * @author Alexei Drummond
 */
public class Taxa implements MutableTaxonList, Identifiable, Comparable {
	public Taxa() {
	}

	public Taxa(String id) {
		this.id = id;
	}

	/** Adds a taxon of the given name and return its index. */
	public int addTaxon(Taxon taxon) {
		taxa.add(taxon);
		fireTaxonAdded(taxon);
		return taxa.size() - 1;
	}

	/** Removes a taxon of the given name and returns true if successful. */
	public boolean removeTaxon(Taxon taxon) {
		boolean success = taxa.remove(taxon);
		if (success) {
			fireTaxonRemoved(taxon);
		}
		return success;
	}

	public void removeAllTaxa() {
		taxa.clear();
		fireTaxonRemoved(null);
	}

	/**
	 * @return a count of the number of taxa in the list.
	 */
	public int getTaxonCount() {
		return taxa.size();
	}

	/**
	 * @return the ith taxon.
	 */
	public Taxon getTaxon(int taxonIndex) {
		return taxa.get(taxonIndex);
	}

	/**
	 * @return the ID of the ith taxon.
	 */
	public String getTaxonId(int taxonIndex) {
		return (taxa.get(taxonIndex)).getId();
	}

	/**
	 * Sets the unique identifier of the ith taxon.
	 */
	public void setTaxonId(int taxonIndex, String id) {
		(taxa.get(taxonIndex)).setId(id);
		fireTaxaChanged();
	}

	/**
	 * returns the index of the taxon with the given id.
	 */
	public int getTaxonIndex(String id) {
		for (int i =0; i < taxa.size(); i++) {
			if (getTaxonId(i).equals(id)) return i;
		}
		return -1;
	}

	/**
	 * returns the index of the given taxon.
	 */
	public int getTaxonIndex(Taxon taxon) {
		for (int i =0; i < taxa.size(); i++) {
			if (getTaxon(i) == taxon) return i;
		}
		return -1;
	}

    /**
     * Returns true if at least 1 member of taxonList is contained in this Taxa.
     * @param taxonList a TaxonList
     * @return true if any of taxonList is in this Taxa
     */
    public boolean containsAny(TaxonList taxonList) {

        for (int i = 0; i < taxonList.getTaxonCount(); i++) {
            Taxon taxon = taxonList.getTaxon(i);
            if (taxa.contains(taxon)) {
                return true;
            }
        }
        return false;
    }

	/**
	 * Returns true if taxonList is a subset of the taxa in this Taxa.
	 * @param taxonList a TaxonList
     * @return true if all of taxonList is in this Taxa
	 */
	public boolean containsAll(TaxonList taxonList) {

		for (int i = 0; i < taxonList.getTaxonCount(); i++) {
			Taxon taxon = taxonList.getTaxon(i);
			if (!taxa.contains(taxon)) {
				return false;
			}
		}
		return true;
	}

	public String getId() { return id; }
	public void setId(String id) { this.id = id; }

	public int compareTo(Object o) {
		return getId().compareTo(((Taxa)o).getId());
	}

	public String toString() { return id; }

	private String id = null;

	/**
	 * Sets an named attribute for a given taxon.
	 * @param taxonIndex the index of the taxon whose attribute is being set.
	 * @param name the name of the attribute.
	 * @param value the new value of the attribute.
	 */
	public void setTaxonAttribute(int taxonIndex, String name, Object value) {
		Taxon taxon = getTaxon(taxonIndex);
		taxon.setAttribute(name, value);
		fireTaxaChanged();
	}

	/**
	 * @return an object representing the named attributed for the given taxon.
	 * @param taxonIndex the index of the taxon whose attribute is being fetched.
	 * @param name the name of the attribute of interest.
	 */
	public Object getTaxonAttribute(int taxonIndex, String name) {
		Taxon taxon = getTaxon(taxonIndex);
		return taxon.getAttribute(name);
	}

	public void addMutableTaxonListListener(MutableTaxonListListener listener) {
		mutableTaxonListListeners.add(listener);
	}

	private void fireTaxonAdded(Taxon taxon) {
        for (MutableTaxonListListener mutableTaxonListListener : mutableTaxonListListeners) {
            mutableTaxonListListener.taxonAdded(this, taxon);
        }
    }

	private void fireTaxonRemoved(Taxon taxon) {
        for (MutableTaxonListListener mutableTaxonListListener : mutableTaxonListListeners) {
            mutableTaxonListListener.taxonRemoved(this, taxon);
        }
    }

	private void fireTaxaChanged() {
        for (MutableTaxonListListener mutableTaxonListListener : mutableTaxonListListeners) {
            mutableTaxonListListener.taxaChanged(this);
        }
    }

	private ArrayList<MutableTaxonListListener> mutableTaxonListListeners = new ArrayList<MutableTaxonListListener>();

	ArrayList<Taxon> taxa = new ArrayList<Taxon>();
	//ArrayList listeners = new ArrayList();
}

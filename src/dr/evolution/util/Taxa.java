/*
 * Taxa.java
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

package dr.evolution.util;

import dr.util.Identifiable;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

/**
 * Class for a list of taxa.
 *
 * @version $Id: Taxa.java,v 1.29 2006/09/05 13:29:34 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class Taxa implements MutableTaxonList, Identifiable, Comparable<Taxa> {

	private final ArrayList<MutableTaxonListListener> mutableTaxonListListeners = new ArrayList<MutableTaxonListListener>();
	ArrayList<Taxon> taxa = new ArrayList<Taxon>();

    private String id = null;

	public Taxa() {
	}

    public Taxa(String id) {
		this.id = id;
	}

    public Taxa(TaxonList list) {
        addTaxa(list);
    }

    public Taxa(Collection<Taxon> taxa) {
        addTaxa(taxa);
    }

    /**
     * Adds the given taxon and returns its index. If the taxon is already in the list then it is not
     * added but it returns its index.
     *
     * @param taxon the taxon to be added
     */
	public int addTaxon(Taxon taxon) {
        int index = getTaxonIndex(taxon);
        if (index == -1) {
            taxa.add(taxon);
            fireTaxonAdded(taxon);
            index = taxa.size() - 1;
        }
        return index;
    }

	/** Removes a taxon of the given name and returns true if successful. */
	public boolean removeTaxon(Taxon taxon) {
		boolean success = taxa.remove(taxon);
		if (success) {
			fireTaxonRemoved(taxon);
		}
		return success;
	}

    public void addTaxa(TaxonList taxa) {
        for(int nt = 0; nt < taxa.getTaxonCount(); ++nt) {
            addTaxon(taxa.getTaxon(nt));
        }
    }

    public void addTaxa(Collection<Taxon> taxa) {
        for(Taxon taxon : taxa) {
            addTaxon(taxon);
        }
    }

    public void removeTaxa(TaxonList taxa) {
        for(int nt = 0; nt < taxa.getTaxonCount(); ++nt) {
            removeTaxon(taxa.getTaxon(nt));
        }
    }

    public void removeTaxa(Collection<Taxon> taxa) {
        for(Taxon taxon : taxa) {
            removeTaxon(taxon);
        }
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
		for (int i = 0; i < taxa.size(); i++) {
			if (getTaxonId(i).equals(id)) return i;
		}
		return -1;
	}

	/**
	 * returns the index of the given taxon.
	 */
	public int getTaxonIndex(Taxon taxon) {
		for (int i = 0; i < taxa.size(); i++) {
			if (getTaxon(i) == taxon) return i;
		}
		return -1;
	}

    public List<Taxon> asList() {
        return new ArrayList<Taxon>(taxa);
    }

    /**
     * returns whether the list contains this taxon
     * @param taxon
     * @return true if taxon is in the list
     */
    public boolean contains(Taxon taxon) {
        return taxa.contains(taxon);
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
     * Returns true if at least 1 member of taxonList is contained in this Taxa.
     * @param taxa a collection of taxa
     * @return true if any of taxa is in this Taxa
     */
    public boolean containsAny(Collection<Taxon> taxa) {

        for (Taxon taxon : taxa) {
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

    /**
     * Returns true if taxonList is a subset of the taxa in this Taxa.
     * @param taxa a collection of taxa
     * @return true if all of taxa is in this Taxa
     */
    public boolean containsAll(Collection<Taxon> taxa) {

        for (Taxon taxon : taxa) {
            if (!taxa.contains(taxon)) {
                return false;
            }
        }
        return true;
    }

	public String getId() { return id; }
	public void setId(String id) { this.id = id; }

	public int compareTo(Taxa o) {
		return getId().compareTo(o.getId());
	}

	public String toString() { return id; }

    public Iterator<Taxon> iterator() {
        return taxa.iterator();
    }

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

    @Override
    public boolean equals(final Object o) {
        if (!(o instanceof TaxonList)) return false;
        return Utils.areTaxaIdentical(this, (TaxonList)o);
    }
}

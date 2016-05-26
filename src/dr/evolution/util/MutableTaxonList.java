/*
 * MutableTaxonList.java
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

/**
 * Interface for a mutable list of taxa.
 *
 * @version $Id: MutableTaxonList.java,v 1.6 2005/05/24 20:25:57 rambaut Exp $
 *
 * @author Alexei Drummond
 */
public interface MutableTaxonList extends TaxonList {

	/**
	 * @return the index of the newly added taxon.
	 * Add the given taxon to the end of the taxonList;
	 */
	int addTaxon(Taxon taxon);
	
	/**
	 * @return true if the taxon was removed, false if it wasn't in the taxonList
	 * Removes the given taxon and renumbers the taxa above it in the list.
	 */
	boolean removeTaxon(Taxon taxon);

	/**
	 * Sets the ID of the ith taxon.
	 */
	public void setTaxonId(int taxonIndex, String id);

	/**
	 * Sets an named attribute for a given taxon.
	 * @param taxonIndex the index of the taxon whose attribute is being set.
	 * @param name the name of the attribute.
	 * @param value the new value of the attribute.
	 */
	public void setTaxonAttribute(int taxonIndex, String name, Object value);

	/**
	 * Adds a listener to this taxon list.
	 */
	void addMutableTaxonListListener(MutableTaxonListListener listener);
}

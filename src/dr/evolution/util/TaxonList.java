/*
 * TaxonList.java
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

import java.util.*;

/**
 * Interface for a list of taxa.
 *
 * @version $Id: TaxonList.java,v 1.16 2006/09/05 13:29:34 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public interface TaxonList extends Identifiable {

	/**
	 * @return a count of the number of taxa in the list.
	 */
	public int getTaxonCount();

	/**
	 * @return the ith taxon.
	 */
	public Taxon getTaxon(int taxonIndex);

	/**
	 * @return the ID of the ith taxon.
	 */
	public String getTaxonId(int taxonIndex);

	/**
	 * returns the index of the taxon with the given id.
	 */
	int getTaxonIndex(String id);

	/**
	 * returns the index of the given taxon.
	 */
	int getTaxonIndex(Taxon taxon);

	/**
	 * @return an object representing the named attributed for the given taxon.
	 * @param taxonIndex the index of the taxon whose attribute is being fetched.
	 * @param name the name of the attribute of interest.
	 */
	public Object getTaxonAttribute(int taxonIndex, String name);

	class Utils {

		public static boolean hasAttribute(TaxonList taxa, int index, String name) {
			return taxa.getTaxonAttribute(index, name) != null;
		}

		public static Set getTaxonListIdSet(TaxonList taxa) {
			Set taxaSet = new HashSet();
			for (int i =0; i < taxa.getTaxonCount(); i++) {
				taxaSet.add(taxa.getTaxonId(i));
			}
			return taxaSet;
		}
	}

	class MissingTaxonException extends Exception {
		public MissingTaxonException(String message) { super(message); }
	}
}

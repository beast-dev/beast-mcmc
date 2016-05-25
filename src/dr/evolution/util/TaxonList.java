/*
 * TaxonList.java
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

import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Interface for a list of taxa.
 *
 * @version $Id: TaxonList.java,v 1.16 2006/09/05 13:29:34 rambaut Exp $
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public interface TaxonList extends Identifiable, Iterable<Taxon> {

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
     * returns the taxa as a Java list
     * @return
     */
    List<Taxon> asList();

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

		public static Set<String> getTaxonListIdSet(TaxonList taxa) {
			Set<String> taxaSet = new HashSet<String>();
			for (int i =0; i < taxa.getTaxonCount(); i++) {
				taxaSet.add(taxa.getTaxonId(i));
			}
			return taxaSet;
		}

        public static int findDuplicateTaxon(TaxonList taxonList) {
            Set<String> taxaSet = new HashSet<String>();
                        for (int i = 0; i < taxonList.getTaxonCount(); i++) {
                Taxon taxon = taxonList.getTaxon(i);
                if (taxaSet.contains(taxon.getId())) {
                    return i;
                }
                taxaSet.add(taxon.getId());
            }
            return -1;
        }

        public static boolean areTaxaIdentical(TaxonList taxa1, TaxonList taxa2) {
            if (taxa1.getTaxonCount() != taxa2.getTaxonCount()) {
                return false;
            }
            for (int i =0; i < taxa1.getTaxonCount(); i++) {
                if (taxa2.getTaxonIndex(taxa1.getTaxon(i)) == -1) {
                    return false;
                }
            }
            return true;
        }


	}

	class MissingTaxonException extends Exception {
		/**
		 *
		 */
		private static final long serialVersionUID = 1864895946392309485L;

		public MissingTaxonException(String message) { super(message); }
	}
}

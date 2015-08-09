/*
 * LinkedGroup.java
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

package dr.evolution;

import dr.evolution.util.TaxonList;
import dr.util.Identifiable;

/**
 * @author Aaron Darling (koadman)
 */
public class LinkedGroup implements Identifiable {

	protected String id = null;
	TaxonList taxa;
	double probability;	// probability that reads in this group are linked to each other

	public LinkedGroup(TaxonList taxa, double probability){
		this.taxa = taxa;
		this.probability = probability;		
	}

	public LinkedGroup(TaxonList taxa){
		// default constructor assumes probability of linkage is 1
		this(taxa, 1.0);
	}
	
	public double getLinkageProbability(){
		return probability;
	}
	public TaxonList getLinkedReads(){
		return taxa;
	}


	public String getId() {
		return id;
	}

	
	public void setId(String id) {
		this.id = id;
	}

}

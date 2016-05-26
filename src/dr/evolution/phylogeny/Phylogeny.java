/*
 * Phylogeny.java
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

package dr.evolution.phylogeny;

public abstract class Phylogeny {
	
	public static final int INCREMENT=10000;
	
	Lineage[] lineages = new Lineage[INCREMENT];
	int lineageCount = 0;
	int availableSize = 0;
	
	int[] extantLineages;
	
	Phylogeny() {
	}

    public void lineageBirth(Lineage lineage, double birthTime) {

        int descendantIndex1 = createLineage(lineage, birthTime);
        int descendantIndex2 = createLineage(lineage, birthTime);

        int lineageIndex = extantLineages[lineage.getExtantIndex()];
        extantLineages[lineage.getExtantIndex()] = descendantIndex1;
    }

    public void lineageDeath(Lineage lineage, double deathTime) {

        lineage.setDeathTime(deathTime);
    }

	private int createLineage(Lineage parent, double birthTime) {
		int index = lineageCount;
		
		if (availableSize > 0) {
			// there are some unused lineage available so return one..
			availableSize--;
		} else {
		
			if (lineageCount == lineages.length) {
				// run out of space in the array so reallocate it...
				Lineage[] newLineages = new Lineage[lineages.length + INCREMENT];
				for (int i = 0; i < lineages.length; i++) {
					newLineages[i] = lineages[i];
				}
				
				lineages = newLineages;
			}
			lineages[index] = newLineage(parent, birthTime);
			lineages[index].setIndex(index);
		}
		
		lineageCount++;
			
		return index;
	}
	
	/**
	 * Override this to return a new class that implements the Lineage interface
	 */
	protected abstract Lineage newLineage(Lineage parent, double birthTime);
	
}
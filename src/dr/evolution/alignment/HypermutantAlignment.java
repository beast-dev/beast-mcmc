/*
 * HypermutantAlignment.java
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

package dr.evolution.alignment;

import dr.evolution.datatype.*;
import dr.evolution.sequence.Sequence;
import dr.evolution.util.Taxon;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * An alignment class that takes another alignment and converts it on the fly
 * to a different dataType.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 *
 * @version $Id: ConvertAlignment.java,v 1.29 2005/05/24 20:25:55 rambaut Exp $
 */
public class HypermutantAlignment extends WrappedAlignment
{
    public enum APOBECType {
        ALL("all"),
        BOTH("both"),
        HA3G("hA3G"),
        HA3F("hA3F");


        APOBECType(String label) {
            this.label = label;
        }

        public String toString() {
            return label;
        }

        final String label;
    }


	/**
	 * Constructor.
	 */
	public HypermutantAlignment(APOBECType type, Alignment alignment) {
        super(alignment);
        this.type = type;

        if (alignment.getDataType().getType() != DataType.NUCLEOTIDES) {
            throw new RuntimeException("HypermutantAlignment can only convert nucleotide alignments");
        }

        mutatedContextCounts = new int[getTaxonCount()];
        unmutatedContextCounts = new int[getTaxonCount()];

        countContexts();

	}

	/**
	 * @return the sequence state at (taxon, site)
	 */
	public int getState(int taxonIndex, int siteIndex) {
		int state = alignment.getState(taxonIndex, siteIndex);
        if (state == Nucleotides.A_STATE) {
            int nextState = getNextContextState(taxonIndex, siteIndex);

            if (    (type == APOBECType.ALL) || // consider all As as G->As 
                    (type == APOBECType.HA3G && nextState == Nucleotides.G_STATE) ||
                    (type == APOBECType.HA3F && nextState == Nucleotides.A_STATE) ||
                    (type == APOBECType.BOTH && (nextState == Nucleotides.G_STATE || nextState == Nucleotides.A_STATE))
                    ) {
                state = Nucleotides.R_STATE;
            }
        }

		return state;
	}

    /**
     * This finds the next state downstream of siteIndex, skipping over gaps
     */
    private int getNextContextState(int taxonIndex, int siteIndex) {
        int nextState = Nucleotides.GAP_STATE;

        int i = siteIndex + 1;
        while (nextState == Nucleotides.GAP_STATE && i < getSiteCount()) {
            nextState = alignment.getState(taxonIndex, i);
            i++;
        }

        return nextState;
    }

    private void countContexts() {
        for (int i = 0; i < getTaxonCount(); i++) {
            for (int j = 0; j < getSiteCount(); j++) {
                int state = alignment.getState(i, j);
                if (state == Nucleotides.A_STATE || state == Nucleotides.G_STATE) {
                    int nextState = getNextContextState(i, j);

                    if (    (type == APOBECType.ALL) || // consider all As as G->As
                            (type == APOBECType.HA3G && nextState == Nucleotides.G_STATE) ||
                            (type == APOBECType.HA3F && nextState == Nucleotides.A_STATE) ||
                            (type == APOBECType.BOTH && (nextState == Nucleotides.G_STATE || nextState == Nucleotides.A_STATE))
                            ) {
                        if (state == Nucleotides.A_STATE) {
                            // a mutated context
                            mutatedContextCounts[i] ++;
                        } else {
                            // an unmutated context
                            unmutatedContextCounts[i] ++;
                        }
                    }
                }
            }
        }
    }

    public int[] getMutatedContextCounts() {
        return mutatedContextCounts;
    }

    public int[] getUnmutatedContextCounts() {
        return unmutatedContextCounts;
    }

    public int getMutatedContextCount() {
        int total = 0;
        for (int count: mutatedContextCounts) {
            total += count;
        }
        return total;
    }

    public int getUnmutatedContextCount() {
        int total = 0;
        for (int count: unmutatedContextCounts) {
            total += count;
        }
        return total;
    }

    // **************************************************************
	// INSTANCE VARIABLES
	// **************************************************************

	private APOBECType type = null;
    private int[] mutatedContextCounts;
    private int[] unmutatedContextCounts;
}
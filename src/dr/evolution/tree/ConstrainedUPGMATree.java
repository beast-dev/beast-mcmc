/*
 * ConstrainedUPGMATree.java
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

package dr.evolution.tree;

import dr.evolution.distance.DistanceMatrix;
import dr.evolution.util.TaxonList;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

/**
 * UPGMA that keeps given taxon sets monophyletic: two clusters are joined only if the union is inside, contains or
 * misses each set, and among those joins the closest pair is taken. Every cluster then stays compatible with every
 * set, and for nested or disjoint sets some join is always allowed, so the tree has each set as a clade.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public class ConstrainedUPGMATree extends UPGMATree {

    public ConstrainedUPGMATree(DistanceMatrix distanceMatrix, List<TaxonList> monophyletic) {
        super(withConstraints(distanceMatrix, monophyletic));
    }

    // the constraints are needed during the superclass constructor, which runs the clustering
    private static final ThreadLocal<List<BitSet>> PENDING = new ThreadLocal<List<BitSet>>();

    private static DistanceMatrix withConstraints(DistanceMatrix distanceMatrix, List<TaxonList> monophyletic) {
        List<BitSet> sets = new ArrayList<BitSet>();
        for (TaxonList taxa : monophyletic) {
            BitSet set = new BitSet();
            for (int i = 0; i < taxa.getTaxonCount(); i++) {
                int index = distanceMatrix.getTaxonIndex(taxa.getTaxon(i).getId());
                if (index < 0) {
                    throw new IllegalArgumentException("Taxon " + taxa.getTaxon(i).getId() + " of a monophyletic " +
                            "set is not in the distance matrix");
                }
                set.set(index);
            }
            sets.add(set);
        }
        for (int a = 0; a < sets.size(); a++) {
            for (int b = a + 1; b < sets.size(); b++) {
                if (!compatible(sets.get(a), sets.get(b))) {
                    throw new IllegalArgumentException("Monophyletic sets " + monophyletic.get(a).getId() + " and " +
                            monophyletic.get(b).getId() + " overlap without one containing the other");
                }
            }
        }
        PENDING.set(sets);
        return distanceMatrix;
    }

    private List<BitSet> constraints;
    private BitSet[] members; // the taxa of the cluster at each storage index (alias)

    private static boolean compatible(BitSet a, BitSet b) {
        BitSet intersection = (BitSet) a.clone();
        intersection.and(b);
        return intersection.isEmpty() || intersection.equals(a) || intersection.equals(b);
    }

    private boolean allowed(int i, int j) {
        BitSet union = (BitSet) members[alias[i]].clone();
        union.or(members[alias[j]]);
        for (BitSet constraint : constraints) {
            if (!compatible(union, constraint)) {
                return false;
            }
        }
        return true;
    }

    @Override
    protected void init(DistanceMatrix distanceMatrix) {
        super.init(distanceMatrix);
        constraints = PENDING.get();
        PENDING.remove();
        members = new BitSet[numClusters];
        for (int i = 0; i < numClusters; i++) {
            members[i] = new BitSet();
            members[i].set(i);
        }
    }

    @Override
    protected void findNextPair() {
        if (numClusters < 2) {
            return; // the clustering loop asks once more after the last join
        }
        double dmin = Double.POSITIVE_INFINITY;
        besti = -1;
        for (int i = 0; i < numClusters - 1; i++) {
            for (int j = i + 1; j < numClusters; j++) {
                if ((besti < 0 || getDist(i, j) < dmin) && allowed(i, j)) {
                    dmin = getDist(i, j);
                    besti = i;
                    bestj = j;
                }
            }
        }
        if (besti < 0) {
            throw new IllegalStateException("No join keeps the monophyletic sets together");
        }
        abi = alias[besti];
        abj = alias[bestj];
    }

    @Override
    protected void newCluster() {
        members[abi].or(members[abj]);
        members[abj] = new BitSet();
        super.newCluster();
    }
}

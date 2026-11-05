/*
 * CandidateSet.java
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

package dr.evomodel.tree.search;

import dr.evolution.tree.Tree;
import dr.math.MathUtils;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * The best distinct topologies the search has reached (plan 7 §7.4, after IQ-TREE's CandidateSet): up to a fixed
 * number, sorted by F; a topology already present keeps its higher F; when full, the worst is dropped. Perturbation
 * draws its parents uniformly from the best few.
 *
 * @author Xiang Ji
 * @author Marc A. Suchard
 */
public final class CandidateSet {

    public static final class Candidate {
        public final long hash;
        public final double f;
        public final Tree tree; // a copy, for the output
        final TreeStateSnapshot snapshot;
        public final int climb; // the hill-climb that reached it

        Candidate(long hash, double f, Tree tree, TreeStateSnapshot snapshot, int climb) {
            this.hash = hash;
            this.f = f;
            this.tree = tree;
            this.snapshot = snapshot;
            this.climb = climb;
        }

        public TreeStateSnapshot getSnapshot() {
            return snapshot;
        }
    }

    private final int capacity;
    private final List<Candidate> candidates = new ArrayList<Candidate>();

    public CandidateSet(int capacity) {
        this.capacity = Math.max(1, capacity);
    }

    /**
     * @return whether the set changed
     */
    public boolean add(long hash, double f, Tree tree, TreeStateSnapshot snapshot, int climb) {
        for (int k = 0; k < candidates.size(); ++k) {
            if (candidates.get(k).hash == hash) {
                if (!(f > candidates.get(k).f)) {
                    return false;
                }
                candidates.remove(k);
                break;
            }
        }
        int at = 0;
        while (at < candidates.size() && candidates.get(at).f >= f) {
            at++;
        }
        if (at >= capacity) {
            return false;
        }
        candidates.add(at, new Candidate(hash, f, tree, snapshot, climb));
        if (candidates.size() > capacity) {
            candidates.remove(candidates.size() - 1);
        }
        return true;
    }

    /**
     * @return a candidate drawn uniformly from the best pool (MathUtils)
     */
    public Candidate drawParent(int pool) {
        if (candidates.isEmpty()) {
            throw new IllegalStateException("The candidate set is empty");
        }
        return candidates.get(MathUtils.nextInt(Math.min(Math.max(1, pool), candidates.size())));
    }

    public List<Candidate> getCandidates() {
        return Collections.unmodifiableList(candidates);
    }

    public int size() {
        return candidates.size();
    }

    public void clear() {
        candidates.clear();
    }
}

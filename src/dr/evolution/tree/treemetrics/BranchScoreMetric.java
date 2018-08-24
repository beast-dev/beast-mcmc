/*
 * BranchScoreMetric.java
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

/**
 *
 */
package dr.evolution.tree.treemetrics;

import java.util.*;

import dr.evolution.tree.Clade;
import dr.evolution.tree.Tree;

import static dr.evolution.tree.treemetrics.TreeMetric.Utils.checkTreeTaxa;

/**
 * @author Andrew Rambaut
 * @author Sebastian Hoehna
 *
 */
public class BranchScoreMetric implements TreeMetric {
    public static Type TYPE = Type.BRANCH_SCORE;

    private Tree focalTree;
    private List<Clade> focalClades;
    private final boolean fixedFocalTree;

    public BranchScoreMetric() {
        this.fixedFocalTree = false;
    }

    public BranchScoreMetric(Tree focalTree) {
        this.focalTree = focalTree;
        this.fixedFocalTree = true;
        focalClades = Clade.getCladeList(focalTree);
    }


    @Override
    public double getMetric(Tree tree1, Tree tree2) {

        checkTreeTaxa(tree1, tree2);

        if (tree1 != focalTree) {
            if (fixedFocalTree) {
                // If we set a focal tree in the constructor then it makes sense to check it is the same
                // as the one set here.
                throw new RuntimeException("Focal tree is different from that set in the constructor.");
            }

            // cache tree1 and the pre-computed path for future calls
            focalTree = tree1;
            focalClades = Clade.getCladeList(focalTree);
        }

        List<Clade> clades2 = Clade.getCladeList(tree2);

        return getDistance(focalClades, clades2);
    }

    protected double getDistance(List<Clade> clades1, List<Clade> clades2) {


        Collections.sort(clades1);
        Collections.sort(clades2);
        double distance = 0.0;
        int indexClade2 = 0;
        Clade clade2 = null;
        Clade parent1, parent2 = null;
        double height1, height2;

        for (Clade clade1 : clades1) {

            parent1 = findParent(clade1, clades1);
            height1 = parent1.getHeight() - clade1.getHeight();

            if (indexClade2 < clades2.size()) {
                clade2 = clades2.get(indexClade2);
                parent2 = findParent(clade2, clades2);
            }
            while (clade1.compareTo(clade2) > 0 && indexClade2 < clades2.size()) {
                height2 = parent2.getHeight() - clade2.getHeight();
                distance += height2 * height2;
                indexClade2++;
                if (indexClade2 < clades2.size()) {
                    clade2 = clades2.get(indexClade2);
                    parent2 = findParent(clade2, clades2);
                }
            }
            if (clade1.compareTo(clade2) == 0) {
                height2 = parent2.getHeight() - clade2.getHeight();
                distance += (height1 - height2) * (height1 - height2);
                indexClade2++;
            } else {
                distance += height1 * height1;
            }
        }

        return Math.sqrt(distance);
    }

    private Clade findParent(Clade clade1, List<Clade> clades) {
        Clade parent = null;
        for (Clade clade2 : clades) {
            if (isParent(clade2, clade1)) {
                if (parent == null || parent.getSize() > clade2.getSize())
                    parent = clade2;
            }
        }

        if (parent == null){
            //the case that this clade is the whole tree
            return clade1;
        }

        return parent;
    }

    private boolean isParent(Clade parent, Clade child) {
        if (parent.getSize() <= child.getSize()) {
            return false;
        }

        tmpBits.clear();
        tmpBits.or(parent.getBits());
        tmpBits.xor(child.getBits());

        return tmpBits.cardinality() < parent.getSize();
    }

    @Override
    public Type getType() {
        return TYPE;
    }

    @Override
    public String toString() {
        return getType().getShortName();
    }

    BitSet tmpBits = new BitSet();
}

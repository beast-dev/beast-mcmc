/*
 * TopologyTracer.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evolution.tree.treemetrics;

import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;

import java.io.IOException;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

import static dr.evolution.tree.treemetrics.TreeMetric.Utils.checkTreeTaxa;

/**
 * @author Guy Baele
 * @author Andrew Rambaut
 * Path difference metric according to Steel & Penny (1993)
 */
public class SteelPennyPathDifferenceMetric implements TreeMetric {

    public static Type TYPE = Type.STEEL_PENNY;

    private Tree focalTree;
    private int dim;
    private double[] focalPath;
    private final boolean fixedFocalTree;

    public SteelPennyPathDifferenceMetric() {
        this.fixedFocalTree = false;
    }

    public SteelPennyPathDifferenceMetric(Tree focalTree) {
        this.focalTree = focalTree;
        this.fixedFocalTree = true;
        this.dim = focalTree.getExternalNodeCount() * focalTree.getExternalNodeCount();
        this.focalPath = new double[dim];

        traverse(this.focalTree, this.focalTree.getRoot(), focalPath);
    }

    /**
     * Compute the metric between two trees. If tree1 is not the focal tree provided to the
     * constructor then it will store this as the new focal tree. If the focal tree is constant
     * from call to call then a cached set of precomputation will be used, increasing efficiency.
     * @param tree1
     * @param tree2
     * @return
     */
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
            if (focalPath == null) {
                this.dim = focalTree.getExternalNodeCount() * focalTree.getExternalNodeCount();
                this.focalPath = new double[dim];
            }
            traverse(focalTree, focalTree.getRoot(), focalPath);
        }

        double[] pathTwo = new double[dim];
        traverse(tree2, tree2.getRoot(), pathTwo);

        int n = tree1.getExternalNodeCount();

        double metric = 0.0;
        for (int i = 0; i < n; i++) {
            for (int j = i + 1; j < n; j++) {
                int index = (i * n) + j;
                metric += Math.pow(focalPath[index] - pathTwo[index], 2);
            }
        }
        metric = Math.sqrt(metric);

        return metric;
    }

    private Set<NodeRef> traverse(Tree tree, NodeRef node, double[] lengths) {
        NodeRef left = tree.getChild(node, 0);
        NodeRef right = tree.getChild(node, 1);


        Set<NodeRef> leftSet = null;
        Set<NodeRef> rightSet = null;

        if (!tree.isExternal(left)) {
            leftSet = traverse(tree, left, lengths);
        } else {
            leftSet = Collections.singleton(left);
        }
        if (!tree.isExternal(right)) {
            rightSet = traverse(tree, right, lengths);
        } else {
            rightSet = Collections.singleton(right);
        }

        for (NodeRef tip1 : leftSet) {
            for (NodeRef tip2 : rightSet) {
                int index;
                if (tip1.getNumber() < tip2.getNumber()) {
                    index = (tip1.getNumber() * tree.getExternalNodeCount()) + tip2.getNumber();
                } else {
                    index = (tip2.getNumber() * tree.getExternalNodeCount()) + tip1.getNumber();
                }
                lengths[index] = tree.getNodeHeight(node) * 2
                        - tree.getNodeHeight(tip1)
                        - tree.getNodeHeight(tip2);

            }
        }

        Set<NodeRef> tips = new HashSet<NodeRef>();
        tips.addAll(leftSet);
        tips.addAll(rightSet);

        return tips;
    }

    @Override
    public Type getType() {
        return TYPE;
    }

    @Override
    public String toString() {
        return getType().getShortName();
    }

}

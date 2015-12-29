/*
 * CoalescentSimulator.java
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

package dr.evomodel.coalescent;

import dr.evolution.tree.*;
import dr.evolution.util.Date;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evolution.util.TimeScale;
import dr.inference.distribution.ParametricDistributionModel;
import dr.math.UnivariateFunction;

/**
 * Simulates a set of coalescent intervals given a demographic model.
 *
 * @author Alexei Drummond
 * @version $Id: CoalescentSimulator.java,v 1.43 2005/10/27 10:40:48 rambaut Exp $
 */
public class CoalescentSimulator {

    dr.evolution.coalescent.CoalescentSimulator simulator = new dr.evolution.coalescent.CoalescentSimulator();
    /**
     * Simulates a coalescent tree from a set of subtrees.
     */
    public CoalescentSimulator() {
    }

    /**
     * Simulates a coalescent tree from a set of subtrees.
     *
     * @param subtrees                an array of tree to be used as subtrees
     * @param model                   the demographic model to use
     * @param rootHeight              an optional root height with which to scale the whole tree
     * @param preserveSubtreesHeights true if heights of subtrees should be preserved
     * @return a simulated coalescent tree
     */
    public SimpleTree simulateTree(Tree[] subtrees, DemographicModel model, double rootHeight,
                                   boolean preserveSubtreesHeights) {

        SimpleNode[] roots = new SimpleNode[subtrees.length];
        SimpleTree tree;

        for (int i = 0; i < roots.length; i++) {
            roots[i] = new SimpleNode(subtrees[i], subtrees[i].getRoot());
        }

        // if just one taxonList then finished
        if (roots.length == 1) {
            tree = new SimpleTree(roots[0]);
        } else {
            tree = new SimpleTree(simulator.simulateCoalescent(roots, model.getDemographicFunction()));
        }

        if (!Double.isNaN(rootHeight) && rootHeight > 0.0) {
            if (preserveSubtreesHeights) {
                limitNodes(tree, rootHeight - 1e-12);
                tree.setRootHeight(rootHeight); 
            } else {
                attemptToScaleTree(tree, rootHeight);
            }
        }

        return tree;
    }


    /**
     * Simulates a coalescent tree, given a taxon list.
     *
     * @param taxa  the set of taxa to simulate a coalescent tree between
     * @param model the demographic model to use
     * @return a simulated coalescent tree
     */
    public SimpleTree simulateTree(TaxonList taxa, DemographicModel model) {
        return simulator.simulateTree(taxa, model.getDemographicFunction());
    }

    public void attemptToScaleTree(MutableTree tree, double rootHeight) {
        // avoid empty tree
        if (tree.getRoot() == null) return;

        double scale = rootHeight / tree.getNodeHeight(tree.getRoot());
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            NodeRef n = tree.getInternalNode(i);
            tree.setNodeHeight(n, tree.getNodeHeight(n) * scale);
        }
        MutableTree.Utils.correctHeightsForTips(tree);
    }

    public int sizeOfIntersection(TaxonList tl1, TaxonList tl2) {
        int nIn = 0;
        for (int j = 0; j < tl1.getTaxonCount(); ++j) {
            if (tl2.getTaxonIndex(tl1.getTaxon(j)) >= 0) {
                ++nIn;
            }
        }
        return nIn;
    }

    public boolean contained(TaxonList taxons, TaxonList taxons1) {
        return sizeOfIntersection(taxons, taxons1) == taxons.getTaxonCount();
    }

    /**
     * Clip nodes height above limit.
     *
     * @param tree  to clip
     * @param limit height limit
     */
    private void limitNodes(MutableTree tree, double limit) {
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            final NodeRef n = tree.getInternalNode(i);
            if (tree.getNodeHeight(n) > limit) {
                tree.setNodeHeight(n, limit);
            }
        }
        MutableTree.Utils.correctHeightsForTips(tree);
    }

    public static class TaxaConstraint {
        public final TaxonList taxons;
        public final double lower;
        public final boolean isMonophyletic;
        public double upper;

        public TaxaConstraint(TaxonList taxons, ParametricDistributionModel p, boolean isMono) {
            this.taxons = taxons;
            this.isMonophyletic = isMono;

            if (p != null) {
                final UnivariateFunction univariateFunction = p.getProbabilityDensityFunction();
                lower = univariateFunction.getLowerBound();
                upper = univariateFunction.getUpperBound();
            } else {
                lower = 0;
                upper = Double.POSITIVE_INFINITY;
            }
        }

        public TaxaConstraint(TaxonList taxons, double low, double high, boolean isMono) {
            this.taxons = taxons;
            this.isMonophyletic = isMono;
            upper = high;
            lower = low;
        }

        public boolean realLimits() {
            return lower != 0 || upper != Double.POSITIVE_INFINITY;
        }
    }

}

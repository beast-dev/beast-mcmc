/*
 * SubtreeRateStatistic.java
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

package dr.evomodel.tree;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodelxml.tree.RateStatisticParser;
import dr.inference.model.Statistic;
import dr.math.matrixAlgebra.Vector;
import dr.stats.DiscreteStatistics;

import java.util.*;

/**
 * A statistic that tracks the mean, variance and coefficent of variation of the rates in a subtree (or the complement thereof).
 *
 * @author Andy Magee
 */
public class SubtreeRateStatistic extends TreeStatistic {

    public SubtreeRateStatistic(String name,
                                Tree tree,
                                BranchRateModel branchRateModel,
                                String mode,
                                TaxonList mrcaTaxa,
                                boolean complement,
                                boolean includeStem,
                                boolean returnGeneticDistance) throws TreeUtils.MissingTaxonException {
        super(name);
        this.tree = tree;
        this.branchRateModel = branchRateModel;
        this.mode = mode;
        this.complement = complement;
        this.includeStem = includeStem;
        this.returnGeneticDistance = returnGeneticDistance;

        this.taxa = new HashSet<>();
        for (Taxon taxon : mrcaTaxa) {
            this.taxa.add(taxon);
        }

        this.leafSet = TreeUtils.getLeavesForTaxa(tree, mrcaTaxa);
    }

    public void setTree(Tree tree) {
        this.tree = tree;
    }

    public Tree getTree() {
        return tree;
    }

    public int getDimension() {
        return 1;
    }

    private void addBranchesToClade(boolean[] use, NodeRef node) {
        for (int i = 0; i < tree.getChildCount(node); i++) {
            NodeRef child = tree.getChild(node, i);
            use[child.getNumber()] = true;
            if ( !tree.isExternal(node) ) {
                addBranchesToClade(use, child);
            }
        }
    }

    private boolean[] getBranchesToUse() {
        boolean[] use = new boolean[tree.getNodeCount()];
        NodeRef node = TreeUtils.getCommonAncestorNode(tree, leafSet);

        if ( includeStem ) {
            use[node.getNumber()] = true;
        }
        addBranchesToClade(use, node);

        if (complement) {
            for (int i = 0; i < use.length; i++) {
                use[i] = !use[i];
            }
            use[tree.getRoot().getNumber()] = false;
        }

        return use;
    }

    // Recursively add rates and branch lengths
    private void getRates(List<Double> rates, List<Double> branchLengths) {
        boolean[] use = getBranchesToUse();
        for (int i = 0; i < use.length; i++) {
            if (use[i]) {
                NodeRef node = tree.getNode(i);
                rates.add(branchRateModel.getBranchRate(tree,node));
                branchLengths.add(tree.getBranchLength(node));
            }
        }
    }

    /**
     * @return the statistic
     */
    public double getStatisticValue(int dim) {

        List<Double> rateList = new ArrayList<Double>();
        List<Double> branchLengthList = new ArrayList<Double>();

        // Get rates and branch lengths for just this subtree
        getRates(rateList, branchLengthList);

        double[] rates = new double[rateList.size()];
        double[] branchLengths = new double[rateList.size()];
        for ( int i = 0; i < rates.length; i++) {
            rates[i] = rateList.get(i);
            branchLengths[i] = branchLengthList.get(i);
        }

        if ( returnGeneticDistance ) {
            double distance = 0;
            for (int i = 0; i < rates.length; i++) {
                distance += rates[i] * branchLengths[i];
            }
            return distance;
        // as RateStatistic.java
        } else if (mode.equals(RateStatisticParser.MEAN)) {
            double totalWeightedRate = 0.0;
            double totalTreeLength = 0.0;
            for (int i = 0; i < rates.length; i++) {
                totalWeightedRate += rates[i] * branchLengths[i];
                totalTreeLength += branchLengths[i];
            }
            return totalWeightedRate / totalTreeLength;
        } else if (mode.equals(RateStatisticParser.VARIANCE)) {
            return DiscreteStatistics.variance(rates);
        } else if (mode.equals(RateStatisticParser.COEFFICIENT_OF_VARIATION)) {
            // don't compute mean twice
            final double mean = DiscreteStatistics.mean(rates);
            return Math.sqrt(DiscreteStatistics.variance(rates, mean)) / mean;
        }

        throw new IllegalArgumentException();
    }

    private Tree tree;
    private BranchRateModel branchRateModel;
    private boolean complement;
    private final boolean includeStem;
    private final boolean returnGeneticDistance;
    private String mode;
    private final Set<Taxon>  taxa;
    private final Set<String> leafSet;
}

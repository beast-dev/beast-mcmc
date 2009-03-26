/*
 * DiscretizedBranchRates.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.DiscretizedBranchRatesParser;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @author Michael Defoin Platel
 * @version $Id: DiscretizedBranchRates.java,v 1.11 2006/01/09 17:44:30 rambaut Exp $
 */
public class DiscretizedBranchRates extends AbstractModel implements BranchRateModel {

    private final ParametricDistributionModel distributionModel;
    private final TreeModel tree;

    // The rate categories of each branch
    private final Parameter rateCategoryParameter;

    // the index of the root node.
    private int rootNodeNumber;
    private int storedRootNodeNumber;

    private final int categoryCount;
    private final double step;
    private final double[] rates;

    // my attempt to normalized rates stymied by store/restore issues - AJD
    //private boolean normalize = false;
    //private double normalizedRate = 1.0;
    //private boolean meanIsCalculated = false;
    //private double meanRate = 1.0;

    //overSampling control the number of effective categories

    public DiscretizedBranchRates(TreeModel tree, Parameter rateCategoryParameter, ParametricDistributionModel model,
                                  int overSampling) {

        super(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES);
        this.tree = tree;

        // todo we should set this automatically if this is a requirment
        int nBranches = tree.getNodeCount() - 1;
        if (rateCategoryParameter.getDimension() != nBranches) {
            throw new IllegalArgumentException("The rate category parameter must be of length " + nBranches + " (nodes in the tree - 1)");
        }

        categoryCount = nBranches * overSampling;
        step = 1.0 / (double) categoryCount;

        rates = new double[categoryCount];

        this.distributionModel = model;

        //Force the boundaries of rateCategoryParameter to match the category count
        Parameter.DefaultBounds bound = new Parameter.DefaultBounds(categoryCount - 1, 0, rateCategoryParameter.getDimension());
        rateCategoryParameter.addBounds(bound);

        this.rateCategoryParameter = rateCategoryParameter;

        for (int i = 0; i < rateCategoryParameter.getDimension(); i++) {
            int index = (int) Math.floor((i + 0.5) * overSampling);
            rateCategoryParameter.setParameterValue(i, index);
        }

        addModel(model);
        addModel(tree);

        addParameter(rateCategoryParameter);

        rootNodeNumber = tree.getRoot().getNumber();
        storedRootNodeNumber = rootNodeNumber;

        setupRates();
        shuffleIndices();
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == distributionModel) {
            setupRates();
            fireModelChanged();
        } else if (model == tree) {
            shuffleIndices();
            //meanIsCalculated = false;
        }
    }

    protected final void handleParameterChangedEvent(Parameter parameter, int index, Parameter.ChangeType type) {
        fireModelChanged(null, getNodeNumberFromCategoryIndex(index));
    }

    protected void storeState() {
        storedRootNodeNumber = rootNodeNumber;
    }

    protected void restoreState() {
        setupRates();
        rootNodeNumber = storedRootNodeNumber;
    }

    protected void acceptState() {
    }

    public double getBranchRate(final Tree tree, final NodeRef node) {

        assert !tree.isRoot(node) : "root node doesn't have a rate!";

//       final double scale = 1.0;
//        if (normalize) {
//            calculateMeanRate(tree);
//            scale = normalizedRate/meanRate;
//        }

        final int nodeNumber = node.getNumber();

        assert nodeNumber != rootNodeNumber :
                "INTERNAL ERROR! node with number " + rootNodeNumber + " should be the root node.";

        final int dim = getCategoryIndexFromNodeNumber(nodeNumber);
        final int rateCategory = (int) Math.round(rateCategoryParameter.getParameterValue(dim));


        if (rateCategory < rateCategoryParameter.getBounds().getLowerLimit(0) || rateCategory > rateCategoryParameter.getBounds().getUpperLimit(0)) {
            throw new IllegalArgumentException("INTERNAL ERROR! invalid category number " + rateCategory);
        }

//        // todo this should be nodeNumber, no?
//        assert (rateCategory >= rateCategoryParameter.getBounds().getLowerLimit(nodeNumber) &&
//                rateCategory <= rateCategoryParameter.getBounds().getUpperLimit(nodeNumber)) :
//            ("INTERNAL ERROR! invalid category number " + rateCategory);

        return rates[rateCategory] /** scale*/;
    }

    /*   private void calculateMeanRate(Tree tree) {

        if (meanIsCalculated) return;

        int length = tree.getNodeCount();
        int rootIndex = 0;

        double[] branchLengths = new double[length];

        for (int i = 0; i < length; i++) {
            NodeRef child = tree.getNode(i);
            if (!tree.isRoot(child)) {
                NodeRef parent = tree.getParent(child);
                branchLengths[child.getNumber()] = tree.getNodeHeight(parent) - tree.getNodeHeight(child);
            } else {
                rootIndex = i;
            }
        }

        double totalWeightedRate = 0.0;
        double totalTreeLength = 0.0;
        for (int i = 0; i < length; i++) {

            if (i != rootIndex) {
                int nodeNumber = tree.getNode(i).getNumber();
                int rateIndex = (int) Math.round(
                        rateCategoryParameter.getParameterValue(
                                getCategoryIndexFromNodeNumber(nodeNumber)));

                totalWeightedRate += rates[rateIndex] * branchLengths[nodeNumber];
                totalTreeLength += branchLengths[nodeNumber];
            }
        }

        meanRate = totalWeightedRate / totalTreeLength;
        meanIsCalculated = true;

    }*/

/*    public NodeRef getNodeForParameter(Parameter parameter, int index) {

        if (parameter != rateCategoryParameter) {
            throw new RuntimeException("Expecting " + rateCategoryParameter + ", but got " + parameter);
        }

        if (index == -1) {
            throw new RuntimeException("Expecting non-negative index!");
        }

        NodeRef node;
        if (index < rootNodeNumber) {
            node = tree.getNode(index);
            if (node.getNumber() != index) throw new RuntimeException();
        } else {
            node = tree.getNode(index - 1);
            if (node.getNumber() != index - 1) throw new RuntimeException();
        }
        return node;
    }*/

    public String getBranchAttributeLabel() {
        return "rate";
    }

    public String getAttributeForBranch(Tree tree, NodeRef node) {
        return Double.toString(getBranchRate(tree, node));
    }

    public int getNodeNumberFromCategoryIndex(int categoryIndex) {
        if (categoryIndex >= rootNodeNumber) return categoryIndex + 1;
        return categoryIndex;
    }

    public int getCategoryIndexFromNodeNumber(int nodeNumber) {
        if (nodeNumber > rootNodeNumber) return nodeNumber - 1;
        return nodeNumber;
    }

    /**
     * Calculates the actual rates corresponding to the category indices.
     */
    protected void setupRates() {

        double z = step / 2.0;
        for (int i = 0; i < categoryCount; i++) {
            rates[i] = distributionModel.quantile(z);
            z += step;
        }
        //meanIsCalculated = false;
    }

    private void shuffleIndices() {
        final int newRootNodeNumber = tree.getRoot().getNumber();

        //if (newRootNodeNumber != rootNodeNumber) {
        //    System.out.println("old root node number =" + rootNodeNumber);
        //    System.out.println("new root node number =" + newRootNodeNumber);
        //}

        if (rootNodeNumber > newRootNodeNumber) {

            //for (int i = 0; i < rateCategoryParameter.getDimension(); i++) {
            //    System.out.print((int)Math.round(rateCategoryParameter.getParameterValue(i)) + "\t");
            //}
            //System.out.println();

//            int oldRateIndex = (int) Math.round(
            final double oldRateIndex =
                    rateCategoryParameter.getParameterValue(newRootNodeNumber);

            final int end = Math.min(rateCategoryParameter.getDimension() - 1, rootNodeNumber);
            for (int i = newRootNodeNumber; i < end; i++) {
                rateCategoryParameter.setParameterValue(i, rateCategoryParameter.getParameterValue(i + 1));
            }

            rateCategoryParameter.setParameterValue(end, oldRateIndex);

            //for (int i = 0; i < rateCategoryParameter.getDimension(); i++) {
            //    System.out.print((int)Math.round(rateCategoryParameter.getParameterValue(i)) + "\t");
            //}
            //System.out.println();

        } else if (rootNodeNumber < newRootNodeNumber) {

            //System.out.println("old root node number =" + rootNodeNumber);
            //System.out.println("new root node number =" + newRootNodeNumber);

            //for (int i = 0; i < rateCategoryParameter.getDimension(); i++) {
            //    System.out.print((int)Math.round(rateCategoryParameter.getParameterValue(i)) + "\t");
            //}
            //System.out.println();

            final int end = Math.min(rateCategoryParameter.getDimension() - 1, newRootNodeNumber);

            final int oldRateIndex = (int) Math.round(rateCategoryParameter.getParameterValue(end));

            for (int i = end; i > rootNodeNumber; i--) {
                rateCategoryParameter.setParameterValue(i, rateCategoryParameter.getParameterValue(i - 1));
            }

            rateCategoryParameter.setParameterValue(rootNodeNumber, oldRateIndex);

            //for (int i = 0; i < rateCategoryParameter.getDimension(); i++) {
            //    System.out.print((int)Math.round(rateCategoryParameter.getParameterValue(i)) + "\t");
            //}
            //System.out.println();
        }
        rootNodeNumber = newRootNodeNumber;
    }

    /*public void setNormalizedMean(double normalizedMean) {
        this.normalizedRate = normalizedMean;
        normalize = true;
        meanIsCalculated = false;
    }*/
}

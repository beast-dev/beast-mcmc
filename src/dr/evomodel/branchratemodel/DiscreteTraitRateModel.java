/*
 * ColouredTreeRateModel.java
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

import dr.evolution.alignment.PatternList;
import dr.evolution.parsimony.FitchParsimony;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * This Branch Rate Model takes a ancestral state likelihood and
 * gives the rate for each branch of the tree based on the child state (for now).
 *
 * @author Alexei Drummond
 * @author Marc Suchard
 */
public class DiscreteTraitRateModel extends AbstractBranchRateModel {
    enum Mode {
        NODE_STATES,
        DWELL_TIMES,
        PARSIMONY
    }

    public static final String DISCRETE_TRAIT_RATE_MODEL = "discreteTraitRateModel";
    public static final String RATES = "rates";
    public static final String TRAIT_INDEX = "traitIndex";
    public static final String TRAIT_NAME = "traitName";

    private TreeTrait trait = null;
    private Parameter ratesParameter;
    private int traitIndex;
    private boolean normKnown = false;
    private boolean storedNormKnown = false;
    private double norm = 1.0;
    private double storedNorm = 1.0;

    private FitchParsimony fitchParsimony;

    private boolean treeChanged = true;
    private boolean shouldRestoreTree = false;

    private Mode mode;

//    private int treeInitializeCounter = 0;

    /**
     * A constructor for the (crude) parsimony reconstruction form of this class.
     * @param treeModel
     * @param patternList
     * @param traitIndex
     * @param ratesParameter
     */
    public DiscreteTraitRateModel(TreeModel treeModel, PatternList patternList, int traitIndex, Parameter ratesParameter) {

        this(treeModel, traitIndex, ratesParameter);

        if (!TaxonList.Utils.getTaxonListIdSet(treeModel).equals(TaxonList.Utils.getTaxonListIdSet(patternList))) {
            throw new IllegalArgumentException("Tree model and pattern list must have the same list of taxa!");
        }

        fitchParsimony = new FitchParsimony(patternList, false);
        mode = Mode.PARSIMONY;
    }

    /**
     * A constructor for a node-sampled discrete trait
     * @param treeModel
     * @param trait
     * @param traitIndex
     * @param ratesParameter
     */
    public DiscreteTraitRateModel(TreeModel treeModel, TreeTrait trait, int traitIndex, Parameter ratesParameter) {

        this(treeModel, traitIndex, ratesParameter);

//        if (trait.getTreeModel() != treeModel)
//            throw new IllegalArgumentException("Tree Models for ancestral state tree likelihood and target model of these rates must match!");

        this.trait = trait;

        if (trait.getClass().equals(Integer.class)) {
            // Assume the trait is one or more discrete traits reconstructed at nodes
            mode = Mode.NODE_STATES;
            if (traitIndex < 0 || traitIndex >= trait.getDimension()) {
                throw new IllegalArgumentException("The trait index must be within the dimension of the trait.");
            }
        } else if (trait.getClass().equals(double[].class)) {
            // Assume the trait itself is the dwell times for the individual states on the branch above the node
            mode = Mode.DWELL_TIMES;
            if (trait.getDimension() != ratesParameter.getDimension()) {
                throw new IllegalArgumentException("The dwell times must have same dimension as rates parameter.");
            }
        } else {
            throw new IllegalArgumentException("The trait class type is not suitable for use in this class.");
        }

        if (trait instanceof Model) {
            addModel((Model)trait);
        }
    }

    private DiscreteTraitRateModel(TreeModel treeModel, int traitIndex, Parameter ratesParameter) {
        super(DISCRETE_TRAIT_RATE_MODEL);
        addModel(treeModel);
        this.traitIndex = traitIndex;
        this.ratesParameter = ratesParameter;
        addVariable(ratesParameter);
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        // TreeModel has changed...
        normKnown = false;
        if (model instanceof TreeModel) {
            treeChanged = true;
            shouldRestoreTree = true;
        }
        fireModelChanged();
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // Rate Parameters have changed
        //ratesCalculated = false;
        normKnown = false;
        fireModelChanged();
    }

    protected void storeState() {
        storedNormKnown = normKnown;
        storedNorm = norm;
        shouldRestoreTree = false;
    }

    protected void restoreState() {
        normKnown = storedNormKnown;
        norm = storedNorm;
        treeChanged = shouldRestoreTree;
    }

    protected void acceptState() {
        // nothing to do
    }

    public double getBranchRate(final Tree tree, final NodeRef node) {

//        if (!normKnown) {
//            norm = calculateNorm(tree);
//            normKnown = true;
//        }
//        return getRawBranchRate(tree, node) / norm;

        // AR - I am not sure the normalization is required here?
        return getRawBranchRate(tree, node);
    }

    double calculateNorm(Tree tree) {

        double time = 0.0;
        double rateTime = 0.0;
        for (int i = 0; i < tree.getNodeCount(); i++) {

            NodeRef node = tree.getNode(i);

            if (!tree.isRoot(node)) {

                double branchTime = tree.getBranchLength(node);

                rateTime += getRawBranchRate(tree, node) * branchTime;
                time += branchTime;
            }

        }
        return rateTime / time;
    }

    double getRawBranchRate(final Tree tree, final NodeRef node) {

        double rate = 0.0;
        double[] dwellTimes;
        if (mode == Mode.DWELL_TIMES) {
            dwellTimes = ((double[][])trait.getTrait(tree, node))[0];
        } else {
            dwellTimes = getDwellTimes(tree, node);
        }

        double totalTime = 0;
        for (int i = 0; i < ratesParameter.getDimension(); i++) {
            rate = ratesParameter.getParameterValue(i) * dwellTimes[i];
            totalTime += dwellTimes[i];
        }
        rate /= totalTime;

        return rate;
    }

    private double[] getDwellTimes(final Tree tree, final NodeRef node) {

        double[] dwellTimes = new double[ratesParameter.getDimension()];
        double branchTime = tree.getBranchLength(node);

        if (mode == Mode.PARSIMONY) {
            // an approximation to dwell times using parsimony, assuming
            // the state changes midpoint on the tree. Does a weighted
            // average of the equally parsimonious state reconstructions
            // at the top and bottom of each branch.

            if (treeChanged) {
                fitchParsimony.initialize(tree);
                // Debugging test to count work
//                treeInitializeCounter += 1;
//                if (treeInitializeCounter % 10 == 0) {
//                    System.err.println("Cnt: "+treeInitializeCounter);
//                }
                treeChanged = false;
            }
            int[] states = fitchParsimony.getStates(tree, node);
            int[] parentStates = fitchParsimony.getStates(tree, tree.getParent(node));

            for (int state : states) {
                dwellTimes[state] += branchTime / 2;
            }
            for (int state : parentStates) {
                dwellTimes[state] += branchTime / 2;
            }

            for (int i = 0; i < dwellTimes.length; i++) {
                // normalize by the number of equally parsimonious states at each end of the branch
                // dwellTimes should add up to the total branch length
                dwellTimes[i] /= (states.length + parentStates.length) / 2;
            }
        } else if (mode == Mode.NODE_STATES) {
            // if the states are being sampled - then there is only one possible state at each
            // end of the branch.
            int state = (Integer)trait.getTrait(tree, node)[traitIndex];
            dwellTimes[state] += branchTime / 2;
            int parentState = (Integer)trait.getTrait(tree, node)[traitIndex];
            dwellTimes[parentState] += branchTime / 2;
        }

        return dwellTimes;
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return DISCRETE_TRAIT_RATE_MODEL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

            PatternList patternList = (PatternList) xo.getChild(PatternList.class);

            TreeTraitProvider traitProvider = (TreeTraitProvider) xo.getChild(TreeTraitProvider.class);

            Parameter ratesParameter = (Parameter) xo.getElementFirstChild(RATES);

            int traitIndex = xo.getAttribute(TRAIT_INDEX, 1) - 1;
            String traitName = xo.getAttribute(TRAIT_NAME, "states");

            Logger.getLogger("dr.evomodel").info("Using discrete trait branch rate model.\n" +
                    "\tIf you use this model, please cite:\n" +
                    "\t\tDrummond and Suchard (in preparation)");

            if (traitProvider == null) {
                // Use the version that reconstructs the trait using parsimony:
                return new DiscreteTraitRateModel(treeModel, patternList, traitIndex, ratesParameter);
            } else {
                TreeTrait trait = traitProvider.getTreeTrait(traitName);
                if (trait == null) {
                    throw new XMLParseException("A trait called, " + traitName + ", was not available from the TreeTraitProvider supplied to " + getParserName() + ", with ID " + xo.getId());
                }

                return new DiscreteTraitRateModel(treeModel, trait, traitIndex, ratesParameter);
            }
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return
                    "This Branch Rate Model takes a discrete trait reconstruction (provided by a TreeTraitProvider) and " +
                            "gives the rate for each branch of the tree based on the child trait of " +
                            "that branch. The rates for each trait value are specified in a multidimensional parameter.";
        }

        public Class getReturnType() {
            return DiscreteTraitRateModel.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(TreeModel.class, "The tree model"),
                new XORRule(
                        new ElementRule(TreeTraitProvider.class, "The trait provider"),
                        new ElementRule(PatternList.class)),
                new ElementRule(RATES, Parameter.class, "The rates of the different trait values", false),
                AttributeRule.newIntegerRule(TRAIT_INDEX, true),
                AttributeRule.newStringRule(TRAIT_NAME, true)
        };
    };

}
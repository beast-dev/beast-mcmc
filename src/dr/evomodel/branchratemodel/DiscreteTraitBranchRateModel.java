/*
 * DiscreteTraitBranchRateModel.java
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

package dr.evomodel.branchratemodel;

import dr.evolution.alignment.PatternList;
import dr.evolution.datatype.DataType;
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

/**
 * This Branch Rate Model takes a ancestral state likelihood and
 * gives the rate for each branch of the tree based on the ancestor
 * state and child state or based on the Markov Rewards (i.e., time
 * spent in each state).
 *
 * @author Alexei Drummond
 * @author Marc Suchard
 * @author Andrew Rambaut
 * @author Trevor Bedford
 */
public class DiscreteTraitBranchRateModel extends AbstractBranchRateModel {
    private static final boolean CACHING_RATES = true;

    enum Mode {
        NODE_STATES,
        MARKOV_JUMP_PROCESS,
        MARKOV_JUMP_COUNT,
        PARSIMONY
    }

    public static final String DISCRETE_TRAIT_BRANCH_RATE_MODEL = "discreteTraitRateModel";

    protected TreeTrait trait = null;
    private Parameter rateParameter;
    private Parameter relativeRatesParameter;
    private Parameter indicatorParameter;
    protected int traitIndex;

    private double[] rates;
    private double[] storedRates;
    private boolean[] rateKnown;

//    private boolean normKnown = false;
//    private boolean storedNormKnown = false;
//    private double norm = 1.0;
//    private double storedNorm = 1.0;

    private TreeTrait[] traits;

    private FitchParsimony fitchParsimony;

    private boolean treeChanged = true;

    private Mode mode;
    private DataType dataType;
//    private int treeInitializeCounter = 0;

    /**
     * A constructor for the (crude) parsimony reconstruction form of this class.
     * @param treeModel
     * @param patternList
     * @param traitIndex
     * @param ratesParameter
     */
    public DiscreteTraitBranchRateModel(TreeModel treeModel, PatternList patternList, int traitIndex, Parameter ratesParameter) {

        this(treeModel, traitIndex, ratesParameter, null, null);

        if (!TaxonList.Utils.getTaxonListIdSet(treeModel).equals(TaxonList.Utils.getTaxonListIdSet(patternList))) {
            throw new IllegalArgumentException("Tree model and pattern list must have the same list of taxa!");
        }

        ratesParameter.setDimension(patternList.getDataType().getStateCount());

        fitchParsimony = new FitchParsimony(patternList, false);
        mode = Mode.PARSIMONY;
    }

    /**
     * A constructor for a node-sampled discrete trait
     * @param treeModel
     * @param trait
     * @param traitIndex
     * @param rateParameter
     * @param relativeRatesParameter
     * @param indicatorParameter
     */
    public DiscreteTraitBranchRateModel(TreeTraitProvider traitProvider, DataType dataType, TreeModel treeModel,
                                        TreeTrait trait, int traitIndex, Parameter rateParameter, Parameter relativeRatesParameter, Parameter indicatorParameter) {

        this(treeModel, traitIndex, rateParameter, relativeRatesParameter, indicatorParameter);

//        if (trait.getTreeModel() != treeModel)
//            throw new IllegalArgumentException("Tree Models for ancestral state tree likelihood and target model of these rates must match!");

        this.trait = trait;
        this.dataType = dataType;

        if (trait.getTraitName().equals("states")) {
            // Assume the trait is one or more discrete traits reconstructed at nodes
            mode = Mode.NODE_STATES;
        } else /*if (double[].class.isAssignableFrom(trait.getClass()))*/ {
            // Assume the trait itself is the dwell times for the individual states on the branch above the node
            mode = Mode.MARKOV_JUMP_PROCESS;
        } /* else {
            throw new IllegalArgumentException("The trait class type is not suitable for use in this class.");
        } */

        relativeRatesParameter.setDimension(dataType.getStateCount());

        if (traitProvider instanceof Model) {
            addModel((Model)traitProvider);
        }

        if (trait instanceof Model) {
            addModel((Model)trait); // MAS: Does this ever occur?
        }
    }

    public DiscreteTraitBranchRateModel(TreeTraitProvider traitProvider, DataType dataType, TreeModel treeModel,
                                        TreeTrait trait, int traitIndex, Parameter ratesParameter) {

        this(treeModel, traitIndex, ratesParameter, null, null);

//        if (trait.getTreeModel() != treeModel)
//            throw new IllegalArgumentException("Tree Models for ancestral state tree likelihood and target model of these rates must match!");

        this.trait = trait;
        this.dataType = dataType;

        if (trait.getTraitName().equals("states")) {
            // Assume the trait is one or more discrete traits reconstructed at nodes
            mode = Mode.NODE_STATES;
        } else /*if (double[].class.isAssignableFrom(trait.getClass()))*/ {
            // Assume the trait itself is the dwell times for the individual states on the branch above the node
            mode = Mode.MARKOV_JUMP_PROCESS;
        } /* else {
            throw new IllegalArgumentException("The trait class type is not suitable for use in this class.");
        } */

        ratesParameter.setDimension(dataType.getStateCount());

        if (traitProvider instanceof Model) {
            addModel((Model)traitProvider);
        }

        if (trait instanceof Model) {
            addModel((Model)trait); // MAS: Does this ever occur?
        }
    }

    public DiscreteTraitBranchRateModel(TreeTraitProvider traitProvider, TreeTrait[] traits, TreeModel treeModel, Parameter ratesParameter) {

        this(treeModel, 0, ratesParameter, null, null);

        this.traits = traits;
        mode = Mode.MARKOV_JUMP_PROCESS;

        ratesParameter.setDimension(traits.length);


        if (traitProvider instanceof Model) {
            addModel((Model)traitProvider);
        }
    }

    private DiscreteTraitBranchRateModel(TreeModel treeModel, int traitIndex,
                                         Parameter rateParameter, Parameter relativeRatesParameter, Parameter indicatorParameter) {
        super(DISCRETE_TRAIT_BRANCH_RATE_MODEL);
        addModel(treeModel);
        this.traitIndex = traitIndex;

        this.rateParameter = rateParameter;
        addVariable(rateParameter);

        this.relativeRatesParameter = relativeRatesParameter;
        if (relativeRatesParameter != null) {
            addVariable(relativeRatesParameter);
        }

        this.indicatorParameter = indicatorParameter;
        if (indicatorParameter != null) {
            addVariable(indicatorParameter);
        }

        rates = new double[treeModel.getNodeCount()];
        storedRates = new double[treeModel.getNodeCount()];
        rateKnown = new boolean[treeModel.getNodeCount()];

    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        // TreeModel has changed...
        for (int i = 0; i < rateKnown.length; i++) {
            rateKnown[i] = false;
        }
        treeChanged = true;
        fireModelChanged();
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // Rate Parameters have changed
        for (int i = 0; i < rateKnown.length; i++) {
            rateKnown[i] = false;
        }
        fireModelChanged();
    }

    protected void storeState() {
        if (CACHING_RATES) {
            System.arraycopy(rates, 0, storedRates, 0, rates.length);
        }
    }

    protected void restoreState() {
        if (CACHING_RATES) {
            double[] tmp = rates;
            rates = storedRates;
            storedRates = tmp;

            for (int i = 0; i < rateKnown.length; i++) {
                rateKnown[i] = true;
            }
        } else {
            for (int i = 0; i < rateKnown.length; i++) {
                rateKnown[i] = false;
            }
            treeChanged = true;
        }
    }

    protected void acceptState() {
        // nothing to do
    }

    protected int getStateCount() {
        int dimen = 0;
        if (mode == Mode.NODE_STATES || mode == Mode.MARKOV_JUMP_PROCESS) {
            dimen = dataType.getStateCount();
        }
        else if (mode == Mode.PARSIMONY) {
            dimen = fitchParsimony.getPatterns().getStateCount();
        }
        return dimen;
    }

    public double getBranchRate(final Tree tree, final NodeRef node) {

        if (CACHING_RATES) {
            if (!rateKnown[node.getNumber()]) {
                rates[node.getNumber()] = getRawBranchRate(tree, node);
                rateKnown[node.getNumber()] = true;
            }

            return rates[node.getNumber()];
        } else {
            return getRawBranchRate(tree, node);
        }
    }

    // produce weighted mean of rate for a branch
    // rate = absRate * branchWeight[0] * relativeRates[0] + absRate * branchWeight[1] * relativeRates[1]
    protected double getRawBranchRate(final Tree tree, final NodeRef node) {

        double rate = 0.0;
        int stateCount = getStateCount();

        double[] processValues = getProcessValues(tree, node);
        //    double[] processValues = {1.0, 1.0};
        double[] branchWeights = new double[stateCount];
        double totalTime = 0;

        for (int i = 0; i < stateCount; i++) {
            branchWeights[i] += processValues[i];
            totalTime += processValues[i];
        }

        for (int i = 0; i < stateCount; i++) {
            branchWeights[i] /= totalTime;
        }

        if (relativeRatesParameter != null && indicatorParameter == null) {
            double absRate = rateParameter.getParameterValue(0);
            for (int i = 0; i < stateCount; i++) {
                rate += absRate * relativeRatesParameter.getParameterValue(i) * branchWeights[i];
            }
        }
        else if (relativeRatesParameter != null && indicatorParameter != null) {
            double absRate = rateParameter.getParameterValue(0);
            for (int i = 0; i < stateCount; i++) {
                rate += absRate * relativeRatesParameter.getParameterValue(i) * branchWeights[i] * (double) indicatorParameter.getParameterValue(i);
            }
        } else {
            for (int i = 0; i < stateCount; i++) {
                rate += rateParameter.getParameterValue(i) * processValues[i];
                totalTime += processValues[i];
            }
        }

        return rate;
    }

    /**
     *
     * @param tree
     * @param node
     * @return and array of the total amount of time spent in each of the discrete states along the branch above the given node.
     */
    private double[] getProcessValues(final Tree tree, final NodeRef node) {

        double[] processValues = null;
        int stateCount = getStateCount();
        double branchTime = tree.getBranchLength(node);

        if (mode == Mode.MARKOV_JUMP_PROCESS) {
            processValues = new double[stateCount];
            for (int i = 0; i < stateCount; i++) {
                processValues[i] = ((TreeTrait.DA)traits[i]).getTrait(tree, node)[0];
            }
        } else if (mode == Mode.PARSIMONY) {
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

            processValues = new double[fitchParsimony.getPatterns().getStateCount()];

            for (int state : states) {
                processValues[state] += branchTime / 2;
            }
            for (int state : parentStates) {
                processValues[state] += branchTime / 2;
            }

            for (int i = 0; i < processValues.length; i++) {
                // normalize by the number of equally parsimonious states at each end of the branch
                // processValues should add up to the total branch length
                processValues[i] /= (states.length + parentStates.length) / 2;
            }
        } else if (mode == Mode.NODE_STATES) {
            processValues = new double[stateCount];

            // if the states are being sampled - then there is only one possible state at each
            // end of the branch.
            int state = ((int[])trait.getTrait(tree, node))[traitIndex];
            processValues[state] += branchTime / 2;
            NodeRef parent = tree.getParent(node);
            int parentState = ((int[])trait.getTrait(tree, parent))[traitIndex];
            processValues[parentState] += branchTime / 2;
        }

        return processValues;
    }

}
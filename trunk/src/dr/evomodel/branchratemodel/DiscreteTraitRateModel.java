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
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treelikelihood.AncestralStateTreeLikelihood;
import dr.inference.model.AbstractModel;
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

    public static final String DISCRETE_TRAIT_RATE_MODEL = "discreteTraitRateModel";
    public static final String RATES = "rates";
    public static final String TRAIT_INDEX = "traitIndex";

    private AncestralStateTreeLikelihood ancestors;
    private Parameter ratesParameter;
    private int traitIndex;
    private boolean normKnown = false;
    private boolean storedNormKnown = false;
    private double norm = 1.0;
    private double storedNorm = 1.0;

    private FitchParsimony fitchParsimony;

    private boolean treeChanged = true;
    private boolean shouldRestoreTree = false;

    private boolean fitch;

//    private int treeInitializeCounter = 0;

    public DiscreteTraitRateModel(TreeModel treeModel, PatternList patternList, int traitIndex, Parameter ratesParameter) {

        this(treeModel, traitIndex, ratesParameter);

        if (!TaxonList.Utils.getTaxonListIdSet(treeModel).equals(TaxonList.Utils.getTaxonListIdSet(patternList))) {
            throw new IllegalArgumentException("Tree model and pattern list must have the same list of taxa!");
        }

        fitchParsimony = new FitchParsimony(patternList, false);
        fitch = true;
    }

    private DiscreteTraitRateModel(TreeModel treeModel, int traitIndex, Parameter ratesParameter) {
        super(DISCRETE_TRAIT_RATE_MODEL);
        addModel(treeModel);
        this.traitIndex = traitIndex;
        this.ratesParameter = ratesParameter;
        addVariable(ratesParameter);
    }

    public DiscreteTraitRateModel(TreeModel treeModel, AncestralStateTreeLikelihood like, int traitIndex, Parameter ratesParameter) {

        this(treeModel, traitIndex, ratesParameter);

        if (like.getTreeModel() != treeModel)
            throw new IllegalArgumentException("Tree Models for ancestral state tree likelihood and target model of these rates must match!");

        ancestors = like;

        if (like.getDataType().getStateCount() != ratesParameter.getDimension())
            throw new IllegalArgumentException("ancestral likelihood datatype must have same size as rates parameter!");

        addModel(like);

        fitch = false;
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

        if (!normKnown) {
            norm = calculateNorm(tree);
            normKnown = true;
        }
        return getRawBranchRate(tree, node) / norm;
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
        if (fitch) {
            if (treeChanged) {
                fitchParsimony.initialize(tree);
                // Debugging test to count work
//                treeInitializeCounter += 1;
//                if (treeInitializeCounter % 10 == 0) {
//                    System.err.println("Cnt: "+treeInitializeCounter);
//                }
                treeChanged = false;
            }
            int[] traits = fitchParsimony.getStates(tree, node);
            for (int trait : traits) {
                rate += ratesParameter.getParameterValue(trait);
            }
            rate /= traits.length;

        } else {
            int trait = ancestors.getStatesForNode(tree, node)[traitIndex];
            rate = ratesParameter.getParameterValue(trait);
        }
        return rate;
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return DISCRETE_TRAIT_RATE_MODEL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

            PatternList patternList = (PatternList) xo.getChild(PatternList.class);

            AncestralStateTreeLikelihood ancestors = (AncestralStateTreeLikelihood) xo.getChild(AncestralStateTreeLikelihood.class);

            Parameter ratesParameter = (Parameter) xo.getElementFirstChild(RATES);

            int traitIndex = xo.getAttribute(TRAIT_INDEX, 1) - 1;

            Logger.getLogger("dr.evomodel").info("Using discrete trait branch rate model.\n" +
                    "\tIf you use this model, please cite:\n" +
                    "\t\tDrummond and Suchard (in preparation)");

            if (ancestors == null) {
                return new DiscreteTraitRateModel(treeModel, patternList, traitIndex, ratesParameter);
            } else {
                return new DiscreteTraitRateModel(treeModel, ancestors, traitIndex, ratesParameter);
            }
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return
                    "This Branch Rate Model takes a discrete trait reconstruction (provided by AncestralStateTreeLikelihood) and " +
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
                        new ElementRule(AncestralStateTreeLikelihood.class, "The colour sampler model"),
                        new ElementRule(PatternList.class)),
                new ElementRule(RATES, Parameter.class, "The rates of the different trait values", false),
                AttributeRule.newIntegerRule(TRAIT_INDEX, true)
        };
    };

}
/*
 * RandomLocalYuleModel.java
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

package dr.evomodel.speciation;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.randomlocalmodel.RandomLocalTreeVariable;
import dr.evomodelxml.speciation.RandomLocalYuleModelParser;
import dr.inference.model.Parameter;

import java.text.NumberFormat;
import java.util.Locale;
import java.util.logging.Logger;

/**
 * This class contains methods that describe a Yule speciation model whose rate of birth changes
 * at different points in the tree.
 *
 * @author Alexei Drummond
 */
public class RandomLocalYuleModel extends UltrametricSpeciationModel implements TreeTraitProvider, RandomLocalTreeVariable {
    private boolean calculateAllBirthRates = false;

    public RandomLocalYuleModel(Parameter birthRates, Parameter indicators, Parameter meanRate,
                                boolean ratesAsMultipliers, Type units, int dp) {

        super(RandomLocalYuleModelParser.YULE_MODEL, units);

        addVariable(birthRates);
        birthRates.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, birthRates.getDimension()));

        for (int i = 0; i < indicators.getDimension(); i++) {
            indicators.setParameterValueQuietly(i, 0.0);
        }

        addVariable(indicators);

        this.meanRate = meanRate;
        addVariable(meanRate);

        birthRatesAreMultipliers = ratesAsMultipliers;

        format.setMaximumFractionDigits(dp);

        birthRatesName = birthRates.getParameterName();
        Logger.getLogger("dr.evomodel").info("  birth rates parameter is named '" + birthRatesName + "'");
        indicatorsName = indicators.getParameterName();
        Logger.getLogger("dr.evomodel").info("  indicator parameter is named '" + indicatorsName + "'");

        this.birthRates = new double[birthRates.getDimension() + 1];

        treeTraits.addTrait(new TreeTrait.I() {
            public String getTraitName() {
                return "I";
            }

            public Intent getIntent() {
                return Intent.NODE;
            }

            public Integer getTrait(Tree tree, NodeRef node) {
                return (isVariableSelected(tree, node) ? 1 : 0);
            }

        });

        treeTraits.addTrait(new TreeTrait.D() {
            public String getTraitName() {
                return "b";
            }

            public Intent getIntent() {
                return Intent.NODE;
            }

            public Double getTrait(Tree tree, NodeRef node) {
                return RandomLocalYuleModel.this.birthRates[node.getNumber()];
            }

        });
    }

    public final double getVariable(Tree tree, NodeRef node) {
        return ((TreeModel)tree).getNodeTrait(node, birthRatesName);
    }

    public final boolean isVariableSelected(Tree tree, NodeRef node) {
        return ((TreeModel)tree).getNodeTrait(node, indicatorsName) > 0.5;
    }

    //
    // functions that define a speciation model
    //
    public final double logTreeProbability(int taxonCount) {

        // calculate all nodes birth rates
        calculateAllBirthRates = true;

        return 0.0;

    }

    //
    // functions that define a speciation model
    //
    public final double logNodeProbability(Tree tree, NodeRef node) {

        if (calculateAllBirthRates) {
            calculateBirthRates((TreeModel) tree, tree.getRoot(), 0.0);
            calculateAllBirthRates = false;
        }

        if (tree.isRoot(node)) {
            return 0.0;
        } else {

            double lambda = birthRates[node.getNumber()];
            double branchLength = tree.getNodeHeight(tree.getParent(node)) - tree.getNodeHeight(node);
            double logP = -lambda * branchLength;

            if (tree.isExternal(node)) logP += Math.log(lambda);

            return logP;
        }
    }

    private void calculateBirthRates(TreeModel tree, NodeRef node, double rate) {

        int nodeNumber = node.getNumber();

        if (tree.isRoot(node)) {
            rate = meanRate.getParameterValue(0);
        } else {
            if (isVariableSelected(tree, node)) {
                if (birthRatesAreMultipliers) {
                    rate *= getVariable(tree, node);
                } else {
                    rate = getVariable(tree, node);
                }
            }
        }
        birthRates[nodeNumber] = rate;

        int childCount = tree.getChildCount(node);
        for (int i = 0; i < childCount; i++) {
            calculateBirthRates(tree, tree.getChild(node, i), rate);
        }
    }

//    /**
//     * @param tree the tree
//     * @param node the node to retrieve the birth rate of
//     * @return the birth rate of the given node;
//     */
//    private double getBirthRate(TreeModel tree, NodeRef node) {
//
//        double birthRate;
//        if (!tree.isRoot(node)) {
//
//            double parentRate = getBirthRate(tree, tree.getParent(node));
//            if (isVariableSelected(tree, node)) {
//                birthRate = getVariable(tree, node);
//                if (birthRatesAreMultipliers) {
//                    birthRate *= parentRate;
//                } else {
//                    throw new RuntimeException("Rates must be multipliers in current implementation! " +
//                            "Otherwise root rate might be ignored");
//                }
//            } else {
//                birthRate = parentRate;
//            }
//        } else {
//            birthRate = meanRate.getParameterValue(0);
//        }
//        return birthRate;
//    }

    protected TreeTraitProvider.Helper treeTraits = new Helper();

    public TreeTrait[] getTreeTraits() {
        return treeTraits.getTreeTraits();
    }

    public TreeTrait getTreeTrait(String key) {
        return treeTraits.getTreeTrait(key);
    }

    public boolean includeExternalNodesInLikelihoodCalculation() {
        return true;
    }

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

    public org.w3c.dom.Element createElement(org.w3c.dom.Document d) {
        throw new RuntimeException("createElement not implemented");
    }

    private double[] birthRates;

    private String birthRatesName = "birthRates";
    private String indicatorsName = "birthRateIndicator";

    private Parameter meanRate;
    private boolean birthRatesAreMultipliers = false;
    private NumberFormat format = NumberFormat.getNumberInstance(Locale.ENGLISH);
}

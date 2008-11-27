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

import dr.evolution.colouring.BranchColouring;
import dr.evolution.colouring.TreeColouring;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.coalescent.structure.ColourSamplerModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * This Branch Rate Model takes a tree colouring (provided by a ColourSamplerModel) and
 * gives the rate for each branch of the tree based on the time spent in each colour along
 * that branch. The rates for each colour are specified in a multidimensional parameter.
 *
 * @author Andrew Rambaut
 * @version $Id: ColouredTreeRateModel.java,v 1.2 2006/01/10 16:48:27 rambaut Exp $
 */
public class ColouredTreeRateModel extends AbstractModel implements BranchRateModel {
    public static final String COLOURED_TREE_RATE_MODEL = "colouredTreeRateModel";
    public static final String SUBSTITUTION_RATES = "substitutionRates";

    private final ColourSamplerModel colourSamplerModel;
    private final Parameter substitutionRatesParameter;

    private final double[] rates;
    private boolean ratesCalculated = false;

    public ColouredTreeRateModel(TreeModel treeModel, ColourSamplerModel colourSamplerModel, Parameter substitutionRatesParameter) {

        super(COLOURED_TREE_RATE_MODEL);

        rates = new double[treeModel.getNodeCount()];

        addModel(treeModel);

        this.colourSamplerModel = colourSamplerModel;
        addModel(colourSamplerModel);

        this.substitutionRatesParameter = substitutionRatesParameter;

        addParameter(substitutionRatesParameter);
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        // TreeModel or ColourSamplerModel has changed...
        ratesCalculated = false;
        fireModelChanged();
    }

    protected final void handleParameterChangedEvent(Parameter parameter, int index, Parameter.ChangeType type) {
        // Rate Parameters have changed
        ratesCalculated = false;
        fireModelChanged();
    }

    protected void storeState() {
        // nothing to do
    }

    protected void restoreState() {
        ratesCalculated = false;
        // nothing to do
    }

    protected void acceptState() {
        // nothing to do
    }

    public double getBranchRate(Tree tree, NodeRef node) {

        if (!ratesCalculated) {
            TreeColouring colouring = colourSamplerModel.getTreeColouring();
            double[] substitutionRates = substitutionRatesParameter.getParameterValues();
            calculateNodeRates(tree, tree.getRoot(), colouring, substitutionRates);
            ratesCalculated = true;
        }

        return rates[node.getNumber()];
    }

    public String getBranchAttributeLabel() {
        return "rate";
    }

    public String getAttributeForBranch(Tree tree, NodeRef node) {
        return Double.toString(getBranchRate(tree, node));
    }

    /**
     * Traverse the tree calculating rates.
     */
    private final void calculateNodeRates(Tree tree, NodeRef node, TreeColouring colouring, double[] substitutionRates) {

        if (!tree.isExternal(node)) {

            // Traverse down the two child nodes
            NodeRef child1 = tree.getChild(node, 0);
            calculateNodeRates(tree, child1, colouring, substitutionRates);

            NodeRef child2 = tree.getChild(node, 1);
            calculateNodeRates(tree, child2, colouring, substitutionRates);

        }

        NodeRef parent = tree.getParent(node);

        // don't bother doing anything further at the root because rate at root is ignored
        if (parent == null) return;

        BranchColouring branchColouring = colouring.getBranchColouring(node);

        double totalTime = tree.getBranchLength(node);

        double rate = 0.0;

        for (int i = 0; i < colouring.getColourCount(); i++) {
            double colourTime = branchColouring.getTimeInColour(i, tree.getNodeHeight(parent), tree.getNodeHeight(node));
            rate += substitutionRates[i] * colourTime;

        }

        rates[node.getNumber()] = rate / totalTime;
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return COLOURED_TREE_RATE_MODEL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

            ColourSamplerModel colourSamplerModel = (ColourSamplerModel) xo.getChild(ColourSamplerModel.class);

            Parameter substitutionRatesParameter = (Parameter) xo.getElementFirstChild(SUBSTITUTION_RATES);

            Logger.getLogger("dr.evomodel").info("Using coloured tree clock model.");

            return new ColouredTreeRateModel(treeModel, colourSamplerModel, substitutionRatesParameter);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return
                    "This Branch Rate Model takes a tree colouring (provided by a ColourSamplerModel) and " +
                            "gives the rate for each branch of the tree based on the time spent in each colour along " +
                            "that branch. The rates for each colour are specified in a multidimensional parameter.";
        }

        public Class getReturnType() {
            return ColouredTreeRateModel.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(TreeModel.class, "The tree model"),
                new ElementRule(ColourSamplerModel.class, "The colour sampler model"),
                new ElementRule(SUBSTITUTION_RATES, Parameter.class, "The substitution rates of the different colours", false)
        };
    };

}
/*
 * TraitBranchRateModel.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
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
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * Takes the log rates at each node provided by a specified rate and give the branch rate as the average.
 *
 * @author Andrew Rambaut
 */
public class TraitBranchRateModel extends AbstractModel implements BranchRateModel {


    public static final String TRAIT_BRANCH_RATES = "traitBranchRates";
    public static final String TRAIT = "trait";
    public static final String DIMENSION = "dimension";
    public static final String RATE = "rate";
    public static final String RATIO = "ratio";

    private final String trait;
    private final int dimension;
    private final Parameter rateParameter;
    private final Parameter ratioParameter;

    public TraitBranchRateModel(String trait, int dimension) {
        super(TRAIT_BRANCH_RATES);

        this.trait = trait;
        this.dimension = dimension;

        this.rateParameter = null;
        this.ratioParameter = null;
    }

    public TraitBranchRateModel(String trait, Parameter rateParameter, Parameter ratioParameter) {
        super(TRAIT_BRANCH_RATES);

        this.trait = trait;
        dimension = 0;
        this.rateParameter = rateParameter;
        this.ratioParameter = ratioParameter;

        if (rateParameter != null) {
            addVariable(rateParameter);
        }

        if (ratioParameter != null) {
            addVariable(ratioParameter);
        }
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        fireModelChanged();
    }

    protected void storeState() {
    }

    protected void restoreState() {
    }

    protected void acceptState() {
    }


    public double getBranchRate(Tree tree, NodeRef node) {
        NodeRef parent = tree.getParent(node);
        if (parent == null) {
            throw new IllegalArgumentException("Root does not have a valid rate");
        }

        double rate = 1.0;
        TreeModel treeModel = (TreeModel) tree;

        if (rateParameter != null) {
            double scale = 1.0;
            double ratio = 1.0;

            if (rateParameter != null) {
                scale = rateParameter.getParameterValue(0);
            }

            if (ratioParameter != null) {
                ratio = ratioParameter.getParameterValue(0);
            }


            // get the log rate for the node and its parent
            double rate1 = ratio * treeModel.getMultivariateNodeTrait(node, trait)[0];
            double rate2 = ratio * treeModel.getMultivariateNodeTrait(parent, trait)[0];

            if (rate1 == rate2) {
                return scale * Math.exp(rate1);
            }

            rate = scale * (Math.exp(rate2) - Math.exp(rate1)) / (rate2 - rate1);
        } else {
            double rate1 =  treeModel.getMultivariateNodeTrait(node, trait)[dimension];
            double rate2 =  treeModel.getMultivariateNodeTrait(parent, trait)[dimension];

            if (rate1 == rate2) {
                return Math.exp(rate1);
            }

            rate = (Math.exp(rate2) - Math.exp(rate1)) / (rate2 - rate1);
        }

        return rate;
    }

    public String getBranchAttributeLabel() {
        return "rate";
    }

    public String getAttributeForBranch(Tree tree, NodeRef node) {
        return Double.toString(getBranchRate(tree, node));
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return TRAIT_BRANCH_RATES;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String trait = xo.getStringAttribute(TRAIT);
            int dimension = 0;
            if (xo.hasAttribute(DIMENSION)) {
                dimension = xo.getIntegerAttribute(DIMENSION) - 1;
            }


            Logger.getLogger("dr.evomodel").info("Using trait, " + trait + ", as log rate estimates.");

            if (xo.hasChildNamed(RATE)) {
                Parameter rateParameter = (Parameter) xo.getElementFirstChild(RATE);
                Parameter ratioParameter = (Parameter) xo.getElementFirstChild(RATIO);

                return new TraitBranchRateModel(trait, rateParameter, ratioParameter);
            } else {
                return new TraitBranchRateModel(trait, dimension);
            }
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return
                    "This element returns an trait rate model." +
                            "The branch rates are an average of the rates provided by a node trait.";
        }

        public Class getReturnType() {
            return TraitBranchRateModel.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newStringRule(TRAIT, false, "The name of the trait that provides the log rates at nodes"),
                AttributeRule.newIntegerRule(DIMENSION, true, "The dimension that supplies the rate"),
                new ElementRule(RATE, Parameter.class, "The rate parameter", true),
                new ElementRule(RATIO, Parameter.class, "The ratio parameter", true),
        };
    };


}
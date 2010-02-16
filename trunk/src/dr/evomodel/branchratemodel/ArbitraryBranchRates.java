/*
 * ArbitraryBranchRates.java
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
import dr.evomodel.tree.TreeParameterModel;
import dr.inference.model.Parameter;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Variable;
import dr.xml.*;

import java.util.logging.Logger;

/**
 * Allows branch rates to take on any double value
 * This is useful for forming a scaled mixture of normals for the continuous diffusion model
 *
 * @author Marc A. Suchard
 * @author Alexei Drummond
 */
public class ArbitraryBranchRates extends AbstractModel implements BranchRateModel {

    public static final String ARBITRARY_BRANCH_RATES = "arbitraryBranchRates";
    public static final String RATES = "rates";
    public static final String RECIPROCAL = "reciprocal";

    // The rates of each branch
    final TreeParameterModel rates;
    final Parameter rateParameter;
    final boolean reciprocal;

    public ArbitraryBranchRates(TreeModel tree, Parameter rateParameter, boolean reciprocal) {

        super(ARBITRARY_BRANCH_RATES);

        for(int i = 0; i < rateParameter.getDimension(); i++) {
            rateParameter.setValue(i, 1.0);
        }
        //Force the boundaries of rate
        Parameter.DefaultBounds bound = new Parameter.DefaultBounds(Double.MAX_VALUE, 0, rateParameter.getDimension());
        rateParameter.addBounds(bound);

        this.rates = new TreeParameterModel(tree, rateParameter, false);
        this.rateParameter = rateParameter;

        addModel(rates);

        this.reciprocal = reciprocal;
    }

    public void setBranchRate(Tree tree, NodeRef node, double value) {
        rates.setNodeValue(tree, node, value);
    }
   
    public double getBranchRate(Tree tree, NodeRef node) {
        // Branch rates are proportional to time.
        // In the traitLikelihoods, time is proportional to variance
        // Fernandez and Steel (2000) shows the sampling density with the scalar proportional to precision 
        final double rate = rates.getNodeValue(tree,node);
        if (reciprocal) {
            return 1.0 / rate;
        }
        return rate;
    }

    public boolean usingReciprocal() {
        return reciprocal;
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        // Should be called by TreeParameterModel
        fireModelChanged(null, index);
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // Changes to rateParameter are handled by model changed events
    }

    protected void storeState() {
    }

    protected void restoreState() {
    }

    protected void acceptState() {
    }

    public String getBranchAttributeLabel() {
        return "rate";
    }

    public String getAttributeForBranch(Tree tree, NodeRef node) {
        return Double.toString(getBranchRate(tree, node));
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return ARBITRARY_BRANCH_RATES;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            TreeModel tree = (TreeModel) xo.getChild(TreeModel.class);

            XMLObject cxo = xo.getChild(RATES);

            Parameter rateCategoryParameter = (Parameter) cxo.getChild(Parameter.class);

            boolean reciprocal = xo.getAttribute(RECIPROCAL,false);

            Logger.getLogger("dr.evomodel").info("Using an scaled mixture of normals model.");
            Logger.getLogger("dr.evomodel").info("  rates = " + rateCategoryParameter.getDimension());
            Logger.getLogger("dr.evomodel").info("  NB: Make sure you have a prior on " + rateCategoryParameter.getId() + " and do not use this model in a treeLikelihood");
            Logger.getLogger("dr.evomodel").info("  reciprocal = "+reciprocal);

            return new ArbitraryBranchRates(tree, rateCategoryParameter, reciprocal);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns an arbitrary rate model." +
                   "The branch rates are drawn from an arbitrary distribution determine by the prior.";
        }

        public Class getReturnType() {
            return ArbitraryBranchRates.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(TreeModel.class),
                new ElementRule(RATES, Parameter.class, "The rate parameter"),
                AttributeRule.newBooleanRule(RECIPROCAL,true),
        };
    };
}

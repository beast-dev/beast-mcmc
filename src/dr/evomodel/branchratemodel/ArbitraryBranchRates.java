/*
 * ArbitraryBranchRates.java
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

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.evomodelxml.branchratemodel.ArbitraryBranchRatesParser;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * Allows branch rates to take on any double value
 * This is useful for forming a scaled mixture of normals for the continuous diffusion model
 *
 * @author Marc A. Suchard
 * @author Alexei Drummond
 */
public class ArbitraryBranchRates extends AbstractBranchRateModel {

    // The rates of each branch
    final TreeParameterModel rates;
    final Parameter rateParameter;
    final boolean reciprocal;
    final boolean exp;

    public ArbitraryBranchRates(TreeModel tree, Parameter rateParameter, boolean reciprocal, boolean exp, boolean setRates) {

        super(ArbitraryBranchRatesParser.ARBITRARY_BRANCH_RATES);

        if (setRates) {
            double value = exp ? 0.0 : 1.0;
            for (int i = 0; i < rateParameter.getDimension(); i++) {
                rateParameter.setValue(i, value);
            }
        }
        //Force the boundaries of rate
        Parameter.DefaultBounds bounds = exp ?
                new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY,
                        rateParameter.getDimension()) :
                new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0,
                        rateParameter.getDimension());
        rateParameter.addBounds(bounds);

        this.rates = new TreeParameterModel(tree, rateParameter, false);
        this.rateParameter = rateParameter;

        addModel(rates);

        this.reciprocal = reciprocal;
        this.exp = exp;
    }

    public void setBranchRate(Tree tree, NodeRef node, double value) {
        rates.setNodeValue(tree, node, value);
    }

    public double getBranchRate(final Tree tree, final NodeRef node) {
        // Branch rates are proportional to time.
        // In the traitLikelihoods, time is proportional to variance
        // Fernandez and Steel (2000) shows the sampling density with the scalar proportional to precision 
        double rate = rates.getNodeValue(tree, node);

        if (exp) {
            rate = Math.exp(rate);
        }

        if (reciprocal) {
            rate = 1.0 / rate;
        }

        return rate;
    }

    public boolean usingReciprocal() {
        return reciprocal;
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        // Should be called by TreeParameterModel
        fireModelChanged(object, index);
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

}

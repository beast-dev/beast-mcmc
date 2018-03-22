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
    private final TreeParameterModel rates;
    private final Parameter rateParameter;
    private final boolean reciprocal;
    private final boolean exp;

    private final Parameter locationParameter;
    private final Parameter scaleParameter;

    public ArbitraryBranchRates(TreeModel tree, Parameter rateParameter,
                                boolean reciprocal, boolean exp, boolean setRates) {
        this(tree, rateParameter, reciprocal, exp, setRates, null, null);
    }

    public ArbitraryBranchRates(TreeModel tree, Parameter rateParameter,
                                boolean reciprocal, boolean exp, boolean setRates, // TODO No multiple-boolean arguments
                                Parameter locationParameter, Parameter scaleParameter) {

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

        this.locationParameter = locationParameter;
        this.scaleParameter = scaleParameter;

        if (locationParameter != null) {
            addVariable(locationParameter);
        }

        if (scaleParameter != null) {
            addVariable(scaleParameter);
        }

        if ((reciprocal || exp) && (locationParameter != null || scaleParameter != null)) {
            throw new RuntimeException("Not yet implemented");
        }
    }

    public void setBranchRate(Tree tree, NodeRef node, double value) {
        rates.setNodeValue(tree, node, value);
    }

    public double getBranchRateDifferential(final Tree tree, final NodeRef node) {

        double differential = exp ? 1.0 : getBranchRate(tree, node);
        
        if (reciprocal) {
            differential = -differential;
        } else {
            differential = 1 / differential;
        }
        
        return differential;
    }

    private final double baseMeasureMu = getMuPhi(1.0); // TODO Depends on prior distribution
    private final double baseMeasureSigma = getSigmaPhi(1.0); // TODO Depends on prior distribution

    private double transformMu;
    private double transformSigma;

    private boolean transformKnown = false;

    private void setupTransform() {
        final double phi = scaleParameter.getParameterValue(0) * scaleParameter.getParameterValue(0);
        transformMu = getMuPhi(phi);
        transformSigma = getSigmaPhi(phi);
    }

    private double transform(double rate) {

        // TODO Delegate depending on distribution
        if (exp) {

            if (scaleParameter != null) {
                rate *= scaleParameter.getParameterValue(0);
            }

            rate = Math.exp(rate);

        } else if (reciprocal) {

            rate = 1.0 / rate;

        } else {

            if (scaleParameter != null) {

                if (!transformKnown) {
                    setupTransform();
                    transformKnown = true;
                }
                
                rate = logNormalTransform(rate, baseMeasureMu, baseMeasureSigma, transformMu, transformSigma);
            }
            
        }

        if (locationParameter != null) {
            rate *= locationParameter.getParameterValue(0);
        }

        return rate;
    }

    private static double logNormalTransform(double baseDraw,
                                             double baseMeasureMu, double baseMeasureSigma,
                                             double transformMu, double transformSigma) {
        return Math.exp(
                (transformSigma / baseMeasureSigma) * (Math.log(baseDraw) - baseMeasureMu) + transformMu
        );
    }

    private static double getMuPhi(double phi) {
        return -0.5 * Math.log(1.0 + phi);
    }

    private static double getSigmaPhi(double phi) {
        return Math.sqrt(Math.log(1.0 + phi));
    }

    public double getBranchRate(final Tree tree, final NodeRef node) {
        // Branch rates are proportional to time.
        // In the traitLikelihoods, time is proportional to variance
        // Fernandez and Steel (2000) shows the sampling density with the scalar proportional to precision 
        double pre = rates.getNodeValue(tree, node);

        double post = transform(pre);

        if (DEBUG) {
            System.err.println("Pre  : " + pre);
            System.err.println("Post : " + post);
            double average = averageRates(tree);
            double variance = varianceRates(tree, average);
            System.err.println("Avg  : " + average);
            System.err.println("Var  : " + variance);
            System.err.println();
            System.exit(-1);
        }
        
        // else
        return post;
    }

    private double averageRates(final Tree tree) {

        double sum = 0.0;

        for (int i = 0; i < tree.getNodeCount(); ++i) {
            NodeRef node = tree.getNode(i);

            if (!tree.isRoot(node)) {
                sum += transform(rates.getNodeValue(tree, node));
            }
        }

        return sum / (tree.getNodeCount() - 1);
    }

    private double varianceRates(final Tree tree, double average) {

        double sum = 0.0;

        for (int i = 0; i < tree.getNodeCount(); ++i) {
            NodeRef node = tree.getNode(i);

            if (!tree.isRoot(node)) {
                double x = transform(rates.getNodeValue(tree, node));
                sum += (x - average) * (x - average);
            }
        }

        return sum / (tree.getNodeCount() - 1);
    }

    private static final boolean DEBUG = false;

//
//    public int getNodeNumberFromParameterIndex(int parameterIndex) {
//        return rates.getNodeNumberFromParameterIndex(parameterIndex);
//    }

    public int getParameterIndexFromNode(final NodeRef node) {
        return rates.getParameterIndexFromNodeNumber(node.getNumber());
    }

    public Parameter getRateParameter() { return rateParameter; }

    public boolean usingReciprocal() {
        return reciprocal;
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        // Should be called by TreeParameterModel
        fireModelChanged(object, index);
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // Changes to rateParameter are handled by model changed events
        if (variable == locationParameter || variable == scaleParameter) {
            transformKnown = false;
            fireModelChanged();
        } else {
            throw new RuntimeException("Unknown variable");
        }
    }

    protected void storeState() {
    }

    protected void restoreState() {
        transformKnown = false;
    }

    protected void acceptState() {
    }

}

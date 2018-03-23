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
import dr.inference.model.AbstractModel;
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

    private final BranchRateTransform transform;

    public ArbitraryBranchRates(TreeModel tree, Parameter rateParameter, BranchRateTransform transform,
                                boolean setRates) {

        super(ArbitraryBranchRatesParser.ARBITRARY_BRANCH_RATES);

        this.transform = transform;
        if (transform instanceof Model) {
            addModel((Model)transform);
        }

        if (setRates) {
            final double value = transform.center();
            for (int i = 0; i < rateParameter.getDimension(); i++) {
                rateParameter.setValue(i, value);
            }
        }

        //Force the boundaries of rate
        double lower = transform.lower();
        double upper = transform.upper();
        Parameter.DefaultBounds bounds = new Parameter.DefaultBounds(upper, lower, rateParameter.getDimension());
        rateParameter.addBounds(bounds);

        this.rates = new TreeParameterModel(tree, rateParameter, false);
        this.rateParameter = rateParameter;

        addModel(rates);
    }

    public void setBranchRate(Tree tree, NodeRef node, double value) {
        rates.setNodeValue(tree, node, value);
    }

    public double getBranchRateDifferential(final Tree tree, final NodeRef node) {
        return transform.differential(getBranchRate(tree, node));
    }

    public double getBranchRate(final Tree tree, final NodeRef node) {
        // Branch rates are proportional to time.
        // In the traitLikelihoods, time is proportional to variance
        // Fernandez and Steel (2000) shows the sampling density with the scalar proportional to precision

        return transform.transform(rates.getNodeValue(tree, node));
    }

    public int getParameterIndexFromNode(final NodeRef node) {
        return rates.getParameterIndexFromNodeNumber(node.getNumber());
    }

    public Parameter getRateParameter() { return rateParameter; }

    public boolean usingReciprocal() {
        return (transform instanceof BranchRateTransform.Reciprocal);
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == rates) {
            fireModelChanged(object, index);
        } else if (model == transform) {
            fireModelChanged();
        } else {
            throw new RuntimeException("Unknown model");
        }
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) { }

    protected void storeState() { }

    protected void restoreState() { }

    protected void acceptState() { }

    public static BranchRateTransform make(boolean reciprocal, boolean exp) {
        return make(reciprocal, exp, null, null);
    }

    public static BranchRateTransform make(boolean reciprocal, boolean exp,
                                    Parameter location,
                                    Parameter scale) {
        final BranchRateTransform transform;

        if ((reciprocal || exp) && (location != null || scale != null)) {
            throw new RuntimeException("Not yet implemented");
        }

        if (exp) {
            transform = new BranchRateTransform.Exponentiate();
        } else if (reciprocal) {
            transform = new BranchRateTransform.Reciprocal();
        } else {
            if (location != null || scale != null) {
                transform = new BranchRateTransform.LocationScaleLogNormal(
                        ArbitraryBranchRatesParser.ARBITRARY_BRANCH_RATES, location, scale);
            } else {
                transform = new BranchRateTransform.None();
            }
        }
        return transform;
    }

    public interface BranchRateTransform {

        double differential(double raw);  // TODO Correct name?  differentialOfInverse

        double transform(double raw);

        double center();

        double lower();

        double upper();

        abstract class Base implements BranchRateTransform {
            @Override
            public double center() {
                return 1.0;
            }

            @Override
            public double lower() {
                return 0.0;
            }

            @Override
            public double upper() {
                return Double.POSITIVE_INFINITY;
            }
        }

        class None extends Base {

            @Override
            public double differential(double raw) {
                return 1.0 / raw;
            }

            @Override
            public double transform(double raw) {
                return raw;
            }
        }

        class Reciprocal extends Base {

            @Override
            public double differential(double raw) {
                return -raw;
            }

            @Override
            public double transform(double raw) {
                return 1.0 / raw;
            }
        }

        class Exponentiate implements BranchRateTransform {

            @Override
            public double differential(double raw) {
                return 1.0;
            }

            @Override
            public double transform(double raw) {
                return Math.exp(raw);
            }

            @Override
            public double center() {
                return 0.0;
            }

            @Override
            public double lower() {
                return Double.NEGATIVE_INFINITY;
            }

            @Override
            public double upper() {
                return Double.POSITIVE_INFINITY;
            }
        }
        
        class LocationScaleLogNormal extends AbstractModel implements BranchRateTransform  {

            private final Parameter location;
            private final Parameter scale;

            private final double baseMeasureMu = getMuPhi(1.0); // TODO Depends on prior distribution
            private final double baseMeasureSigma = getSigmaPhi(1.0); // TODO Depends on prior distribution

            private double transformMu;
            private double transformSigma;

            private boolean transformKnown;

            LocationScaleLogNormal(String name, Parameter location, Parameter scale) {
                super(name);

                this.location = location;
                this.scale = scale;

                if (location != null) {
                    addVariable(location);
                }

                addVariable(scale);

                this.transformKnown = false;
            }

            @Override
            public double differential(double rate) {
                throw new RuntimeException("Not yet implemented");
            }

            @Override
            public double transform(double rate) {

                if (!transformKnown) {
                    setupTransform();
                    transformKnown = true;
                }

                rate = logNormalTransform(rate, baseMeasureMu, baseMeasureSigma, transformMu, transformSigma);

                if (location != null) {
                    rate *= location.getParameterValue(0);
                }

                return rate;
            }

            @Override
            public double center() {
                return 1;
            }

            @Override
            public double lower() {
                return 0.0;
            }

            @Override
            public double upper() {
                return Double.POSITIVE_INFINITY;
            }

            @Override
            protected void handleModelChangedEvent(Model model, Object object, int index) { }

            @Override
            protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
                if (variable == location || variable == scale) {
                    transformKnown = false;
                    fireModelChanged();
                } else {
                    throw new RuntimeException("Unknown variable");
                }
            }

            @Override
            protected void storeState() { }

            @Override
            protected void restoreState() {
                transformKnown = false;
            }

            @Override
            protected void acceptState() { }

            private void setupTransform() {
                final double phi = scale.getParameterValue(0) * scale.getParameterValue(0);
                transformMu = getMuPhi(phi);
                transformSigma = getSigmaPhi(phi);
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
        }
    }
}

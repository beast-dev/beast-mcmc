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
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.Collections;
import java.util.List;

/**
 * Allows branch rates to take on any double value
 * This is useful for forming a scaled mixture of normals for the continuous diffusion model
 *
 * @author Marc A. Suchard
 * @author Alexei Drummond
 */
public class ArbitraryBranchRates extends AbstractBranchRateModel implements Citable {

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

    public double getBranchRateDifferential(double rate) {
        return transform.differential(rate);
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

    public double getBranchRateSecondDifferential(double rate) {
        return transform.secondDifferential(rate);
    }

    public interface BranchRateTransform {

        double differential(double rate);

        double secondDifferential(double rate);

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
            public double differential(double rate) {
                return 1.0;
            }

            @Override
            public double secondDifferential(double rate) {
                return 0.0;
            }

            @Override
            public double transform(double raw) {
                return raw;
            }
        }

        class Reciprocal extends Base {

            @Override
            public double differential(double rate) {
                return -rate * rate;
            }

            @Override
            public double secondDifferential(double rate) {
                return 2.0 * rate * rate * rate;
            }

            @Override
            public double transform(double raw) {
                return 1.0 / raw;
            }
        }

        class Exponentiate implements BranchRateTransform {

            @Override
            public double differential(double rate) {
                return rate;
            }

            @Override
            public double secondDifferential(double rate) {
                return rate;
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

            private final double baseMeasureMu;
            private final double baseMeasureSigma;

            private double transformMu;
            private double transformSigma;

            private boolean transformKnown;

            LocationScaleLogNormal(String name, Parameter location, Parameter scale) {
                this(name, location, scale, getMuPhi(1.0), getSigmaPhi(1.0));
            }

            LocationScaleLogNormal(String name, Parameter location, Parameter scale,
                                   double baseMeasureMu, double baseMeasureSigma) {

                super(name);

                this.baseMeasureMu = baseMeasureMu;
                this.baseMeasureSigma = baseMeasureSigma;

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

                if (!transformKnown) {
                    setupTransform();
                    transformKnown = true;
                }
                
                double multiplier = (location != null) ? location.getParameterValue(0) : 1.0;
                double raw = logNormalTransform(rate / multiplier,
                        transformMu, transformSigma, baseMeasureMu, baseMeasureSigma);

                return raw > 0. ? (rate * transformSigma) / (raw * baseMeasureSigma) : Double.POSITIVE_INFINITY;
            }

            @Override
            public double secondDifferential(double rate) {

                if (!transformKnown) {
                    setupTransform();
                    transformKnown = true;
                }

                double multiplier = (location != null) ? location.getParameterValue(0) : 1.0;
                double raw = logNormalTransform(rate / multiplier,
                        transformMu, transformSigma, baseMeasureMu, baseMeasureSigma);

                if (rate > 0.0) {
                    return (rate * transformSigma) / (raw * raw * baseMeasureSigma) * (transformSigma / baseMeasureSigma - 1.0);
                }else{
                    if (transformSigma > baseMeasureSigma){
                        return Double.POSITIVE_INFINITY;
                    } else if (transformSigma < baseMeasureSigma) {
                        return Double.NEGATIVE_INFINITY;
                    } else {
                        return 0.0;
                    }
                }
            }


            @Override
            public double transform(double raw) {

                if (!transformKnown) {
                    setupTransform();
                    transformKnown = true;
                }

                double rate = logNormalTransform(raw, baseMeasureMu, baseMeasureSigma, transformMu, transformSigma);

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

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.MOLECULAR_CLOCK;
    }

    @Override
    public String getDescription() {
        return "Location-scale relaxed clock";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CITATION);
    }

    public static Citation CITATION = new Citation(
            new Author[]{
                    new Author("X", "Ji"),
                    new Author("P", "Lemey"),
                    new Author("MA", "Suchard")
            },
            Citation.Status.IN_PREPARATION
    );
}

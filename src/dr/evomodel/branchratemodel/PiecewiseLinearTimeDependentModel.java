/*
 * PiecewiseLinearTimeDependentModel.java
 *
 * Copyright (c) 2002-2022 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.inference.loggers.Loggable;
import dr.inference.model.*;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 */
public class PiecewiseLinearTimeDependentModel extends AbstractModel implements ContinuousBranchValueProvider,
        CountableMixtureBranchRates.TimeDependentModel {

    private final static boolean TEST = true;

    private final TreeModel treeModel;
    private final ParameterPack pack;
    private final Scale scale;

    private boolean slopeInterceptKnown;
    private SlopeInterceptPack slopeInterceptPack;

    public PiecewiseLinearTimeDependentModel(TreeModel treeModel, ParameterPack pack, Scale scale) {
        super("piecewiseLinearBranchValues");
        this.treeModel = treeModel;
        this.pack = pack;
        this.scale = scale;

        addModel(treeModel);

        for (Parameter p : pack) {
            addVariable(p);
        }

        slopeInterceptKnown = false;

        addStatistic(new Statistic.Abstract("intercept") {
            @Override
            public int getDimension() {
                return slopeInterceptPack.intercepts.length;
            }

            @Override
            public double getStatisticValue(int dim) {
                checkSlopeIntercept();
                return slopeInterceptPack.intercepts[dim];
            }
        });

        addStatistic(new Statistic.Abstract("slope") {
            @Override
            public int getDimension() {
                return slopeInterceptPack.slopes.length;
            }

            @Override
            public double getStatisticValue(int dim) {
                checkSlopeIntercept();
                return slopeInterceptPack.slopes[dim];
            }
        });

        addStatistic(new Statistic.Abstract("breaks") {
            @Override
            public int getDimension() {
                return slopeInterceptPack.breaks.length;
            }

            @Override
            public double getStatisticValue(int dim) {
                checkSlopeIntercept();
                return slopeInterceptPack.breaks[dim];
            }
        });
    }

    SlopeInterceptPack getSlopeInterceptPack() {
        checkSlopeIntercept();
        return slopeInterceptPack;
    }

    private double integrate(double x0, double x1, double slope, double intercept) {

        if (TEST) {
            return scale.integral(x0, x1, slope, intercept);
        }

        if (slope == 0) {
            return intercept * (x1 - x0);
        } else {
            if (Double.isInfinite(x0) || Double.isInfinite(x1)) {
                throw new IllegalArgumentException("Unbounded function integral");
            } else {
                double halfSlope = slope / 2;
                return (halfSlope * x1 + intercept) * x1 - (halfSlope * x0 + intercept) * x0;
            }
        }
    }

    double computeIntegratedValue(double parent, double child) {

        final double[] slopes = slopeInterceptPack.slopes;
        final double[] intercepts = slopeInterceptPack.intercepts;
        final double[] breaks = slopeInterceptPack.breaks;

//        if (Double.isInfinite(child)) {
//            return slopes[0];
//        }

        int currentEpoch = 0;
        while (child > scale.inverseTransformTime(breaks[currentEpoch])) {
            ++currentEpoch;
        }

        double integral = 0.0;
        double currentTime = child;
        while (scale.inverseTransformTime(breaks[currentEpoch]) <= parent) {
            integral += integrate(currentTime,
                    scale.inverseTransformTime(breaks[currentEpoch]
                    ),
                    slopes[currentEpoch], intercepts[currentEpoch]);
            currentTime = scale.inverseTransformTime(breaks[currentEpoch]);
            ++currentEpoch;
        }

        integral += integrate(currentTime, parent, slopes[currentEpoch], intercepts[currentEpoch]);

        return integral / (parent - child);
    }

    private void checkSlopeIntercept() {
        if (!slopeInterceptKnown) {
            slopeInterceptPack = pack.getSlopesAndIntercepts();
            slopeInterceptKnown = true;
        }
    }

    @Override
    public double getBranchValue(Tree tree, NodeRef node) {
        
        checkSlopeIntercept();

        double parent = tree.getNodeHeight(tree.getParent(node));
        double child = tree.getNodeHeight(node);

        if (!TEST) {
            parent = scale.transformTime(parent);
            child = scale.transformTime(child);
        }

        double value = computeIntegratedValue(parent, child);

        if (TEST) {
            return Math.log(value);
        } else {
            return scale.inverseTransformRate(value);
        }
    }

    @Override
    public double getMidpointValue(Tree tree, NodeRef node, boolean log) {
        double midpoint = getBranchValue(tree, node);

        if (log) {
            return Math.log(midpoint);
        } else {
            return midpoint;
        }
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model != treeModel) {
            throw new IllegalArgumentException("Unknown model");
        }
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (pack.contains(variable)) {
            slopeInterceptKnown = false;
        } else {
            throw new IllegalArgumentException("Unknown variable");
        }
    }

    @Override
    protected void storeState() {
        // Do nothing
    }

    @Override
    protected void restoreState() {
        slopeInterceptKnown = false;
    }

    @Override
    protected void acceptState() { }

    public enum Scale {
        LOG10_UNIT("log10-rate.unit-time") { // rate = 10^(slope * time + intercept) is the correct function
                                                   //      = exp(log10 * slope * time) * 10^(intercept)
            @Override
            double transformTime(double time) {
                return time;
            }

            @Override
            double inverseTransformTime(double time) {
                return time;
            }

            @Override
            double inverseTransformRate(double rate) {
                return log10 * rate;
            }

            @Override
            double integral(double t0, double t1, double slope, double intercept) { // TODO Check
                if (slope == 0.0) {
                    return Math.exp(log10 * intercept) * (t1 - t0);
                } else {
                    final double factor = log10 * slope;
                    return Math.exp(log10 * intercept) / factor * (Math.exp(factor * t1) - Math.exp(factor * t0));
                }
            }

        },
        LOG10_LOG10("log10-rate.log10-time") { // rate = (time)^(slope) * 10^(intercept) is the correct function
            @Override
            double transformTime(double time) {
                return Math.log10(time);
            }

            @Override
            double inverseTransformTime(double time) {
                return Math.pow(10, time);
            }

            @Override
            double inverseTransformRate(double rate) {
                return log10 * rate;
            }

            @Override
            double integral(double t0, double t1, double slope, double intercept) { // TODO Check
                final double power = slope + 1.0;
                return Math.exp(log10 * intercept) / power * (Math.pow(t1, power) - Math.pow(t0, power));
            }
        };

        Scale(String name) {
            this.name = name;
        }

        abstract double transformTime(double time);

        abstract double inverseTransformTime(double time);

        abstract double inverseTransformRate(double rate);

        abstract double integral(double t0, double t1, double slope, double intercept);

        private final String name;

        public static Scale parse(String text) {
            for (Scale s : Scale.values()) {
                if (s.name.equalsIgnoreCase(text)) {
                    return s;
                }
            }
            return null;
        }
    }

    private final static double log10 = Math.log(10.0);

    abstract public static class ParameterPack implements Iterable<Parameter> {

        final Parameter historicValue;
        final Parameter currentValue;
        final Parameter epochStartTime;

        public ParameterPack(Parameter historicValue,
                             Parameter currentValue,
                             Parameter epochStartTime) {
            this.historicValue = historicValue;
            this.currentValue = currentValue;
            this.epochStartTime = epochStartTime;

            parameterList.add(historicValue);
            parameterList.add(currentValue);
            parameterList.add(epochStartTime);
        }

        final List<Parameter> parameterList = new ArrayList<>();

        public boolean contains(Variable variable) {
            return parameterList.contains((Parameter) variable);
        }

        abstract SlopeInterceptPack getSlopesAndIntercepts();

        @Override
        public Iterator<Parameter> iterator() {
            return parameterList.iterator();
        }
    }

    public static class EpochLengthParameterPack extends ParameterPack {

        final Parameter epochLength;

        public EpochLengthParameterPack(Parameter historicValue,
                                        Parameter currentValue,
                                        Parameter epochStartTime,
                                        Parameter epochLength) {
            super(historicValue, currentValue, epochStartTime);

            this.epochLength = epochLength;
            parameterList.add(epochLength);
        }

        @Override
        SlopeInterceptPack getSlopesAndIntercepts() {
            return new SlopeInterceptPack(this);
        }
    }

    public static class SlopeParameterPack extends ParameterPack {

        final Parameter slope;

        public SlopeParameterPack(Parameter historicValue,
                                  Parameter currentValue,
                                  Parameter epochStartTime,
                                  Parameter slope) {
            super(historicValue, currentValue, epochStartTime);

            this.slope = slope;
            parameterList.add(slope);
        }

        @Override
        SlopeInterceptPack getSlopesAndIntercepts() {
            return new SlopeInterceptPack(this);
        }
    }

    static class SlopeInterceptPack {

        final double[] slopes;
        final double[] intercepts;
        final double[] breaks;

        SlopeInterceptPack(EpochLengthParameterPack pack) {
            double y0 = pack.currentValue.getParameterValue(0);
            double y1 = pack.historicValue.getParameterValue(0);

            double x0 = pack.epochStartTime.getParameterValue(0);
            double x1 = x0 + pack.epochLength.getParameterValue(0);

            double midSlope = (y1 - y0) / (x1 - x0);
            double midIntercept = y0 - midSlope * x0;

            slopes = new double[] { 0.0, midSlope, 0.0 };
            intercepts = new double[] { y0, midIntercept, y1 };
            breaks = new double[] { x0, x1, Double.POSITIVE_INFINITY };

//            for (int i = 0; i < breaks.length - 1; ++i) {
//                breaks[i] = scale.transformTime(breaks[i]);
//            }
        }

        SlopeInterceptPack(SlopeParameterPack pack) {
            double y0 = pack.currentValue.getParameterValue(0);
            double y1 = pack.historicValue.getParameterValue(0);

            double x0 = pack.epochStartTime.getParameterValue(0);
            double midSlope = pack.slope.getParameterValue(0);

            double midIntercept = y0 - midSlope * x0;
            double x1 = (y1 - midIntercept) / midSlope;

            slopes = new double[] { 0.0, midSlope, 0.0 };
            intercepts = new double[] { y0, midIntercept, y1 };
            breaks = new double[] { x0, x1, Double.POSITIVE_INFINITY };

//            for (int i = 0; i < breaks.length - 1; ++i) {
//                breaks[i] = scale.transformTime(breaks[i]);
//            }
        }
    }
}

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
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 */
public class PiecewiseLinearTimeDependentModel extends AbstractModel implements ContinuousBranchValueProvider,
        CountableMixtureBranchRates.TimeDependentModel {

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
    }

    private double integrate(double x0, double x1, double slope, double intercept) {
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

        if (Double.isInfinite(child)) {
            return slopes[0];
        }

        int currentEpoch = 0;
        while (child > breaks[currentEpoch]) {
            ++currentEpoch;
        }

        double integral = 0.0;
        double currentTime = child;
        while (breaks[currentEpoch] <= parent) {
            integral += integrate(currentTime, breaks[currentEpoch],
                    slopes[currentEpoch], intercepts[currentEpoch]);
            currentTime = breaks[currentEpoch];
            ++currentEpoch;
        }

        integral += integrate(currentTime, parent, slopes[currentEpoch], intercepts[currentEpoch]);

        return integral / (parent - child);
    }

    @Override
    public double getBranchValue(Tree tree, NodeRef node) {

        if (!slopeInterceptKnown) {
            slopeInterceptPack = new SlopeInterceptPack(pack);
            slopeInterceptKnown = true;
        }

        double p = tree.getNodeHeight(tree.getParent(node));
        double c = tree.getNodeHeight(node);

        double parent = scale.transformTime(p);
        double child = scale.transformTime(c);

        double value = computeIntegratedValue(parent, child);
        return scale.inverseTransformRate(value);
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
        LOG10_UNIT("log10-rate.unit-time") {
            @Override
            double transformTime(double time) {
                return time;
            }

            @Override
            double inverseTransformRate(double rate) {
                return log10 * rate;
            }
        },
        LOG10_LOG10("log10-rate.log10-time") {
            @Override
            double transformTime(double time) {
                return Math.log10(time);
            }

            @Override
            double inverseTransformRate(double rate) {
                return log10 * rate; // return natural log rate
            }
        };

        Scale(String name) {
            this.name = name;
        }

        abstract double transformTime(double time);

        abstract double inverseTransformRate(double rate);

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

    public static class ParameterPack implements Iterable<Parameter> {

        final Parameter historicValue;
        final Parameter currentValue;
        final Parameter epochStartTime;
        final Parameter epochLength;

        final List<Parameter> parameterList = new ArrayList<>();

        public ParameterPack(Parameter historicValue,
                             Parameter currentValue,
                             Parameter epochStartTime,
                             Parameter epochLength) {
            this.historicValue = historicValue;
            this.currentValue = currentValue;
            this.epochStartTime = epochStartTime;
            this.epochLength = epochLength;

            parameterList.add(historicValue);
            parameterList.add(currentValue);
            parameterList.add(epochStartTime);
            parameterList.add(epochLength);
        }

        public boolean contains(Variable variable) {
            return parameterList.contains((Parameter) variable);
        }

        @Override
        public Iterator<Parameter> iterator() {
            return parameterList.iterator();
        }
    }

    private static class SlopeInterceptPack {

        final double[] slopes;
        final double[] intercepts;
        final double[] breaks;

        SlopeInterceptPack(ParameterPack pack) {
            double y0 = pack.currentValue.getParameterValue(0);
            double y1 = pack.historicValue.getParameterValue(0);

            double x0 = pack.epochStartTime.getParameterValue(0);
            double x1 = x0 + pack.epochLength.getParameterValue(0);

            double midSlope = (y1 - y0) / (x1 - x0);

            slopes = new double[] { 0.0, midSlope, 0.0 };
            intercepts = new double[] { y0, y0 - midSlope * x0, y1 };
            breaks = new double[] { x0, x1, Double.POSITIVE_INFINITY };
        }
    }
}


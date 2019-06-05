/*
 * AutoCorrelatedGradientWrtIncrements.java
 *
 * Copyright (c) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

/**
 * @author Marc A. Suchard
 * @author Xiang Ji
 */
public class AutoCorrelatedGradientWrtIncrements implements GradientWrtParameterProvider, Reportable {

    private final AutoCorrelatedBranchRatesDistribution distribution;
    private final ArbitraryBranchRates branchRates;
    private final Tree tree;

    private final AutoCorrelatedBranchRatesDistribution.BranchVarianceScaling scaling;
    private final AutoCorrelatedBranchRatesDistribution.BranchRateUnits units;

    private Parameter parameter;
    private double[] cachedIncrements;

    public AutoCorrelatedGradientWrtIncrements(AutoCorrelatedBranchRatesDistribution distribution) {
        this.distribution = distribution;
        this.branchRates = distribution.getBranchRateModel();
        this.tree = distribution.getTree();
        this.scaling = distribution.getScaling();
        this.units = distribution.getUnits();
    }

    @Override
    public Likelihood getLikelihood() {
        return distribution.getLikelihood();
    }

    @Override
    public Parameter getParameter() {
        if (parameter == null) {
            parameter = createParameter();
        }
        return parameter;
    }

    @Override
    public int getDimension() {
        return distribution.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {

        double[] gradientWrtIncrements = distribution.getGradientWrtIncrements();
        if (units.needsIncrementCorrection()) {
            recursePostOrderToCorrectGradient(tree.getRoot(), gradientWrtIncrements);
        }

        return gradientWrtIncrements;
    }

    private int recursePostOrderToCorrectGradient(NodeRef node, double[] gradientWrtIncrements) {

        // On STRICTLY_POSITIVE scale, log-likelihood includes log-Jacobian (\sum_{increments} -> rate)

        int numberDescendents = 1;

        if (!tree.isExternal(node)) {
            numberDescendents += recursePostOrderToCorrectGradient(tree.getChild(node, 0), gradientWrtIncrements);
            numberDescendents += recursePostOrderToCorrectGradient(tree.getChild(node, 1), gradientWrtIncrements);
        }

        if (!tree.isRoot(node)) {
            int index = branchRates.getParameterIndexFromNode(node);
            gradientWrtIncrements[index] -= scaling.inverseRescaleIncrement(
                    1.0 * numberDescendents, tree.getBranchLength(node));  // d / d c_i log-Jacobian
        }

        return numberDescendents;
    }

    @Override
    public String getReport() {
        String report;
        try {
            report = new CheckGradientNumerically(
                    this, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, null
            ).getReport();
        } catch (GradientMismatchException e) {
            throw new RuntimeException(e.getMessage());
        }

        return report;
    }

    private Parameter createParameter() {
        return new Parameter.Proxy("increments", distribution.getDimension()) {

            @Override
            public double getParameterValue(int dim) {
                return distribution.getIncrement(dim);
            }

            @Override
            public void setParameterValue(int dim, double value) {
                throw new RuntimeException("Do not set single value at a time");
            }

            @Override
            public void setParameterValueQuietly(int dim, double value) {
                if (cachedIncrements == null) {
                    cachedIncrements = new double[getDimension()];
                }

                cachedIncrements[dim] = value;
            }

            @Override
            public void setParameterValueNotifyChangedAll(int dim, double value) {
                throw new RuntimeException("Do not set single value at a time");
            }

            public void fireParameterChangedEvent(int index, Parameter.ChangeType type) {

                double[] rates = new double[getDimension()];
                recurse(tree.getRoot(), rates, cachedIncrements, 0.0);

                Parameter rateParameter = distribution.getParameter();
                for (int i = 0; i < rates.length; ++i) {
                    rateParameter.setParameterValueQuietly(i, rates[i]);
                }

                rateParameter.fireParameterChangedEvent();
            }
        };
    }

    private void recurse(NodeRef node, double[] rates, double[] increments, double parentIncrement) {

        double increment = parentIncrement;

        if (!tree.isRoot(node)) {
            int index = branchRates.getParameterIndexFromNode(node);
            increment += scaling.inverseRescaleIncrement(increments[index], tree.getBranchLength(node));

            rates[index] = units.inverseTransform(increment);
        }

        if (!tree.isExternal(node)) {
            recurse(tree.getChild(node, 0), rates, increments, increment);
            recurse(tree.getChild(node, 1), rates, increments, increment);
        }
    }
}

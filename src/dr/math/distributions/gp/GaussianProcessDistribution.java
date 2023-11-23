/*
 * GaussianMarkovRandomField.java
 *
 * Copyright (c) 2002-2023 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.math.distributions.gp;

import dr.inference.distribution.RandomField;
import dr.inference.model.*;
import dr.math.distributions.RandomFieldDistribution;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.LinearSolverFactory;
import org.ejml.interfaces.decomposition.CholeskyDecomposition;
import org.ejml.interfaces.linsol.LinearSolver;
import org.ejml.ops.CommonOps;

import java.util.Arrays;

/**
 * @author Marc Suchard
 * @author Filippo Monti
 */
public class GaussianProcessDistribution extends RandomFieldDistribution {

    public static final String TYPE = "GaussianProcess";

    private final int dim;
    private final Parameter meanParameter;
    private final Kernel kernel;
    private final RandomField.WeightProvider weightProvider;

    private final double[] mean;
    private final double[] tmp;
    private final DenseMatrix64F precision;
    private final DenseMatrix64F variance;
    private final LinearSolver<DenseMatrix64F> solver;

    private double logDeterminant;
    
    private boolean meanKnown;
    private boolean precisionAndDeterminantKnown;
    private boolean varianceKnown;

    public GaussianProcessDistribution(String name,
                                       int dim,
                                       Parameter meanParameter,
                                       Kernel kernel,
                                       RandomField.WeightProvider weightProvider) {
        super(name);

        this.dim = dim;
        this.meanParameter = meanParameter;
        this.kernel = kernel;
        this.weightProvider = weightProvider;

        this.mean = new double[dim];
        this.tmp = new double[dim];
        this.precision = new DenseMatrix64F(dim, dim);
        this.variance = new DenseMatrix64F(dim, dim);

        this.solver = LinearSolverFactory.symmPosDef(dim);

        addVariable(meanParameter);

        if (kernel instanceof AbstractModel) {
            addModel((AbstractModel) kernel);
        }

        if (weightProvider != null) {
            addModel(weightProvider);
        }
    }

    private DenseMatrix64F getPrecision() {
        if (!precisionAndDeterminantKnown) {
            DenseMatrix64F variance = getVariance();
            solver.solve(variance, precision);
            CholeskyDecomposition<DenseMatrix64F> d = solver.getDecomposition();
            logDeterminant = Math.log(d.computeDeterminant().getReal());
            precisionAndDeterminantKnown = true;
        }
        return precision;
    }

    private DenseMatrix64F getVariance() {
        if (!varianceKnown) {
            for (int i = 0; i < dim; ++i) {
                for (int j = 0; j < dim; ++j) {
                    variance.set(i, j, kernel.getCorrelation(0, 0)); // TODO
                }
            }
            varianceKnown = true;
        }
        return variance;
    }

    @Override
    public double[] getMean() {
        if (!meanKnown) {
            if (meanParameter == null) {
                Arrays.fill(mean, 0.0);
            } else if (meanParameter.getDimension() == 1) {
                Arrays.fill(mean, meanParameter.getParameterValue(0));
            } else {
                for (int i = 0; i < mean.length; ++i) {
                    mean[i] = meanParameter.getParameterValue(i);
                }
            }
            meanKnown = true;
        }
        return mean;
    }

    @Override
    public double getIncrement(int i, Parameter field) {
        double[] mean = getMean();
        return (field.getParameterValue(i) - mean[i]) - (field.getParameterValue(i + 1) - mean[i + 1]);
    }

    @Override
    public GradientProvider getGradientWrt(Parameter parameter) {
        throw new RuntimeException("Unknown parameter");
    }

    public String getType() {
        return TYPE;
    }

    @Override
    public double[][] getScaleMatrix() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public Variable<Double> getLocationVariable() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double logPdf(double[] x) {

        final double[] mean = getMean();
        final double[] diff = tmp;
        final double[] precision = getPrecision().getData();

        for (int i = 0; i < dim; ++i) {
            diff[i] = x[i] - mean[i];
        }

        double exponent = 0.0;
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                exponent += diff[i] * precision[i * dim + j] * diff[j];
            }
        }

        

        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public int getDimension() { return dim; }

    @Override
    public double[] getGradientLogDensity(Object x) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[] getDiagonalHessianLogDensity(Object x) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[][] getHessianLogDensity(Object x) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public double[] nextRandom() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        throw new IllegalArgumentException("Unknown model");
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    protected void storeState() {

    }

    @Override
    protected void restoreState() {

    }

    @Override
    protected void acceptState() { }
}

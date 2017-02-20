/*
 * FactorOperator.java
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

package dr.inference.operators;

import dr.inference.model.*;
import dr.math.MathUtils;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.SymmetricMatrix;

/**
 * Created with IntelliJ IDEA.
 * User: max
 * Date: 5/22/14
 * Time: 12:49 PM
 * To change this template use File | Settings | File Templates.
 */
public class FactorOperator extends AbstractCoercableOperator {
    private static final String FACTOR_OPERATOR = "factorOperator";
    private LatentFactorModel LFM;
    private MatrixParameter diffusionPrecision;
    double[][] precision;
    double[] mean;
    double[] midMean;
    private int numFactors;
    private boolean randomScan;
    private double scaleFactor;

    public FactorOperator(LatentFactorModel LFM, double weight, boolean randomScan, DiagonalMatrix diffusionPrecision, double scaleFactor, CoercionMode mode) {
        super(mode);
        this.scaleFactor = scaleFactor;
        this.LFM = LFM;
        setWeight(weight);
        this.randomScan = randomScan;
        this.diffusionPrecision = diffusionPrecision;
        setupParameters();

    }

    private void setupParameters() {
        if (numFactors != LFM.getFactorDimension()) {
            numFactors = LFM.getFactorDimension();
            mean = new double[numFactors];
            midMean = new double[numFactors];
            precision = new double[numFactors][numFactors];
        }
    }

    private void getPrecision(double[][] precision) {
        MatrixParameterInterface Loadings = LFM.getLoadings();
        MatrixParameter Precision = LFM.getColumnPrecision();
        int outerDim = Loadings.getRowDimension();
        int innerDim = Loadings.getColumnDimension();
        for (int i = 0; i < outerDim; i++) {
            for (int j = i; j < outerDim; j++) {
                double sum = 0;
                for (int k = j; k < innerDim; k++) {
                    sum += Loadings.getParameterValue(i, k) * Loadings.getParameterValue(j, k) * Precision.getParameterValue(k, k);
                }
                if (i == j) {
                    precision[i][j] = sum + diffusionPrecision.getParameterValue(i, j);
                } else {
                    precision[i][j] = sum;
                    precision[j][i] = sum;
                }
            }

        }
    }

    private void getMean(int column, double[][] variance, double[] midMean, double[] mean) {
        MatrixParameterInterface scaledData = LFM.getScaledData();
        MatrixParameterInterface Precision = LFM.getColumnPrecision();
        MatrixParameterInterface Loadings = LFM.getLoadings();
        for (int i = 0; i < Loadings.getRowDimension(); i++) {
            double sum = 0;
            for (int j = i; j < Loadings.getColumnDimension(); j++) {
                sum += Loadings.getParameterValue(i, j) * Precision.getParameterValue(j, j) * scaledData.getParameterValue(j, column);
            }
            midMean[i] = sum;
        }
        for (int i = 0; i < numFactors; i++) {
            double sum = 0;
            for (int j = 0; j < numFactors; j++) {
                sum += variance[i][j] * midMean[j];
            }
            mean[i] = sum;
        }

//        try {
//            answer=getPrecision().inverse().product(new Matrix(LFM.getLoadings().getParameterAsMatrix())).product(new Matrix(LFM.getColumnPrecision().getParameterAsMatrix())).product(data);
//        } catch (IllegalDimension illegalDimension) {
//            illegalDimension.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
//        }
    }

    private void copy(double[] put, int i) {
        Parameter working = LFM.getFactors().getParameter(i);
        for (int j = 0; j < working.getSize(); j++) {
            working.setParameterValueQuietly(j, put[j]);
        }
        working.fireParameterChangedEvent();
    }

    public int getStepCount() {
        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String getPerformanceSuggestion() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    @Override
    public String getOperatorName() {
        return FACTOR_OPERATOR;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public void randomDraw(int i, double[][] variance) {
        double[] nextValue;
        nextValue = MultivariateNormalDistribution.nextMultivariateNormalVariance(LFM.getFactors().getParameter(i).getParameterValues(), variance, scaleFactor);
//        System.out.println(nextValue[0]);
//        System.out.println(nextValue[1]);
        copy(nextValue, i);
    }

    @Override
    public double doOperation() {
        setupParameters();
        getPrecision(precision);
        double[][] variance = (new SymmetricMatrix(precision)).inverse().toComponents();
        if (randomScan) {
            int i = MathUtils.nextInt(LFM.getFactors().getColumnDimension());
            randomDraw(i, variance);
        }
        for (int i = 0; i < LFM.getFactors().getColumnDimension(); i++) {
            randomDraw(i, variance);
        }

        LFM.getFactors().fireParameterChangedEvent();
//        LFM.computeResiduals();

        return 0;  //To change body of implemented methods use File | Settings | File Templates.
    }


    @Override
    public double getCoercableParameter() {
        return Math.log(scaleFactor);
    }

    @Override
    public void setCoercableParameter(double value) {
        scaleFactor = Math.exp(value);
    }

    @Override
    public double getRawParameter() {
        return scaleFactor;
    }
}

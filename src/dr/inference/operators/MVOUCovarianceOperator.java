/*
 * MVOUCovarianceOperator.java
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

import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.inferencexml.operators.MVOUCovarianceOperatorParser;
import dr.math.distributions.WishartDistribution;
import dr.math.matrixAlgebra.Matrix;

/**
 * @author Marc Suchard
 */
public class MVOUCovarianceOperator extends AbstractCoercableOperator {

    private double mixingFactor;
    private MatrixParameter varMatrix;
    private int dim;

    private MatrixParameter precisionParam;
    private WishartDistribution priorDistribution;
    private int priorDf;
    private double[][] I;
    private Matrix Iinv;


    public MVOUCovarianceOperator(double mixingFactor,
                                  MatrixParameter varMatrix,
                                  int priorDf,
                                  double weight, CoercionMode mode) {
        super(mode);
        this.mixingFactor = mixingFactor;
        this.varMatrix = varMatrix;
        this.priorDf = priorDf;
        setWeight(weight);
        dim = varMatrix.getColumnDimension();
        I = new double[dim][dim];
        for (int i = 0; i < dim; i++)
            I[i][i] = 1.0;
//			I[i][i] = i;
        Iinv = new Matrix(I).inverse();
    }

    public double doOperation() {

        double[][] draw = WishartDistribution.nextWishart(priorDf, I);
//		double[][] good = varMatrix.getParameterAsMatrix();
//		double[][] saveOld = varMatrix.getParameterAsMatrix();

//		System.err.println("draw:\n"+new Matrix(draw));
        double[][] oldValue = varMatrix.getParameterAsMatrix();
        for (int i = 0; i < dim; i++) {
            Parameter column = varMatrix.getParameter(i);
            for (int j = 0; j < dim; j++)
                column.setParameterValue(j,
                        mixingFactor * oldValue[j][i] + (1.0 - mixingFactor) * draw[j][i]
                );

        }
//        varMatrix.fireParameterChangedEvent();
        // calculate Hastings ratio

//		System.err.println("oldValue:\n"+new Matrix(oldValue).toString());
//		System.err.println("newValue:\n"+new Matrix(varMatrix.getParameterAsMatrix()).toString());

        Matrix forwardDrawMatrix = new Matrix(draw);
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
//				saveOld[i][j] *= - mixingFactor;
//				saveOld[i][j] += varMatrix.getParameterValue(i,j);
//				saveOld[i][j] /= 1.0 - mixingFactor;
                oldValue[i][j] -= mixingFactor * varMatrix.getParameterValue(i, j);
                oldValue[i][j] /= 1.0 - mixingFactor;
            }
        }

//		double[][] saveNew = varMatrix.getParameterAsMatrix();

        Matrix backwardDrawMatrix = new Matrix(oldValue);

//		System.err.println("forward:\n"+forwardDrawMatrix);
//		System.err.println("backward:\n"+backwardDrawMatrix);

//		System.err.println("calc start");

//		if( Math.abs(backwardDrawMatrix.component(0,0) + 0.251) < 0.001 ) {
//			System.err.println("found:\n"+backwardDrawMatrix);
//
//			System.err.println("original:\n"+new Matrix(good));
//			System.err.println("draw:\n"+new Matrix(draw));
//			System.err.println("proposed:\n"+new Matrix(varMatrix.getParameterAsMatrix()));
//			System.err.println("mixing = "+mixingFactor);
//			System.err.println("back[0][0] = "+backwardDrawMatrix.component(0,0));
//			System.err.println("saveOld[0][0] = "+saveOld[0][0]);
//
//
//		}

        double bProb = WishartDistribution.logPdf(backwardDrawMatrix, Iinv, priorDf, dim,
//				WishartDistribution.computeNormalizationConstant(Iinv,priorDf,dim));
                0);

        if (bProb == Double.NEGATIVE_INFINITY) {
//            throw new OperatorFailedException("Not reversible");
            // not clear if this means a HR of -Inf or a RuntimeException
            return Double.NEGATIVE_INFINITY;
        }

        double fProb = WishartDistribution.logPdf(forwardDrawMatrix, Iinv, priorDf, dim,
//				WishartDistribution.computeNormalizationConstant(Iinv,priorDf,dim));
                0);

//		System.err.println("calc end");

//		if( fProb == Double.NEGATIVE_INFINITY ) {
//			System.err.println("forwards is problem");
//			System.exit(-1);
//		}

//		if( bProb == Double.NEGATIVE_INFINITY ) {
//			System.err.println("backwards is problem");
//			System.exit(-1);
//		}

//		System.err.println("fProb = "+fProb);
//		System.err.println("bProb = "+bProb);

//		System.exit(-1);

        return bProb - fProb;
    }

    //MCMCOperator INTERFACE
    public final String getOperatorName() {
        return MVOUCovarianceOperatorParser.MVOU_OPERATOR + "(" +
                varMatrix.getId() + ")";
    }

    public double getCoercableParameter() {
        return Math.log(mixingFactor / (1.0 - mixingFactor));
//		return Math.log((1.0 - mixingFactor) / mixingFactor);
    }

    public void setCoercableParameter(double value) {
        mixingFactor = Math.exp(value) / (1.0 + Math.exp(value));
//		mixingFactor = Math.exp(-value) / (1.0 + Math.exp(-value));
    }

    public double getRawParameter() {
        return mixingFactor;
    }

    public double getMixingFactor() {
        return mixingFactor;
    }

    public double getTargetAcceptanceProbability() {
        return 0.234;
    }

    public double getMinimumAcceptanceLevel() {
        return 0.1;
    }

    public double getMaximumAcceptanceLevel() {
        return 0.4;
    }

    public double getMinimumGoodAcceptanceLevel() {
        return 0.20;
    }

    public double getMaximumGoodAcceptanceLevel() {
        return 0.30;
    }

    public final String getPerformanceSuggestion() {

        double prob = MCMCOperator.Utils.getAcceptanceProbability(this);
        double targetProb = getTargetAcceptanceProbability();
        dr.util.NumberFormatter formatter = new dr.util.NumberFormatter(5);
        double sf = OperatorUtils.optimizeWindowSize(mixingFactor, prob, targetProb);
        if (prob < getMinimumGoodAcceptanceLevel()) {
            return "Try setting mixingFactor to about " + formatter.format(sf);
        } else if (prob > getMaximumGoodAcceptanceLevel()) {
            return "Try setting mixingFactor to about " + formatter.format(sf);
        } else return "";
    }

}

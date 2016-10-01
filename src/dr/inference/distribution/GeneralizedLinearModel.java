/*
 * GeneralizedLinearModel.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inference.distribution;

import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.SingularValueDecomposition;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.*;
import dr.inferencexml.distribution.GeneralizedLinearModelParser;
import dr.math.MultivariateFunction;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * @author Marc Suchard
 */
@Deprecated // GLM stuff is now in inference.glm - this is here for backwards compatibility temporarily
public abstract class GeneralizedLinearModel extends AbstractModelLikelihood implements MultivariateFunction {

    protected Parameter dependentParam;
    protected List<Parameter> independentParam;
    protected List<Parameter> indParamDelta;
    protected List<DesignMatrix> designMatrix; // fixed constants, access as double[][] to save overhead

    //    protected double[][] scaleDesignMatrix;
    protected int[] scaleDesign;
    protected Parameter scaleParameter;

    protected int numIndependentVariables = 0;
    protected int numRandomEffects = 0;
    protected int N;

    protected List<Parameter> randomEffects = null;

    public GeneralizedLinearModel(Parameter dependentParam) {
        super(GeneralizedLinearModelParser.GLM_LIKELIHOOD);
        this.dependentParam = dependentParam;

        if (dependentParam != null) {
            addVariable(dependentParam);
            N = dependentParam.getDimension();
        } else
            N = 0;
    }

    //    public double[][] getScaleDesignMatrix() { return scaleDesignMatrix; }
    public int[] getScaleDesign() {
        return scaleDesign;
    }

    public void addRandomEffectsParameter(Parameter effect) {
        if (randomEffects == null) {
            randomEffects = new ArrayList<Parameter>();
        }
        if (N != 0 && effect.getDimension() != N) {
            throw new RuntimeException("Random effects have the wrong dimension");
        }
        addVariable(effect);
        randomEffects.add(effect);
        numRandomEffects++;
    }

    public void addIndependentParameter(Parameter effect, DesignMatrix matrix, Parameter delta) {
        if (designMatrix == null)
            designMatrix = new ArrayList<DesignMatrix>();
        if (independentParam == null)
            independentParam = new ArrayList<Parameter>();
        if (indParamDelta == null)
            indParamDelta = new ArrayList<Parameter>();

        if (N == 0) {
            N = matrix.getRowDimension();
        }
        designMatrix.add(matrix);
        independentParam.add(effect);
        indParamDelta.add(delta);

        if (designMatrix.size() != independentParam.size())
            throw new RuntimeException("Independent variables and their design matrices are out of sync");
        addVariable(effect);
        addVariable(matrix);
        if (delta != null)
            addVariable(delta);
        numIndependentVariables++;
        Logger.getLogger("dr.inference").info("\tAdding independent predictors '" + effect.getStatisticName() + "' with design matrix '" + matrix.getStatisticName() + "'");
    }

    public boolean getAllIndependentVariablesIdentifiable() {

        int totalColDim = 0;
        for (DesignMatrix mat : designMatrix) {
            totalColDim += mat.getColumnDimension();
        }

        double[][] grandDesignMatrix = new double[N][totalColDim];

        int offset = 0;
        for (DesignMatrix mat : designMatrix) {
            final int length = mat.getColumnDimension();
            for (int i = 0; i < N; ++i) {
                for (int j = 0; j < length; ++j) {
                    grandDesignMatrix[i][offset + j] = mat.getParameterValue(i, j);
                }
            }
            offset += length;
        }

        double[][] mat = grandDesignMatrix;

        if (grandDesignMatrix.length < grandDesignMatrix[0].length) {
            mat = new double[grandDesignMatrix[0].length][grandDesignMatrix.length];

            for (int i = 0; i < grandDesignMatrix.length; ++i) {
                for (int j = 0; j < grandDesignMatrix[i].length; ++j) {
                    mat[j][i] = grandDesignMatrix[i][j];
                }
            }
        }

        SingularValueDecomposition svd = new SingularValueDecomposition(
                new DenseDoubleMatrix2D(mat));

        int rank = svd.rank();
        boolean isFullRank = (totalColDim == rank);
        Logger.getLogger("dr.inference").info("\tTotal # of predictors = " + totalColDim + " and rank = " + rank);
        return isFullRank;
    }

    public int getNumberOfFixedEffects() {
        return numIndependentVariables;
    }

    public int getNumberOfRandomEffects() {
        return numRandomEffects;
    }

    public double[] getXBeta() {

        double[] xBeta = new double[N];

        for (int j = 0; j < numIndependentVariables; j++) {
            Parameter beta = independentParam.get(j);
            Parameter delta = indParamDelta.get(j);
            DesignMatrix X = designMatrix.get(j);
            final int K = beta.getDimension();
            for (int k = 0; k < K; k++) {
                double betaK = beta.getParameterValue(k);
                if (delta != null)
                    betaK *= delta.getParameterValue(k);
                for (int i = 0; i < N; i++)
                    xBeta[i] += X.getParameterValue(i, k) * betaK;
            }
        }

        for (int j = 0; j < numRandomEffects; j++) {
            Parameter effect = randomEffects.get(j);
            for (int i = 0; i < N; i++) {
                xBeta[i] += effect.getParameterValue(i);
            }
        }

        return xBeta;
    }

    public Parameter getFixedEffect(int j) {
        return independentParam.get(j);
    }

    public Parameter getRandomEffect(int j) {
        return randomEffects.get(j);
    }

    public Parameter getDependentVariable() {
        return dependentParam;
    }

    public double[] getXBeta(int j) {

        double[] xBeta = new double[N];

        Parameter beta = independentParam.get(j);
        Parameter delta = indParamDelta.get(j);
        DesignMatrix X = designMatrix.get(j);
        final int K = beta.getDimension();
        for (int k = 0; k < K; k++) {
            double betaK = beta.getParameterValue(k);
            if (delta != null)
                betaK *= delta.getParameterValue(k);
            for (int i = 0; i < N; i++)
                xBeta[i] += X.getParameterValue(i, k) * betaK;
        }

        if (numRandomEffects != 0)
            throw new RuntimeException("Attempting to retrieve fixed effects without controlling for random effects");

        return xBeta;

    }

    public int getEffectNumber(Parameter effect) {
        return independentParam.indexOf(effect);
    }

//	public double[][] getXtScaleX(int j) {
//
//		final Parameter beta = independentParam.get(j);
//		double[][] X = designMatrix.get(j);
//		final int dim = X[0].length;
//
//		if( dim != beta.getDimension() )
//			throw new RuntimeException("should have checked eariler");
//
//		double[] scale = getScale();
//
//
//	}

    public double[][] getX(int j) {
        return designMatrix.get(j).getParameterAsMatrix();
    }


    public double[] getScale() {

        double[] scale = new double[N];

//        final int K = scaleParameter.getDimension();
//        for (int k = 0; k < K; k++) {
//            final double scaleK = scaleParameter.getParameterValue(k);
//            for (int i = 0; i < N; i++)
//                scale[i] += scaleDesignMatrix[i][k] * scaleK;
//        }
        for (int k = 0; k < N; k++)
            scale[k] = scaleParameter.getParameterValue(scaleDesign[k]);

        return scale;
    }


    public double[][] getScaleAsMatrix() {

//        double[][] scale = new double[N][N];
//
//        return scale;
        throw new RuntimeException("Not yet implemented: GeneralizedLinearModel.getScaleAsMatrix()");
    }

//	protected abstract double calculateLogLikelihoodAndGradient(double[] beta, double[] gradient);

    protected abstract double calculateLogLikelihood(double[] beta);

    protected abstract double calculateLogLikelihood();

    protected abstract boolean confirmIndependentParameters();

    public abstract boolean requiresScale();

    public void addScaleParameter(Parameter scaleParameter, Parameter design) {
        this.scaleParameter = scaleParameter;
//        this.scaleDesignMatrix = matrix.getParameterAsMatrix();
        scaleDesign = new int[design.getDimension()];
        for (int i = 0; i < scaleDesign.length; i++)
            scaleDesign[i] = (int) design.getParameterValue(i);
        addVariable(scaleParameter);
    }

/*	// **************************************************************
          // RealFunctionOfSeveralVariablesWithGradient IMPLEMENTATION
	// **************************************************************


	public double eval(double[] beta, double[] gradient) {
		return calculateLogLikelihoodAndGradient(beta, gradient);
	}


	public double eval(double[] beta) {
		return calculateLogLikelihood(beta);
	}


	public int getNumberOfVariables() {
		return independentParam.getDimension();
	}*/

    // ************
    //       Mutlivariate implementation
    // ************


    public double evaluate(double[] beta) {
        return calculateLogLikelihood(beta);
    }

    public int getNumArguments() {
        int total = 0;
        for (Parameter effect : independentParam)
            total += effect.getDimension();
        return total;
    }

    public double getLowerBound(int n) {
        int which = n;
        int k = 0;
        while (which > independentParam.get(k).getDimension()) {
            which -= independentParam.get(k).getDimension();
            k++;
        }
        return independentParam.get(k).getBounds().getLowerLimit(which);
    }

    public double getUpperBound(int n) {
        int which = n;
        int k = 0;
        while (which > independentParam.get(k).getDimension()) {
            which -= independentParam.get(k).getDimension();
            k++;
        }
        return independentParam.get(k).getBounds().getUpperLimit(which);
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
//        fireModelChanged();
    }

    protected void storeState() {
        // No internal states to save
    }

    protected void restoreState() {
        // No internal states to restore
    }

    protected void acceptState() {
        // Nothing to do
    }

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        return calculateLogLikelihood();
    }

    @Override
    public String toString() {
        return super.toString() + ": " + getLogLikelihood();
    }

    public void makeDirty() {
    }

    // **************************************************************
    // Loggable IMPLEMENTATION
    // **************************************************************

//    /**
//     * @return the log columns.
//     */
//    public LogColumn[] getColumns() {
//        return new dr.inference.loggers.LogColumn[]{
//                new LikelihoodColumn(getId())
//        };
//    }
//
//    private class LikelihoodColumn extends dr.inference.loggers.NumberColumn {
//        public LikelihoodColumn(String label) {
//            super(label);
//        }
//
//        public double getDoubleValue() {
//            return getLogLikelihood();
//        }
//    }

    public LogColumn[] getColumns() {
        LogColumn[] output = new LogColumn[N];
        for (int i = 0; i < N; i++)
            output[i] = new NumberArrayColumn(getId() + i, i);
        return output;
    }

    private class NumberArrayColumn extends NumberColumn {

        private final int index;

        public NumberArrayColumn(String label, int index) {
            super(label);
            this.index = index;
        }

        public double getDoubleValue() {
            return getXBeta()[index];
        }
    }

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

    public Element createElement(Document d) {
        throw new RuntimeException("Not implemented yet!");
    }

}

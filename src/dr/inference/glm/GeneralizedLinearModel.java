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

package dr.inference.glm;

import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.SingularValueDecomposition;
import dr.inference.distribution.DensityModel;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.distribution.ParametricMultivariateDistributionModel;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.*;
import dr.inferencexml.glm.GeneralizedLinearModelParser;
import dr.math.MultivariateFunction;
import dr.math.distributions.Distribution;
import dr.util.Transform;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * @author Marc Suchard
 * @author Andrew Rambaut
 */
public final class GeneralizedLinearModel extends AbstractModelLikelihood {


    public enum LinkFunction {
        IDENTITY(new Transform.NoTransform()),
        LOG(new Transform.LogTransform()),
        LOGIT(new Transform.LogitTransform());

        LinkFunction(Transform transform) {
            this.transform = transform;
        }

        public Transform getTransform() {
            return transform;
        }

        private final Transform transform;
    }

    private final Transform linkFunction;
    private final DensityModel density;
    private final boolean isMultivariateDensity;

    private final Parameter dependentParameter;
    private final List<Parameter> independentParameter = new ArrayList<Parameter>();
    private final List<Parameter> independentParameterDelta = new ArrayList<Parameter>();
    private final List<DesignMatrix> designMatrix = new ArrayList<DesignMatrix>();

    private int numIndependentVariables = 0;
    private int numRandomEffects = 0;
    private int N;

    protected List<Parameter> randomEffects = null;

    private double[] transformedXBeta;
    private double[] storedTransformedXBeta;
    private boolean transformedXBetaKnown = false;

    private double[] Y;

    private double storedLogLikelihood;
    private double logLikelihood;
    private boolean likelihoodKnown = false;

    public GeneralizedLinearModel(Parameter dependentParameter, DensityModel density, LinkFunction linkFunction) {
        super(GeneralizedLinearModelParser.GLM_LIKELIHOOD);
        this.dependentParameter = dependentParameter;
        this.linkFunction = linkFunction.getTransform();
        this.density = density;
        isMultivariateDensity = density instanceof ParametricMultivariateDistributionModel;

        addModel(density);

        if (dependentParameter != null) {
            addVariable(dependentParameter);
            N = dependentParameter.getDimension();
        } else {
            N = 0;
        }

        transformedXBeta = new double[N];
        storedTransformedXBeta = new double[N];

        Y = dependentParameter.getParameterValues();

        transformedXBetaKnown = false;
        likelihoodKnown = false;
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
        if (N == 0) {
            N = matrix.getRowDimension();
        }

        designMatrix.add(matrix);
        independentParameter.add(effect);
        independentParameterDelta.add(delta);

        if (designMatrix.size() != independentParameter.size()) {
            throw new RuntimeException("Independent variables and their design matrices are out of sync");
        }

        addVariable(effect);
        addVariable(matrix);

        if (delta != null) {
            addVariable(delta);
        }

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

        SingularValueDecomposition svd = new SingularValueDecomposition(new DenseDoubleMatrix2D(mat));

        int rank = svd.rank();
        boolean isFullRank = (totalColDim == rank);
        Logger.getLogger("dr.inference").info("\tTotal # of predictors = " + totalColDim + " and rank = " + rank);
        return isFullRank;
    }

    public double[] getXBeta() {

        double[] xBeta = new double[N];

        for (int j = 0; j < numIndependentVariables; j++) {
            Parameter beta = independentParameter.get(j);
            Parameter delta = independentParameterDelta.get(j);
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

    public double[] getXBeta(int j) {

        double[] xBeta = new double[N];

        Parameter beta = independentParameter.get(j);
        Parameter delta = independentParameterDelta.get(j);
        DesignMatrix X = designMatrix.get(j);
        final int K = beta.getDimension();
        for (int k = 0; k < K; k++) {
            double betaK = beta.getParameterValue(k);
            if (delta != null) {
                betaK *= delta.getParameterValue(k);
            }
            for (int i = 0; i < N; i++) {
                xBeta[i] += X.getParameterValue(i, k) * betaK;
            }
        }

        if (numRandomEffects != 0) {
            throw new RuntimeException("Attempting to retrieve fixed effects without controlling for random effects");
        }

        return xBeta;
    }

    public int getNumberOfFixedEffects() {
        return numIndependentVariables;
    }

    public int getNumberOfRandomEffects() {
        return numRandomEffects;
    }

    public Parameter getFixedEffect(int j) {
        return independentParameter.get(j);
    }

    public Parameter getRandomEffect(int j) {
        return randomEffects.get(j);
    }

    public Parameter getDependentVariable() {
        return dependentParameter;
    }

    public int getEffectNumber(Parameter effect) {
        return independentParameter.indexOf(effect);
    }

    public double[][] getX(int j) {
        return designMatrix.get(j).getParameterAsMatrix();
    }

    private void calculateTransformedXBeta() {
        double[] xBeta = getXBeta();

        for (int i = 0; i < N; i++) {
            transformedXBeta[i] = linkFunction.inverse(xBeta[i]);
        }

        transformedXBetaKnown = true;
    }

    private double calculateLogLikelihood() {
        if (!transformedXBetaKnown) {
            calculateTransformedXBeta();
        }

        double logL = 0.0;

        if (isMultivariateDensity) {
            // todo - implement
        } else {
            for (int i = 0; i < Y.length; i++) {
                density.getLocationVariable().setValue(0, transformedXBeta[i]);
                logL += density.logPdf(new double[] { Y[i] });
            }
        }

        return logL;
    }

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    /**
     * todo - cache likelihood values
     * @return
     */
    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
        }
        return logLikelihood;
    }

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // some aspect of the density function has changed
        likelihoodKnown = false;
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == dependentParameter) {
            Y = dependentParameter.getParameterValues();
        }

        transformedXBetaKnown = false;
        likelihoodKnown = false;
    }

    @Override
    protected void storeState() {
        storedLogLikelihood = logLikelihood;
        System.arraycopy(transformedXBeta, 0, storedTransformedXBeta, 0, transformedXBeta.length);
    }

    @Override
    protected void restoreState() {
        logLikelihood = storedLogLikelihood;

        // could use double buffering to speed this up for very large matrices...
        System.arraycopy(storedTransformedXBeta, 0, transformedXBeta, 0, transformedXBeta.length);
    }

    @Override
    protected void acceptState() {
        // Nothing to do
    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public String toString() {
        return super.toString() + ": " + getLogLikelihood();
    }

    @Override
    public void makeDirty() {
    }

    // **************************************************************
    // Loggable IMPLEMENTATION
    // **************************************************************

    // probably makes more consistent sense to log the likelihood and expose the xBeta values
    // through a statistic...

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

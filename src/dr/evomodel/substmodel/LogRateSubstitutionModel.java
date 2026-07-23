/*
 * LogRateSubstitutionModel.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.evomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.inference.loggers.LogColumn;
import dr.inference.model.BayesianStochasticSearchVariableSelection;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.WrappedMatrix;
import dr.util.Citation;
import dr.util.CommonCitations;
import dr.util.Transform;

import java.util.*;

/**
 * @author Filippo Monti
 * @author Marc A. Suchard
 */

public class LogRateSubstitutionModel extends ComplexSubstitutionModel implements DifferentiableSubstitutionModel { // implements ParameterReplaceableSubstitutionModel

    double[] Q = new double[stateCount * stateCount];
    double normalizingConstant;
    boolean Qknown = false;

    public LogRateSubstitutionModel(String name, DataType dataType, FrequencyModel rootFreqModel,
                                    LogAdditiveCtmcRateProvider lrm) {

        super(name, dataType, rootFreqModel, null);
        this.lrm = lrm;
        addModel(lrm);
        testProbabilities = new double[stateCount*stateCount];

    }

    @Override
    public LogAdditiveCtmcRateProvider getRateProvider() {
        return lrm;
    }

    protected void setupRelativeRates(double[] rates) {
        double[] transformedRates = lrm.getRates();
        System.arraycopy(transformedRates,0,rates,0,rates.length);
    }

    @Override
    public Set<Likelihood> getLikelihoodSet() {
        return new HashSet<>(Collections.singletonList(this));
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == lrm) {
            updateMatrix = true;
            Qknown = false;
            fireModelChanged();
        } else {
            super.handleModelChangedEvent(model, object, index);
        }
    }

    public LogColumn[] getColumns() {
        //Aggregate columns from ComplexSubstitutionModel with glm.columns
        LogColumn[] aggregated = new LogColumn[lrm.getColumns().length + 2];
        int index = 0;
        for (LogColumn col : lrm.getColumns()) {
            aggregated[index] = col;
            index++;
        }
        aggregated[index++] = new LikelihoodColumn(getId() + ".L");
        aggregated[index] = new NormalizationColumn(getId() + ".Norm");

        return aggregated;
        //return lrm.getColumns();
    }

    public double getLogLikelihood() {
        double logL = super.getLogLikelihood();
        if (logL == 0 &&
                BayesianStochasticSearchVariableSelection.Utils.connectedAndWellConditioned(testProbabilities,this)) { // Also check that graph is connected
            return 0;
        }
        return Double.NEGATIVE_INFINITY;
    }

    @Override
    public String getDescription() { return "Log Rate Substitution Model"; }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CommonCitations.LEMEY_2014_UNIFYING); // TODO change citation
    }

    private final LogAdditiveCtmcRateProvider lrm;
    private final double[] testProbabilities;

    public Transform getTransform() {
        return lrm.getTransform();
    }

//    @Override
//    public WrappedMatrix getInfinitesimalDifferentialMatrix(DifferentialMassProvider.DifferentialWrapper.WrtParameter wrt) {
//        return DifferentiableSubstitutionModelUtil.getInfinitesimalDifferentialMatrix(wrt, this);
//    }
    private void cachingQMatrix() {
        if (!Qknown) {
            normalizingConstant = this.setupMatrix();
            this.getInfinitesimalMatrix(Q);
            Qknown = true;
        }
    }
    @Override
    public WrappedMatrix getInfinitesimalDifferentialMatrix(DifferentialMassProvider.DifferentialWrapper.WrtParameter wrt) { // TODO this is duplicated code from DifferentialMassProvider

        cachingQMatrix();

        final double[] differentialRates = new double[rateCount];
        ((DifferentiableSubstitutionModel) this).setupDifferentialRates(wrt, differentialRates, normalizingConstant);

        final double[] differentialFrequencies = new double[stateCount];
        ((DifferentiableSubstitutionModel) this).setupDifferentialFrequency(wrt, differentialFrequencies);

        double[][] differentialMassMatrix = new double[stateCount][stateCount];
        DifferentiableSubstitutionModelUtil.setupQDerivative(this, differentialRates, differentialFrequencies, differentialMassMatrix);
        this.makeValid(differentialMassMatrix, stateCount);

        final double weightedNormalizationGradient
                = ((DifferentiableSubstitutionModel) this).getWeightedNormalizationGradient(
                wrt, differentialMassMatrix, differentialFrequencies);

        for (int i = 0; i < stateCount; i++) {
            for (int j = 0; j < stateCount; j++) {
                differentialMassMatrix[i][j] -= Q[i * stateCount + j] * weightedNormalizationGradient;
            }
        }

        WrappedMatrix differential = new WrappedMatrix.ArrayOfArray(differentialMassMatrix);

        return differential;
//        return DifferentiableSubstitutionModelUtil.getInfinitesimalDifferentialMatrix(wrt, this);
    }

    public class WrtLogRate implements DifferentialMassProvider.DifferentialWrapper.WrtParameter {

        private final int stateCount;
        private final int dim;
        /**
         * @param stateCount The number of states in the substitution model.
         */
        public WrtLogRate(int stateCount, int dim) {
            this.stateCount = stateCount;
            this.dim = dim;
        }

        @Override
        public double getRate(int switchCase) {
            throw new RuntimeException("Should not be called.");
        }

        @Override
        public double getNormalizationDifferential() {
            return 0; // TODO update this
        }

        @Override
        public void setupDifferentialFrequencies(double[] differentialFrequencies, double[] frequencies) {
            Arrays.fill(differentialFrequencies, 1);
        }

        @Override
        public void setupDifferentialRates(double[] differentialRates, double[] Q, double normalizingConstant) {
            // Initialize all derivatives to zero.
            Arrays.fill(differentialRates, 0.0);

            if (dim < differentialRates.length / 2) {
                int k = 0;
                for (int i = 0; i < stateCount; ++i) {
                    for (int j = i + 1; j < stateCount; ++j) {
                        if (k == dim) {
                            differentialRates[k] = Q[index(i, j)];
                        }
                        k++;
                    }
                }
            } else {
                int k = differentialRates.length / 2;
                for (int j = 0; j < stateCount; ++j) {
                    for (int i = j + 1; i < stateCount; ++i) {
                        if (k == dim) {
                            differentialRates[k] = Q[index(i, j)];
                        }
                        k++;
                    }
                }
            }
        }
        private int index(int i, int j) { return i * stateCount + j; }
    }

    @Override
    public DifferentialMassProvider.DifferentialWrapper.WrtParameter factory(Parameter parameter, int dim) {
        return new LogRateSubstitutionModel.WrtLogRate(stateCount, dim);
    }


    @Override
    public void setupDifferentialRates(DifferentialMassProvider.DifferentialWrapper.WrtParameter wrt,
                                       double[] differentialRates, double normalizingConstant) {
        cachingQMatrix();
//        final double[] Q = new double[stateCount * stateCount];
//        getInfinitesimalMatrix(Q); // TODO These are large; should cache
        wrt.setupDifferentialRates(differentialRates, Q, normalizingConstant);
    }

    @Override
    public void setupDifferentialFrequency(DifferentialMassProvider.DifferentialWrapper.WrtParameter wrt, double[] differentialFrequency) {
        wrt.setupDifferentialFrequencies(differentialFrequency, getFrequencyModel().getFrequencies());
    }

    @Override
    public double getWeightedNormalizationGradient(DifferentialMassProvider.DifferentialWrapper.WrtParameter wrt, double[][] differentialMassMatrix, double[] differentialFrequencies) {
        double derivative = 0;

        if (getNormalization()) {
            for (int i = 0; i < stateCount; ++i) {
                derivative -= differentialMassMatrix[i][i] * getFrequencyModel().getFrequency(i); // TODO CHECK THIS
            }
        }
        return derivative;
    }

}

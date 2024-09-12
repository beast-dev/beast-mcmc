/*
 * GlmSubstitutionModel.java
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
import dr.inference.distribution.GeneralizedLinearModel;
import dr.inference.distribution.LogLinearModel;
import dr.inference.loggers.LogColumn;
import dr.inference.model.BayesianStochasticSearchVariableSelection;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.WrappedMatrix;
import dr.util.Citation;
import dr.util.CommonCitations;

import java.util.*;

/**
 * @author Marc A. Suchard
 */

public class GlmSubstitutionModel extends ComplexSubstitutionModel
        implements ParameterReplaceableSubstitutionModel, DifferentiableSubstitutionModel {

    public GlmSubstitutionModel(String name, DataType dataType, FrequencyModel rootFreqModel,
                                LogAdditiveCtmcRateProvider glm) {

        super(name, dataType, rootFreqModel, null);
        this.glm = glm;
        addModel(glm);
        testProbabilities = new double[stateCount*stateCount];

    }

    @Override
    public LogAdditiveCtmcRateProvider getRateProvider() {
        return glm;
    }

    public GeneralizedLinearModel getGeneralizedLinearModel() {
        if (glm instanceof GeneralizedLinearModel) {
            return (GeneralizedLinearModel) glm;
        }
        throw new RuntimeException("Not yet implemented");
    }

    protected void setupRelativeRates(double[] rates) {
        System.arraycopy(glm.getXBeta(),0,rates,0,rates.length);
    }

    @Override
    public Set<Likelihood> getLikelihoodSet() {
        return new HashSet<>(Arrays.asList(this, glm));
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == glm) {
            updateMatrix = true;
            fireModelChanged();
        } else {
            super.handleModelChangedEvent(model, object, index);
        }
    }

    public LogColumn[] getColumns() {
        //Aggregate columns from ComplexSubstitutionModel with glm.columns
        LogColumn[] aggregated = new LogColumn[glm.getColumns().length + 2];
        int index = 0;
        for (LogColumn col : glm.getColumns()) {
            aggregated[index] = col;
            index++;
        }
        aggregated[index++] = new LikelihoodColumn(getId() + ".L");
        aggregated[index] = new NormalizationColumn(getId() + ".Norm");

        return aggregated;
        //return glm.getColumns();
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
    public String getDescription() {
        return "Generalized linear (model, GLM) substitution model"; // TODO Horrible; fix
    }

    @Override
    public List<Citation> getCitations() {

        return Collections.singletonList(CommonCitations.LEMEY_2014_UNIFYING);
    }

    private final LogAdditiveCtmcRateProvider glm;
    private final double[] testProbabilities;

    @Override
    public ParameterReplaceableSubstitutionModel factory(List<Parameter> oldParameters, List<Parameter> newParameters) {

        LogLinearModel newGLM = ((LogLinearModel)glm).factory(oldParameters, newParameters);
        return new GlmSubstitutionModel(getModelName(), dataType, freqModel, newGLM);
    }

    @Override
    public WrappedMatrix getInfinitesimalDifferentialMatrix(DifferentialMassProvider.DifferentialWrapper.WrtParameter wrt) {
        // TODO all instantiations of this function currently do the same thing; remove duplication
        return DifferentiableSubstitutionModelUtil.getInfinitesimalDifferentialMatrix(wrt, this);
    }

    static class WrtOldGLMSubstitutionModelParameter implements DifferentialMassProvider.DifferentialWrapper.WrtParameter {

        final private int dim;
        final private int fixedEffectIndex;
        final private int stateCount;
        final private LogLinearModel glm;

        public WrtOldGLMSubstitutionModelParameter(LogLinearModel glm, int fixedEffectIndex, int dim, int stateCount) {
            this.glm = glm;
            this.fixedEffectIndex = fixedEffectIndex;
            this.dim = dim;
            this.stateCount = stateCount;
        }

        @Override
        public double getRate(int switchCase) {
            throw new RuntimeException("Should not be called.");
        }

        @Override
        public double getNormalizationDifferential() {
            return 0;
        }

        @Override
        public void setupDifferentialFrequencies(double[] differentialFrequencies, double[] frequencies) {
            Arrays.fill(differentialFrequencies, 1);
        }

        @Override
        public void setupDifferentialRates(double[] differentialRates, double[] Q, double normalizingConstant) {

            final double[] covariate = glm.getDesignMatrix(fixedEffectIndex).getColumnValues(dim);

            int k = 0;
            for (int i = 0; i < stateCount; ++i) {
                for (int j = i + 1; j < stateCount; ++j) {

                    differentialRates[k] = covariate[k] * Q[index(i, j)];
                    k++;

                }
            }

            for (int j = 0; j < stateCount; ++j) {
                for (int i = j + 1; i < stateCount; ++i) {

                    differentialRates[k] = covariate[k] * Q[index(i, j)];
                    k++;

                }
            }

//            final double chainRule = getChainRule();
////            double[][] design = glm.getX(effect);
//
//            for (int i = 0; i < differentialRates.length; ++i) {
//                differentialRates[i] = covariate[i] / normalizingConstant * chainRule;
//            }
        }

        double getChainRule() {
            return Math.exp(glm.getFixedEffect(fixedEffectIndex).getParameterValue(dim));
        }

        private int index(int i, int j) {
            return i * stateCount + j;
        }
    }

    @Override
    public DifferentialMassProvider.DifferentialWrapper.WrtParameter factory(Parameter parameter, int dim) {

        final int effectIndex = ((LogLinearModel)glm).getEffectNumber(parameter);
        if (effectIndex == -1) {
            throw new RuntimeException("Only implemented for single dimensions, break up beta to one for each block for now please.");
        }
        return new WrtOldGLMSubstitutionModelParameter((LogLinearModel) glm, effectIndex, dim, stateCount);
    }

//    @Override
//    public DifferentialMassProvider.DifferentialWrapper.WrtParameter factory(Parameter parameter, int dim) {
//        for (int i = 0; i < glm.getNumberOfFixedEffects(); ++i) {
//            Parameter effect = glm.getFixedEffect(i);
//            if (parameter == effect) {
//                return new WrtGlmCoefficient(effect, i, dim);
//            }
//        }
//        throw new RuntimeException("Parameter not found");
//    }

    @Override
    public void setupDifferentialRates(DifferentialMassProvider.DifferentialWrapper.WrtParameter wrt, double[] differentialRates, double normalizingConstant) {
        final double[] Q = new double[stateCount * stateCount];
        getInfinitesimalMatrix(Q); // TODO These are large; should cache
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
                derivative -= differentialMassMatrix[i][i] * getFrequencyModel().getFrequency(i);
            }
        }
        return derivative;
    }
}

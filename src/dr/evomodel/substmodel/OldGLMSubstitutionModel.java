/*
 * OldGLMSubstitutionModel.java
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

package dr.evomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.inference.distribution.GeneralizedLinearModel;
import dr.inference.distribution.LogLinearModel;
import dr.inference.loggers.LogColumn;
import dr.inference.model.BayesianStochasticSearchVariableSelection;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.WrappedMatrix;
import dr.util.Citation;
import dr.util.CommonCitations;

import java.util.Collections;
import java.util.List;

/**
 * @author Marc A. Suchard
 */
@Deprecated
public class OldGLMSubstitutionModel extends ComplexSubstitutionModel implements  DifferentiableSubstitutionModel {

    public OldGLMSubstitutionModel(String name, DataType dataType, FrequencyModel rootFreqModel,
                                   LogLinearModel glm) {

        super(name, dataType, rootFreqModel, null);
        this.glm = glm;
        addModel(glm);
        testProbabilities = new double[stateCount*stateCount];

    }

    public GeneralizedLinearModel getGeneralizedLinearModel() { return glm; }

    protected void setupRelativeRates(double[] rates) {
        System.arraycopy(glm.getXBeta(),0,rates,0,rates.length);
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == glm) {
            updateMatrix = true;
            fireModelChanged();
        }
        else
            super.handleModelChangedEvent(model,object,index);
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
        aggregated[index++] = new NormalizationColumn(getId() + ".Norm");

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

    final private LogLinearModel glm;
    final private double[] testProbabilities;

    @Override
    public WrappedMatrix getInfinitesimalDifferentialMatrix(DifferentialMassProvider.DifferentialWrapper.WrtParameter wrt) {
        // TODO all instantiations of this function currently do the same thing; remove duplication
        return DifferentiableSubstitutionModelUtil.getInfinitesimalDifferentialMatrix(wrt, this);
    }

    @Override
    public DifferentialMassProvider.DifferentialWrapper.WrtParameter factory(Parameter parameter, int dim) {
        for (int i = 0; i < glm.getNumberOfFixedEffects(); ++i) {
            Parameter effect = glm.getFixedEffect(i);
            if (parameter == effect) {
                return new WrtGlmCoefficient(effect, i, dim);
            }
        }
        throw new RuntimeException("Parameter not found");
    }

    @Override
    public void setupDifferentialRates(DifferentialMassProvider.DifferentialWrapper.WrtParameter wrt,
                                       double[] differentialRates,
                                       double normalizingConstant) {
        double[] relativeRates = new double[rateCount];
        setupRelativeRates(relativeRates); // TODO These are large; should cache
        wrt.setupDifferentialRates(differentialRates, relativeRates, normalizingConstant);
    }

    @Override
    public void setupDifferentialFrequency(DifferentialMassProvider.DifferentialWrapper.WrtParameter wrt,
                                           double[] differentialFrequency) {
        double[] frequencies = freqModel.getFrequencies();
        System.arraycopy(frequencies, 0, differentialFrequency, 0, frequencies.length);
    }

    @Override
    public double getWeightedNormalizationGradient(DifferentialMassProvider.DifferentialWrapper.WrtParameter wrt,
                                                   double[][] differentialMassMatrix,
                                                   double[] differentialFrequencies) {
        double weight = 0.0;
        for (int i = 0; i < stateCount; ++i) {
            weight -= differentialMassMatrix[i][i] * getFrequencyModel().getFrequency(i);
        }
        return weight;
    }

    class WrtGlmCoefficient implements DifferentialMassProvider.DifferentialWrapper.WrtParameter {

        final private Parameter parameter;
        final int effect;
        final int dim;

        public WrtGlmCoefficient(Parameter parameter, int effect, int dim) {
            this.parameter = parameter;
            this.effect = effect;
            this.dim = dim;
        }

        @Override
        public void setupDifferentialRates(double[] differentialRates,
                                           double[] relativeRates,
                                           double normalizingConstant) {

            final double chainRule = getChainRule();
            double[][] design = glm.getX(effect);

            for (int i = 0; i < relativeRates.length; ++i) {
                differentialRates[i] = design[i][dim] / normalizingConstant * chainRule;
            }
        }

        double getChainRule() {
            return Math.exp(parameter.getParameterValue(dim));
        }

        @Override
        public double getRate(int switchCase) {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        public double getNormalizationDifferential() {
            throw new RuntimeException("Not yet implemented");
        }

        @Override
        public void setupDifferentialFrequencies(double[] differentialFrequencies, double[] frequencies) {
            throw new RuntimeException("Not yet implemented");
        }
    }
}

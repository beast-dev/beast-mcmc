/*
 * MG94HKYCodonModel.java
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

package dr.evomodel.substmodel.codon;

import dr.evomodel.substmodel.*;
import dr.evolution.datatype.Codons;
import dr.inference.model.Parameter;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import dr.math.matrixAlgebra.WrappedMatrix;

/**
 * Muse-Gaut model of codon evolution
 *
 * @author Marc A. Suchard
 * @author Guy Baele
 * @author Philippe lemey
 */
public class MG94HKYCodonModel extends MG94CodonModel implements ParameterReplaceableSubstitutionModel {

    protected Parameter kappaParameter;
//    protected Parameter AtoC_Parameter;
//    protected Parameter AtoG_Parameter;
//    protected Parameter AtoT_Parameter;
//    protected Parameter CtoG_Parameter;
//    protected Parameter CtoT_Parameter;
//    protected Parameter GtoT_Parameter;

    public MG94HKYCodonModel(Codons codonDataType, Parameter alphaParameter, Parameter betaParameter, Parameter kappaParameter,
                             FrequencyModel freqModel) {
        this(codonDataType, alphaParameter, betaParameter, kappaParameter, freqModel,
                new DefaultEigenSystem(codonDataType.getStateCount()));
    }

    public MG94HKYCodonModel(Codons codonDataType,
                             Parameter alphaParameter,
                             Parameter betaParameter,
                             Parameter kappaParameter,
//                          Parameter AtoC_Parameter,
//                          Parameter AtoG_Parameter,
//                          Parameter AtoT_Parameter,
//                          Parameter CtoG_Parameter,
//                          Parameter CtoT_Parameter,
//                          Parameter GtoT_Parameter,
                             FrequencyModel freqModel, EigenSystem eigenSystem) {
        super(codonDataType, alphaParameter, betaParameter, freqModel, eigenSystem);

        this.kappaParameter = kappaParameter;
        addVariable(kappaParameter);
        kappaParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0,
        kappaParameter.getDimension()));

    }

    public double getKappa() {
        return kappaParameter.getParameterValue(0);
    }

    protected void setupRelativeRates(double[] rates) {

        double alpha = getAlpha() / numSynTransitions;
        double beta = getBeta() / numNonsynTransitions;

        double kappa = getKappa();
        for (int i = 0; i < rateCount; i++) {
            switch (rateMap[i]) {
                case 0:
                    rates[i] = 0.0;
                    break;            // codon changes in more than one codon position
                case 1:
                    rates[i] = alpha  * kappa;
                    break;        // synonymous transition
                case 2:
                    rates[i] = alpha;
                    break;        // synonymous transversion
                case 3:
                    rates[i] = beta  * kappa;
                    break;         // non-synonymous transition
                case 4:
                    rates[i] = beta;
                    break;            // non-synonymous transversion
            }
        }
    }

    @Override
    public SubstitutionModel replaceParameter(Parameter oldParameter, Parameter newParameter) {
        if (oldParameter == alphaParameter) {
            return new MG94HKYCodonModel(codonDataType, newParameter, betaParameter, kappaParameter, freqModel);
        } else if (oldParameter == betaParameter) {
            return new MG94HKYCodonModel(codonDataType, alphaParameter, newParameter, kappaParameter, freqModel);
        } else if (oldParameter == kappaParameter) {
            return new MG94HKYCodonModel(codonDataType, alphaParameter, betaParameter, newParameter, freqModel);
        } else {
            throw new RuntimeException("Not yet implemented!");
        }
    }

    @Override
    public double[] getDifferentialMassMatrix(double time, Parameter parameter) {
        WrappedMatrix.ArrayOfArray infinitesimalDifferentialMatrix = getInfinitesimalDifferentialMatrix(parameter);
        double[] result = DifferentiableSubstitutionModelUtil.getDifferentialMassMatrix(time, stateCount, infinitesimalDifferentialMatrix, eigenDecomposition);
        return result;
    }

    protected WrappedMatrix.ArrayOfArray getInfinitesimalDifferentialMatrix(Parameter parameter) {
        if (parameter == alphaParameter || parameter == betaParameter) {
            final double alphaPlusBetaInverse = 1.0 / (getAlpha() + getBeta());
            final double normalizingConstant = setupMatrix();
            final double[] Q = new double[stateCount * stateCount];
            getInfinitesimalMatrix(Q);
            final double[] differentialRates = new double[rateCount];
            setupDifferentialRates(parameter, differentialRates, normalizingConstant);


            WrappedMatrix.ArrayOfArray differentialMassMatrix = new WrappedMatrix.ArrayOfArray(new double[stateCount][stateCount]);
            setupQMatrix(differentialRates, freqModel.getFrequencies(), differentialMassMatrix.getArrays());
            makeValid(differentialMassMatrix.getArrays(), stateCount);
            final double weightedNormalizationGradient
                    = super.getNormalizationValue(differentialMassMatrix.getArrays(), freqModel.getFrequencies()) - alphaPlusBetaInverse;

            for (int i = 0; i < stateCount; i++) {
                for (int j = 0; j < stateCount; j++) {
                    final double result = differentialMassMatrix.get(i, j) - Q[i * stateCount + j] * weightedNormalizationGradient;
                    differentialMassMatrix.set(i, j, result);
                }
            }

            return differentialMassMatrix;
        } else {
            throw new RuntimeException("Not yet implemented");
        }
    }

    private void setupDifferentialRates(Parameter parameter, double[] differentialRates, double normalizingConstant) {
        if (parameter == alphaParameter) {
            final double kappa = getKappa();
            for (int i = 0; i < rateCount; i++) {
                switch (rateMap[i]) {
                    case 0:
                        differentialRates[i] = 0.0;
                        break;
                    case 1:
                        differentialRates[i] = kappa / normalizingConstant / numSynTransitions;
                        break;        // synonymous transition
                    case 2:
                        differentialRates[i] = 1.0 / normalizingConstant / numSynTransitions;
                        break;        // synonymous transversion
                    case 3:
                        differentialRates[i] = 0.0;
                        break;
                    case 4:
                        differentialRates[i] = 0.0;
                        break;
                }
            }
        } else if (parameter == betaParameter) {
            final double kappa = getKappa();
            for (int i = 0; i < rateCount; i++) {
                switch (rateMap[i]) {
                    case 0:
                        differentialRates[i] = 0.0;
                        break;
                    case 1:
                        differentialRates[i] = 0.0;
                        break;
                    case 2:
                        differentialRates[i] = 0.0;
                        break;
                    case 3:
                        differentialRates[i] = kappa / normalizingConstant / numNonsynTransitions;
                        break;         // non-synonymous transition
                    case 4:
                        differentialRates[i] = 1.0 / normalizingConstant / numNonsynTransitions;
                        break;            // non-synonymous transversion
                }
            }
        } else {
            throw new RuntimeException("Not yet implemented!");
        }
    }
}
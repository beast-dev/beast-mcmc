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
import dr.evomodel.substmodel.DifferentialMassProvider.DifferentialWrapper.WrtParameter;
import dr.evolution.datatype.Codons;
import dr.inference.model.Parameter;

import java.util.List;

/**
 * Muse-Gaut model of codon evolution
 *
 * @author Marc A. Suchard
 * @author Guy Baele
 * @author Philippe lemey
 */
public class MG94HKYCodonModel extends MG94CodonModel {

    protected Parameter kappaParameter;

    public MG94HKYCodonModel(Codons codonDataType, Parameter alphaParameter, Parameter betaParameter,
                             Parameter kappaParameter, FrequencyModel freqModel) {
        this(codonDataType, alphaParameter, betaParameter, kappaParameter, freqModel,
                new DefaultEigenSystem(codonDataType.getStateCount()));
    }

    private MG94HKYCodonModel(Codons codonDataType,
                             Parameter alphaParameter,
                             Parameter betaParameter,
                             Parameter kappaParameter,
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
    public ParameterReplaceableSubstitutionModel factory(List<Parameter> oldParameters, List<Parameter> newParameters) {

        Parameter alpha = alphaParameter;
        Parameter beta = betaParameter;
        Parameter kappa = kappaParameter;
        FrequencyModel frequencyModel = freqModel;

        assert(oldParameters.size() == newParameters.size());

        for (int i = 0; i < oldParameters.size(); i++) {

            Parameter oldParameter = oldParameters.get(i);
            Parameter newParameter = newParameters.get(i);

            if (oldParameter == alphaParameter) {
                alpha = newParameter;
            } else if (oldParameter == betaParameter) {
                beta = newParameter;
            } else if (oldParameter == kappaParameter) {
                kappa = newParameter;
            } else {
                throw new RuntimeException("Not yet implemented!");
            }
        }
        return new MG94HKYCodonModel(codonDataType, alpha, beta, kappa, frequencyModel);
    }

    @Override
    public WrtParameter factory(Parameter parameter) {
        WrtMG94HKYModelParameter wrt;
        if (parameter == alphaParameter) {
            wrt = WrtMG94HKYModelParameter.ALPHA;
        } else if (parameter == betaParameter) {
            wrt = WrtMG94HKYModelParameter.BETA;
        } else {
            throw new RuntimeException("Not yet implemented!");
        }
        return wrt;
    }

    enum WrtMG94HKYModelParameter implements WrtParameter {
        ALPHA {
            @Override
            public double getRate(int switchCase, double normalizingConstant,
                                  DifferentiableSubstitutionModel substitutionModel) {
                MG94HKYCodonModel thisSubstitutionModel = (MG94HKYCodonModel) substitutionModel;
                final double numSynTransitions = thisSubstitutionModel.getNumSynTransitions();
                final double kappa = thisSubstitutionModel.getKappa();
                switch (switchCase) {
                    case 0: return 0.0;
                    case 1: return kappa / normalizingConstant / numSynTransitions; // synonymous transition
                    case 2: return 1.0 / normalizingConstant / numSynTransitions; // synonymous transversion
                    case 3: return 0.0;
                    case 4: return 0.0;
                }
                throw new IllegalArgumentException("Invalid switch case");
            }
        },
        BETA {
            @Override
            public double getRate(int switchCase, double normalizingConstant, DifferentiableSubstitutionModel substitutionModel) {
                MG94HKYCodonModel thisSubstitutionModel = (MG94HKYCodonModel) substitutionModel;
                final double numNonsynTransitions = thisSubstitutionModel.getNumNonsynTransitions();
                final double kappa = thisSubstitutionModel.getKappa();
                switch (switchCase) {
                    case 0: return 0.0;
                    case 1: return 0.0;
                    case 2: return 0.0;
                    case 3: return kappa / normalizingConstant / numNonsynTransitions; // non-synonymous transversion
                    case 4: return 1.0 / normalizingConstant / numNonsynTransitions; // non-synonymous transversion
                }
                throw new IllegalArgumentException("Invalid switch case");
            }
        }
    }
}
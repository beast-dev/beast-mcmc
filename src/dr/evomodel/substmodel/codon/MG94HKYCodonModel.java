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

import dr.evolution.datatype.Codons;
import dr.evomodel.substmodel.*;
import dr.evomodel.substmodel.DifferentialMassProvider.DifferentialWrapper.WrtParameter;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.WrappedMatrix;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.Collections;
import java.util.List;

/**
 * Muse-Gaut model of codon evolution
 *
 * @author Marc A. Suchard
 * @author Guy Baele
 * @author Philippe lemey
 */
public class MG94HKYCodonModel extends AbstractCodonModel implements Citable,
        ParameterReplaceableSubstitutionModel, DifferentiableSubstitutionModel {

    protected Parameter alphaParameter;
    protected Parameter betaParameter;
    protected Parameter kappaParameter;

    private final int numSynTransitions;
    private final int numNonsynTransitions;
    private final CodonOptions options;

    public MG94HKYCodonModel(Codons codonDataType,
                             Parameter alphaParameter,
                             Parameter betaParameter,
                             Parameter kappaParameter,
                             FrequencyModel freqModel,
                             CodonOptions options) {
        this(codonDataType, alphaParameter, betaParameter, kappaParameter, freqModel, options,
                new DefaultEigenSystem(codonDataType.getStateCount()));
    }

    MG94HKYCodonModel(Codons codonDataType,
                      Parameter alphaParameter,
                      Parameter betaParameter,
                      Parameter kappaParameter,
                      FrequencyModel freqModel,
                      CodonOptions options,
                      EigenSystem eigenSystem) {
        super("MG94", codonDataType, freqModel, eigenSystem);

        this.alphaParameter = alphaParameter;
        addVariable(alphaParameter);
        alphaParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0,
                alphaParameter.getDimension()));

        this.betaParameter = betaParameter;
        addVariable(betaParameter);
        betaParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0,
                betaParameter.getDimension()));

        this.kappaParameter = kappaParameter;
        addVariable(kappaParameter);
        kappaParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0,
                kappaParameter.getDimension()));

        numSynTransitions = getNumSynTransitions();
        numNonsynTransitions = getNumNonsynTransitions();

        this.options = options;
    }

    private int countRates(int i, int j) {
        int count = 0;

        for (byte rate : rateMap) {
            if (rate == i || rate == j) {
                count++;
            }
        }
        return count;
    }

    private int getNumSynTransitions() {
        return 2 * countRates(1, 2);
    }

    private int getNumNonsynTransitions() {
        return 2 * countRates(3, 4);
    }

    private double getNormalizationRatioForParameterization() {
        return getTotalS() + getTotalN();
        //(numSynTransitions * getAlpha() + numNonsynTransitions * getBeta()) / (numSynTransitions + numNonsynTransitions);
    }

    protected double getNormalizationValue(double[][] matrix, double[] pi) {
        double norm = 1.0;
        if (doNormalization) {
            double ratio = getNormalizationRatioForParameterization();
            norm = super.getNormalizationValue(matrix, pi) / ratio;
        }
        return norm;
    }

    private double getDS() {
        double rate = alphaParameter.getParameterValue(0);
        if (options.isParameterTotalRate) {
            rate /= numSynTransitions;
        }
        return rate;
    }

    private double getDN() {
        double rate = betaParameter.getParameterValue(0);
        if (options.isParameterTotalRate) {
            rate /= numNonsynTransitions;
        }
        return rate;
    }

    private double getTotalS() {
        double rate = alphaParameter.getParameterValue(0);
        if (!options.isParameterTotalRate) {
            rate *= numSynTransitions;
        }
        return rate;
    }

    private double getTotalN() {
        double rate = betaParameter.getParameterValue(0);
        if (!options.isParameterTotalRate) {
            rate *= numNonsynTransitions;
        }
        return rate;
    }

    public double getKappa() {
        return kappaParameter.getParameterValue(0);
    }

    protected void setupRelativeRates(double[] rates) {

        double alpha = getDS();
        double beta = getDN();
        double kappa = getKappa();

        for (int i = 0; i < rateCount; i++) {
            switch (rateMap[i]) {
                case 0:
                    rates[i] = 0.0;
                    break;            // codon changes in more than one codon position
                case 1:
                    rates[i] = alpha  * kappa;
                    break;            // synonymous transition
                case 2:
                    rates[i] = alpha;
                    break;            // synonymous transversion
                case 3:
                    rates[i] = beta  * kappa;
                    break;            // non-synonymous transition
                case 4:
                    rates[i] = beta;
                    break;            // non-synonymous transversion
            }
        }
    }

    public void setNormalization(boolean normalize) {
        this.doNormalization = normalize;
    }

    private boolean doNormalization = true;

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.SUBSTITUTION_MODELS;
    }

    @Override
    public String getDescription() {
        return "Muse-Gaut codon substitution model";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CITATION);
    }

    public static Citation CITATION = new Citation(
            new Author[]{
                    new Author("SV", "Muse"),
                    new Author("BS", "Gaut")
            },
            "A likelihood approach for comparing synonymous and non-synonymous nucleotide substitution rates, with application to the chloroplast genome",
            1994,
            "Molecular Biology and Evolution",
            11, 715, 724
    );

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
                throw new RuntimeException("Unknown parameter");
            }
        }
        return new MG94HKYCodonModel(codonDataType, alpha, beta, kappa, frequencyModel, options);
    }

    public void setupDifferentialRates(WrtParameter wrt, double[] differentialRates, double normalizingConstant) {

        for (int i = 0; i < rateCount; ++i) {
            differentialRates[i] = wrt.getRate(rateMap[i]) / normalizingConstant;
        }
    }

    @Override
    public double getWeightedNormalizationGradient(WrtParameter wrtParameter, double[][] differentialMassMatrix, double[] frequencies) {
        double normalizationRatio = getNormalizationRatioForParameterization();
        double normalizationRatioInverse = normalizationRatio == 1.0 ?
                0.0 :
                1.0 / getNormalizationRatioForParameterization() * wrtParameter.getNormalizationDifferential();
        return getNormalizationValue(differentialMassMatrix, frequencies) - normalizationRatioInverse;
        // TODO This is not correct when doNormalization == false
    }

    @Override
    public WrappedMatrix getInfinitesimalDifferentialMatrix(WrtParameter wrt) {
        return DifferentiableSubstitutionModelUtil.getInfinitesimalDifferentialMatrix(wrt, this);
    }

    @Override
    public WrtParameter factory(Parameter parameter) {

        if (parameter == alphaParameter) {
            return new Alpha(getNumSynTransitions(), options.isParameterTotalRate);
        } else if (parameter == betaParameter) {
            return new Beta(getNumNonsynTransitions(), options.isParameterTotalRate);
        } else {
            throw new RuntimeException("Not yet implemented");
        }
    }

    abstract class WrtMG94ModelParameter implements WrtParameter {

        private final double normalizationDifferential;
        final double perEventRateScalar;

        WrtMG94ModelParameter(int eventCount, boolean isParameterTotalRate) {
            if (isParameterTotalRate) {
                normalizationDifferential = 1.0;
                perEventRateScalar = 1.0 / eventCount;
            } else {
                normalizationDifferential = eventCount;
                perEventRateScalar = 1.0;
            }
        }

        @Override
        public double getNormalizationDifferential() {
            return normalizationDifferential;
        }
    }

    class Alpha extends WrtMG94ModelParameter {

        Alpha(int eventCount, boolean isParameterTotalRate) {
            super(eventCount, isParameterTotalRate);
        }

        @Override
        public double getRate(int switchCase) {

            final double kappa = getKappa();
            switch (switchCase) {
                case 0:
                    return 0.0;
                case 1:
                    return kappa * perEventRateScalar; // synonymous transition
                case 2:
                    return 1.0 * perEventRateScalar;   // synonymous transversion
                case 3:
                    return 0.0;
                case 4:
                    return 0.0;
            }
            throw new IllegalArgumentException("Invalid switch case");
        }
    }

    class Beta extends WrtMG94ModelParameter {

        Beta(int eventCount, boolean isParameterTotalRate) {
            super(eventCount, isParameterTotalRate);
        }

        @Override
        public double getRate(int switchCase) {

            final double kappa = getKappa();
            switch (switchCase) {
                case 0:
                    return 0.0;
                case 1:
                    return 0.0;
                case 2:
                    return 0.0;
                case 3:
                    return kappa * perEventRateScalar; // non-synonymous transversion
                case 4:
                    return 1.0 * perEventRateScalar;   // non-synonymous transversion
            }
            throw new IllegalArgumentException("Invalid switch case");
        }
    }
}
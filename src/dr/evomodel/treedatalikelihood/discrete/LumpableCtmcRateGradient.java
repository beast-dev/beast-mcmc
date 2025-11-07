/*
 * LogCtmcRateGradient.java
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

package dr.evomodel.treedatalikelihood.discrete;

import dr.evomodel.substmodel.*;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.loggers.LogColumn;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.util.Citation;
import dr.util.CommonCitations;
import dr.util.Transform;

import java.util.Collections;
import java.util.List;

/**
 * @author Xinghua Tao
 * @author Marc A. Suchard
 */

public class LumpableCtmcRateGradient extends AbstractLogAdditiveSubstitutionModelGradient {

    private final StronglyLumpableCtmcRates rateProvider;
    private final CompoundParameter parameter;
    private final int[][] mapEffectToIndices;

    private final ParameterDimensionLink[] link;

    public LumpableCtmcRateGradient(String traitName,
                                    TreeDataLikelihood treeDataLikelihood,
                                    BeagleDataLikelihoodDelegate likelihoodDelegate,
                                    ComplexSubstitutionModel substitutionModel,
                                    CompoundParameter compoundParameter) {
        super(traitName, treeDataLikelihood, likelihoodDelegate, substitutionModel,
                ApproximationMode.FIRST_ORDER);
        this.rateProvider = extractRateProvider(substitutionModel);
        this.parameter = compoundParameter;

        this.mapEffectToIndices = makeAsymmetricMap();
        this.link = createLink(compoundParameter);

    }

    private ParameterDimensionLink[] createLink(CompoundParameter compoundParameter) {
        ParameterDimensionLink[] link = new ParameterDimensionLink[compoundParameter.getDimension()];
        int index = 0;
        for (int i = 0; i < compoundParameter.getParameterCount(); ++i) {
            Parameter parameter = compoundParameter.getParameter(i);
            for (int j = 0; j < parameter.getDimension(); ++j) {
                link[index] = new ParameterDimensionLink(parameter, j);
                ++index;
            }
        }
        return link;
    }

    private StronglyLumpableCtmcRates extractRateProvider(ComplexSubstitutionModel substitutionModel) {
        if (substitutionModel.getRateProvider() instanceof StronglyLumpableCtmcRates) {
            return (StronglyLumpableCtmcRates) substitutionModel.getRateProvider();
        } else {
            throw new IllegalArgumentException("Invalid substitution model");
        }
    }

    @Override
    protected double preProcessNormalization(double[] differentials, double[] generator,
                                             boolean normalize) {
        double total = 0.0;
        if (normalize) {
            for (int i = 0; i < stateCount; ++i) {
                for (int j = 0; j < stateCount; ++j) {
                    final int ij = i * stateCount + j;
                    total += differentials[ij] * generator[ij];
                }
            }
        }
        return total;
    }
    
    private int[][] makeAsymmetricMap() {
        int[][] map = new int[stateCount * (stateCount - 1)][];

        int k = 0;
        for (int i = 0; i < stateCount; ++i) {
            for (int j = i + 1; j < stateCount; ++j) {
                map[k++] = new int[]{i, j};
            }
        }

        for (int j = 0; j < stateCount; ++j) {
            for (int i = j + 1; i < stateCount; ++i) {
                map[k++] = new int[]{i, j};
            }
        }

        return map;
    }

    @Override
    double processSingleGradientDimension(int k, double[] differentials, double[] generator, double[] pi,
                                          boolean normalize, double normalizationGradientContribution,
                                          double normalizationScalar,
                                          Transform transform, boolean scaleByFrequencies) {

        final int i = mapEffectToIndices[k][0], j = mapEffectToIndices[k][1];
        final int ii = i * stateCount + i;
        final int ij = i * stateCount + j;

        double element;
        if (transform == null) {
            element = generator[ij]; // Default is exp()
        } else {
            final Parameter transformedParameter = rateProvider.getLogRateParameter();
            element = transform.gradient(transformedParameter.getParameterValue(k));
            if (normalize) {
                element *= normalizationScalar;
            }
            if (scaleByFrequencies) {
                element *= pi[i];
            }
        }

        double total = (differentials[ij]  - differentials[ii]) * element;

        if (normalize) {
            total -= element * pi[i] * normalizationGradientContribution;
        }

        return total;
    }

    @Override
    public Parameter getParameter() {
        return parameter;
    }

    @Override
    public LogColumn[] getColumns() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.SUBSTITUTION_MODELS;
    }

    @Override
    public String getDescription() {
        return null; // TODO
    }

    @Override
    public List<Citation> getCitations() {
        // TODO Update
        return Collections.singletonList(CommonCitations.MONTI_GENERIC_RATES_2024);
    }

    static class ParameterDimensionLink {

        private final Parameter parameter;
        private final int index;

        ParameterDimensionLink(Parameter parameter, int index) {
            this.parameter = parameter;
            this.index = index;
        }

        public double getDifferential(double element, StronglyLumpableCtmcRates.SuperInfo info) {
            return element / parameter.getParameterValue(index);
//            return info.getRate() / parameter.getParameterValue(index);
        }
    }
}

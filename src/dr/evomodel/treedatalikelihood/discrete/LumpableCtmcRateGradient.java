/*
 * LumpableCtmcRateGradient.java
 *
 * Copyright © 2002-2025 the BEAST Development Team
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

    private final CompoundParameter parameter;
    private final ParameterDimensionLink[] links;

    public LumpableCtmcRateGradient(String traitName,
                                    TreeDataLikelihood treeDataLikelihood,
                                    BeagleDataLikelihoodDelegate likelihoodDelegate,
                                    ComplexSubstitutionModel substitutionModel,
                                    CompoundParameter compoundParameter) {
        super(traitName, treeDataLikelihood, likelihoodDelegate, substitutionModel,
                ApproximationMode.FIRST_ORDER);
        this.parameter = compoundParameter;
        this.links = createLinks(compoundParameter,
                (StronglyLumpableCtmcRates) substitutionModel.getRateProvider());
    }

    private static ParameterDimensionLink[] createLinks(CompoundParameter compoundParameter,
                                                 StronglyLumpableCtmcRates rateProvider) {
        ParameterDimensionLink[] link = new ParameterDimensionLink[compoundParameter.getDimension()];
        int index = 0;
        for (int i = 0; i < compoundParameter.getParameterCount(); ++i) {
            Parameter parameter = compoundParameter.getParameter(i);
            for (int j = 0; j < parameter.getDimension(); ++j) {
                link[index] = new ParameterDimensionLink(parameter, j,
                        rateProvider.searchForParameterAndDimension(parameter, j));
                ++index;
            }
        }
        return link;
    }

    @Override
    double processSingleGradientDimension(int k, double[] differentials, double[] generator, double[] pi,
                                          boolean normalize, double normalizationGradientContribution,
                                          double normalizationScalar,
                                          Transform transform, boolean scaleByFrequencies) {

        ParameterDimensionLink lk = links[k];

        if (lk.entries.size() == 0) {
            throw new RuntimeException("Should not get here");
        }

        double total = 0;
        for (int[] index : lk.entries) {

            final int i = index[0], j = index[1];
            final int ii = i * stateCount + i;
            final int ij = i * stateCount + j;

            double chain = lk.getDifferential(generator[ij]);

            double partial = (differentials[ij] - differentials[ii]) * chain;

            if (normalize) {
                partial -= chain * pi[i] * normalizationGradientContribution;
            }

            total += partial;
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
    public String getDescription() {
        return null; // TODO
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CommonCitations.TAO_LUMPING_2025);
    }

    private static class ParameterDimensionLink {

        private final Parameter parameter;
        private final int index;
        private final List<int[]> entries;

        ParameterDimensionLink(Parameter parameter, int index, List<int[]> entries) {
            this.parameter = parameter;
            this.index = index;
            this.entries = entries;
        }

        public double getDifferential(double element) {
            return element / parameter.getParameterValue(index);
        }
    }
}

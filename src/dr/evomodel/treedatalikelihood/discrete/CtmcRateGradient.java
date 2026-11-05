/*
 * CtmcRateGradient.java
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

import dr.evomodel.substmodel.SVSComplexSubstitutionModel;
import dr.evomodel.treedatalikelihood.GradientDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.loggers.LogColumn;
import dr.inference.model.Parameter;
import dr.util.Citation;
import dr.util.CommonCitations;
import dr.util.Transform;

import java.util.Collections;
import java.util.List;

/**
 * Gradient of the log-likelihood with respect to the (natural-scale) rates of an SVSComplexSubstitutionModel,
 * i.e. a generalSubstitutionModel with a rateIndicator.
 * <p>
 * The normalized generator is Q_ij = indicator_k * rate_k * pi_j / c, so
 * d logL / d rate_k = (d logL / d log rate_k) / rate_k, where the log-rate term is the one LogCtmcRateGradient
 * forms from Q_ij. A rate switched off by its indicator has Q_ij = 0 and so a zero gradient.
 *
 * @author Marc A. Suchard
 */
public class CtmcRateGradient extends AbstractLogAdditiveSubstitutionModelGradient {

    private final Parameter rates;
    private final int[][] mapEffectToIndices;

    public CtmcRateGradient(String traitName,
                            TreeDataLikelihood treeDataLikelihood,
                            GradientDataLikelihoodDelegate likelihoodDelegate,
                            SVSComplexSubstitutionModel substitutionModel,
                            ApproximationMode mode) {
        super(traitName, treeDataLikelihood, likelihoodDelegate, substitutionModel, mode, false);
        this.rates = substitutionModel.getRatesParameter();
        this.mapEffectToIndices = makeAsymmetricMap();
    }

    @Override
    double processSingleGradientDimension(int k, double[] differentials, double[] generator, double[] pi,
                                          boolean normalize, double normalizationGradientContribution,
                                          double normalizationScalar,
                                          Transform transform, boolean scaleByFrequencies) {

        final int i = mapEffectToIndices[k][0], j = mapEffectToIndices[k][1];
        final int ii = i * stateCount + i;
        final int ij = i * stateCount + j;

        double total = (differentials[ij] - differentials[ii]) * generator[ij];

        if (normalize) {
            total -= generator[ij] * pi[i] * normalizationGradientContribution;
        }

        return total / rates.getParameterValue(k);
    }

    @Override
    public Parameter getParameter() {
        return rates;
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
        return Collections.singletonList(CommonCitations.MONTI_GENERIC_RATES_2024);
    }
}

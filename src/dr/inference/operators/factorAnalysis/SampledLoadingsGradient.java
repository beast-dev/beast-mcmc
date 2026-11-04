/*
 * SampledLoadingsGradient.java
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

package dr.inference.operators.factorAnalysis;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.*;
import dr.inferencexml.operators.factorAnalysis.LoadingsOperatorParserUtilities;
import dr.xml.*;


public class SampledLoadingsGradient implements GradientWrtParameterProvider {

    private final MatrixParameterInterface loadings;

    private final FactorAnalysisStatisticsProvider statisticsProvider;
    private final FactorAnalysisOperatorAdaptor adaptor;

    private final double[][] scaledFactorTraitProducts;
    private final double[][][] precisions;
    private final int nFactors;
    private final int nTraits;
    private boolean statisticsKnown = false;

    private Likelihood likelihood;

    SampledLoadingsGradient(FactorAnalysisStatisticsProvider statisticsProvider) {

        this.statisticsProvider = statisticsProvider;
        this.adaptor = statisticsProvider.getAdaptor();
        this.loadings = adaptor.getLoadings();

        this.nFactors = adaptor.getNumberOfFactors();
        this.nTraits = adaptor.getNumberOfTraits();
        this.scaledFactorTraitProducts = new double[nTraits][nFactors];
        this.precisions = new double[nTraits][nFactors][nFactors];

        this.likelihood = new CompoundLikelihood(adaptor.getLikelihoods());
    }

    @Override
    public Likelihood getLikelihood() {
        return likelihood;
    }

    @Override
    public Parameter getParameter() {
        return loadings;
    }

    @Override
    public int getDimension() {
        return loadings.getDimension();
    }

    private void updateStatistics() {
        adaptor.drawFactors();
        for (int i = 0; i < nTraits; i++) {
            statisticsProvider.getScaledFactorInnerProduct(i, nFactors, precisions[i]);
            statisticsProvider.getScaledFactorTraitProduct(i, nFactors, scaledFactorTraitProducts[i]);
        }

        statisticsKnown = true;
    }


    @Override
    public double[] getGradientLogDensity() {
        updateStatistics();

        double[] gradient = new double[getDimension()];

        for (int i = 0; i < nTraits; i++) {
            for (int j = 0; j < nFactors; j++) {
                int index = j * nTraits + i;
                gradient[index] = scaledFactorTraitProducts[i][j];
                for (int k = 0; k < nFactors; k++) {
                    gradient[index] -= precisions[i][j][k] * loadings.getParameterValue(i, k);
                }
            }
        }
        return gradient;
    }


    private static final String SAMPLED_LOADINGS_GRADIENT = "sampledLoadingsGradient";


    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            FactorAnalysisStatisticsProvider statisticsProvider =
                    LoadingsOperatorParserUtilities.parseAdaptorAndStatistics(xo);
            return new SampledLoadingsGradient(statisticsProvider);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return LoadingsOperatorParserUtilities.statisticsProviderRules;
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return SampledLoadingsGradient.class;
        }

        @Override
        public String getParserName() {
            return SAMPLED_LOADINGS_GRADIENT;
        }
    };


}

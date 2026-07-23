/*
 * LoadingsOperatorParserUtilities.java
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

package dr.inferencexml.operators.factorAnalysis;

import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.IntegratedFactorAnalysisLikelihood;
import dr.inference.model.LatentFactorModel;
import dr.inference.operators.factorAnalysis.FactorAnalysisOperatorAdaptor;
import dr.inference.operators.factorAnalysis.FactorAnalysisStatisticsProvider;
import dr.xml.*;


public class LoadingsOperatorParserUtilities {

    private final static String USE_CACHE = "cacheInnerProducts";

    public static FactorAnalysisStatisticsProvider parseAdaptorAndStatistics(XMLObject xo) throws XMLParseException {


        final FactorAnalysisStatisticsProvider.CacheProvider cacheProvider;
        boolean useCache = xo.getAttribute(USE_CACHE, false);
        if (useCache) {
            cacheProvider = FactorAnalysisStatisticsProvider.CacheProvider.USE_CACHE;
        } else {
            cacheProvider = FactorAnalysisStatisticsProvider.CacheProvider.NO_CACHE;
        }

        FactorAnalysisOperatorAdaptor adaptor = parseFactorAnalsysisOperatorAdaptor(xo);

        return new FactorAnalysisStatisticsProvider(adaptor, cacheProvider);
    }


    public static FactorAnalysisOperatorAdaptor parseFactorAnalsysisOperatorAdaptor(XMLObject xo) {
        LatentFactorModel factorModel = (LatentFactorModel) xo.getChild(LatentFactorModel.class);

        if (factorModel == null) {
            IntegratedFactorAnalysisLikelihood integratedModel =
                    (IntegratedFactorAnalysisLikelihood) xo.getChild(IntegratedFactorAnalysisLikelihood.class);
            TreeDataLikelihood treeLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);

            return new FactorAnalysisOperatorAdaptor.IntegratedFactors(integratedModel, treeLikelihood);

        } else {
            return new FactorAnalysisOperatorAdaptor.SampledFactors(factorModel);
        }
    }

    public static final XMLSyntaxRule[] adaptorRules = new XMLSyntaxRule[]{
            new XORRule(
                    new ElementRule(LatentFactorModel.class),
                    new AndRule(
                            new ElementRule(IntegratedFactorAnalysisLikelihood.class),
                            new ElementRule(TreeDataLikelihood.class)
                    )
            )
    };


    public static final XMLSyntaxRule[] statisticsProviderRules =
            XMLSyntaxRule.Utils.concatenate(adaptorRules,
                    new XMLSyntaxRule[]{AttributeRule.newBooleanRule(USE_CACHE, true)});
}

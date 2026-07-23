/*
 * IntegratedLoadingsAndPrecisionGradientParser.java
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

package dr.evomodelxml.continuous.hmc;

import dr.evomodel.continuous.hmc.IntegratedLoadingsAndPrecisionGradient;
import dr.evomodel.continuous.hmc.IntegratedLoadingsGradient;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.evomodel.treedatalikelihood.continuous.IntegratedFactorAnalysisLikelihood;
import dr.inference.model.CompoundParameter;
import dr.util.TaskPool;
import dr.xml.ElementRule;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class IntegratedLoadingsAndPrecisionGradientParser extends IntegratedLoadingsGradientParser {

    public static final String PARSER_NAME = "integratedFactorAnalysisLoadingsAndPrecisionGradient";

    protected IntegratedLoadingsGradient factory(TreeDataLikelihood treeDataLikelihood,
                                                 ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                 IntegratedFactorAnalysisLikelihood factorAnalysisLikelihood,
                                                 ContinuousTraitPartialsProvider jointPartialsProvider,
                                                 TaskPool taskPool,
                                                 IntegratedLoadingsGradient.ThreadUseProvider threadUseProvider,
                                                 IntegratedLoadingsGradient.RemainderCompProvider remainderCompProvider,
                                                 CompoundParameter parameter)
            throws XMLParseException {

        return new IntegratedLoadingsAndPrecisionGradient(
                parameter,
                treeDataLikelihood,
                likelihoodDelegate,
                factorAnalysisLikelihood,
                jointPartialsProvider,
                taskPool,
                threadUseProvider,
                remainderCompProvider);

    }


    @Override
    public String getParserDescription() {
        return "Generates a gradient provider for the loadings matrix & precision when factors are integrated out";
    }

    @Override
    public Class getReturnType() {
        return IntegratedLoadingsAndPrecisionGradient.class;
    }

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        XMLSyntaxRule[] newRules = new XMLSyntaxRule[rules.length + 1];
        newRules[0] = new ElementRule(CompoundParameter.class);
        System.arraycopy(rules, 0, newRules, 1, rules.length);
        return newRules;
    }
}

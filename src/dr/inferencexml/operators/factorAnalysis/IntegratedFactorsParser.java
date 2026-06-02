/*
 * IntegratedFactorsParser.java
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
import dr.inference.operators.factorAnalysis.FactorAnalysisOperatorAdaptor;
import dr.xml.*;

/**
 * @author Gabriel Hassler
 */

public class IntegratedFactorsParser extends AbstractXMLObjectParser {

    private static String INTEGRATED_FACTORS = "integratedFactors";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        TreeDataLikelihood treeDataLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
        IntegratedFactorAnalysisLikelihood factorLikelihood =
                (IntegratedFactorAnalysisLikelihood) xo.getChild(IntegratedFactorAnalysisLikelihood.class);
        return new FactorAnalysisOperatorAdaptor.IntegratedFactors(factorLikelihood, treeDataLikelihood);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(TreeDataLikelihood.class),
                new ElementRule(IntegratedFactorAnalysisLikelihood.class)
        };
    }

    @Override
    public String getParserDescription() {
        return "Class used internally for drawing (latent) factors.";
    }

    @Override
    public Class getReturnType() {
        return FactorAnalysisOperatorAdaptor.IntegratedFactors.class;
    }

    @Override
    public String getParserName() {
        return INTEGRATED_FACTORS;
    }
}

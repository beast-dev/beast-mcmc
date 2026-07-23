/*
 * IntegratedFactorGibbsOperator.java
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

import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.IntegratedFactorAnalysisLikelihood;
import dr.inference.model.CompoundParameter;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.xml.*;

public class IntegratedFactorGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {

    private final CompoundParameter factors;
    private final FactorAnalysisOperatorAdaptor adaptor;

    IntegratedFactorGibbsOperator(CompoundParameter factors, TreeDataLikelihood treeLikelihood,
                                  IntegratedFactorAnalysisLikelihood integratedFactorLikelihood, double weight) {
        setWeight(weight);
        this.factors = factors;
        this.adaptor = new FactorAnalysisOperatorAdaptor.IntegratedFactors(integratedFactorLikelihood, treeLikelihood);
    }

    @Override
    public String getOperatorName() {
        return INTEGRATED_FACTORS_GIBBS;
    }

    @Override
    public double doOperation() {
        adaptor.fireLoadingsChanged();
        adaptor.drawFactors();
        for (int i = 0; i < adaptor.getNumberOfTaxa(); i++) {
            for (int j = 0; j < adaptor.getNumberOfFactors(); j++) {
                factors.getParameter(i).setParameterValueQuietly(j, adaptor.getFactorValue(j, i));
//                factors.setParameterValueQuietly(j, i, adaptor.getFactorValue(j, i));
            }
        }
        factors.fireParameterChangedEvent();
        return 0;
    }

    private static final String INTEGRATED_FACTORS_GIBBS = "integratedFactorsGibbsOperator";

    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            CompoundParameter factors = (CompoundParameter) xo.getChild(CompoundParameter.class);
            TreeDataLikelihood treeDataLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
            IntegratedFactorAnalysisLikelihood factorLikelihood =
                    (IntegratedFactorAnalysisLikelihood) xo.getChild(IntegratedFactorAnalysisLikelihood.class);
            double weight = xo.getDoubleAttribute(WEIGHT);
            return new IntegratedFactorGibbsOperator(factors, treeDataLikelihood, factorLikelihood, weight);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[] {
                    new ElementRule(CompoundParameter.class),
                    new ElementRule(TreeDataLikelihood.class),
                    new ElementRule(IntegratedFactorAnalysisLikelihood.class),
                    AttributeRule.newDoubleRule(WEIGHT)
            };
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return IntegratedFactorGibbsOperator.class;
        }

        @Override
        public String getParserName() {
            return INTEGRATED_FACTORS_GIBBS;
        }
    };
}

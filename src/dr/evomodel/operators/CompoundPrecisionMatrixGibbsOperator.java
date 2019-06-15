/*
 * CompoundPrecisionMatrixGibbsOperator.java
 *
 * Copyright (c) 2002-2018 Alexei Drummond, Andrew Rambaut and Marc Suchard
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


package dr.evomodel.operators;

import dr.evomodel.treedatalikelihood.continuous.RepeatedMeasuresWishartStatistics;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Gabriel Hassler
 */


public class CompoundPrecisionMatrixGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {

    private final PrecisionMatrixGibbsOperator diffusionOperator;
    private final PrecisionMatrixGibbsOperator residualOperator;


    CompoundPrecisionMatrixGibbsOperator(PrecisionMatrixGibbsOperator diffusionOperator,
                                         PrecisionMatrixGibbsOperator residualOperator) {

        this.diffusionOperator = diffusionOperator;
        this.residualOperator = residualOperator;

    }


    @Override
    public double doOperation() {
        diffusionOperator.doOperationDontFireChange();
        residualOperator.doOperationDontFireChange();
        diffusionOperator.getPrecisionParam().fireParameterChangedEvent();
        residualOperator.getPrecisionParam().fireParameterChangedEvent();
        return 0;
    }

    @Override
    public String getOperatorName() {
        return COMPOUND_OPERATOR;
    }


    private static final String DIFFUSION_OPERATOR = "diffusionOperator";
    private static final String RESIDUAL_OPERATOR = "residualOperator";
    private static final String COMPOUND_OPERATOR = "compoundPrecisionOperator";

    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {


        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            XMLObject diffXO = xo.getChild(DIFFUSION_OPERATOR);

            PrecisionMatrixGibbsOperator diffusionOperator =
                    (PrecisionMatrixGibbsOperator) diffXO.getChild(PrecisionMatrixGibbsOperator.class);

            XMLObject resXO = xo.getChild(RESIDUAL_OPERATOR);

            PrecisionMatrixGibbsOperator residualOperator =
                    (PrecisionMatrixGibbsOperator) resXO.getChild(PrecisionMatrixGibbsOperator.class);

            RepeatedMeasuresWishartStatistics residualProvider =
                    (RepeatedMeasuresWishartStatistics) residualOperator.getConjugateWishartProvider();

            residualProvider.setForceResample(false);


            return new CompoundPrecisionMatrixGibbsOperator(diffusionOperator, residualOperator);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(DIFFUSION_OPERATOR, new XMLSyntaxRule[]{
                            new ElementRule(PrecisionMatrixGibbsOperator.class)
                    }),

                    new ElementRule(RESIDUAL_OPERATOR, new XMLSyntaxRule[]{
                            new ElementRule(PrecisionMatrixGibbsOperator.class)
                    })

            };
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return CompoundPrecisionMatrixGibbsOperator.class;
        }

        @Override
        public String getParserName() {
            return COMPOUND_OPERATOR;
        }
    };

}

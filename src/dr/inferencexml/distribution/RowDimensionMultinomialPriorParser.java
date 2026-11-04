/*
 * RowDimensionMultinomialPriorParser.java
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

package dr.inferencexml.distribution;

import dr.inference.distribution.RowDimensionMultinomialPrior;
import dr.inference.model.AdaptableSizeFastMatrixParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Max Tolkoff
 */
public class RowDimensionMultinomialPriorParser extends AbstractXMLObjectParser {
    public static final String ROW_DIMENSION_MULTINOMIAL_PRIOR = "rowDimensionMultinomialPrior";
    public static final String PROBABILITIES = "probabilities";
    public static final String TRANSPOSE = "transpose";


    @Override
    public String getParserName() {
        return ROW_DIMENSION_MULTINOMIAL_PRIOR;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        AdaptableSizeFastMatrixParameter data = (AdaptableSizeFastMatrixParameter) xo.getChild(AdaptableSizeFastMatrixParameter.class);
        Parameter probabilities = (Parameter) xo.getChild(PROBABILITIES).getChild(Parameter.class);
        boolean transpose = false;
        if(xo.hasAttribute(TRANSPOSE))
            transpose = xo.getBooleanAttribute(TRANSPOSE);
        String id = xo.getId();


        return new RowDimensionMultinomialPrior(id, data, probabilities, transpose);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(PROBABILITIES,
                    new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    }
            ),
            new ElementRule(AdaptableSizeFastMatrixParameter.class),
            AttributeRule.newBooleanRule(TRANSPOSE, true),
    };


    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return RowDimensionMultinomialPrior.class;
    }
}

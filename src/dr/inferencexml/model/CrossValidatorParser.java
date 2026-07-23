/*
 * CrossValidatorParser.java
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

package dr.inferencexml.model;

import dr.inference.model.CrossValidationProvider;
import dr.xml.*;

/**
 * @author Gabriel Hassler
 */

public class CrossValidatorParser extends AbstractXMLObjectParser {

    public final static String LOG_SUM = "logSum";
    private final static String TYPE = "type";
    public final static String CROSS_VALIDATION = "crossValidation";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        CrossValidationProvider provider = (CrossValidationProvider) xo.getChild(CrossValidationProvider.class);
        boolean logSum = xo.getAttribute(LOG_SUM, false);

        CrossValidationProvider.ValidationType validationType;

        String validation = xo.getAttribute(TYPE, CrossValidationProvider.ValidationType.SQUARED_ERROR.getName());

        if (validation.equalsIgnoreCase(CrossValidationProvider.ValidationType.SQUARED_ERROR.getName())) {
            validationType = CrossValidationProvider.ValidationType.SQUARED_ERROR;
        } else if (validation.equalsIgnoreCase(CrossValidationProvider.ValidationType.BIAS.getName())) {
            validationType = CrossValidationProvider.ValidationType.BIAS;
        } else if (validation.equalsIgnoreCase(CrossValidationProvider.ValidationType.VALUE.getName())) {
            validationType = CrossValidationProvider.ValidationType.VALUE;
        } else {
            throw new XMLParseException("The attribute '" + TYPE + "' can only take values '" +
                    CrossValidationProvider.ValidationType.SQUARED_ERROR.getName() + "', " +
                    CrossValidationProvider.ValidationType.BIAS.getName() + "', or" +
                    CrossValidationProvider.ValidationType.VALUE.getName() + "'.");
        }

        if (logSum) return new CrossValidationProvider.CrossValidatorSum(provider, validationType);
        return new CrossValidationProvider.CrossValidator(provider, validationType);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newBooleanRule(LOG_SUM, true),
                AttributeRule.newStringRule(TYPE, true),
                new ElementRule(CrossValidationProvider.class)
        };
    }

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return CrossValidationProvider.CrossValidator.class;
    }

    @Override
    public String getParserName() {
        return CROSS_VALIDATION;
    }
}

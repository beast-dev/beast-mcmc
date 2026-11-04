/*
 * ElementaryVectorDataModelParser.java
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

package dr.evomodelxml.continuous;

import dr.evomodel.treedatalikelihood.continuous.ElementaryVectorDataModel;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */

public class ElementaryVectorDataModelParser extends AbstractXMLObjectParser {

    private static final String ELEMENTARY_DATA_MODEL = "elementaryVectorDataModel";
    private static final String PRECISION_TYPE = "precisionType";
    private static final String DIMENSION = "dim";
    private static final String TIPS = "tips";

    private static final String TIP_INDICATOR = "tipIndicator";
    private static final String DIM_INDICATOR = "dimIndicator";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter tipIndicator = (Parameter) xo.getElementFirstChild(TIP_INDICATOR);

        Parameter dimIndicator = null;
        if (xo.hasChildNamed(DIM_INDICATOR)) {
            dimIndicator = (Parameter) xo.getElementFirstChild(DIM_INDICATOR);
        }

        int tips = parsePositiveInteger(xo, TIPS);
        int dim = parsePositiveInteger(xo, DIMENSION);

        PrecisionType precisionType = parsePrecisionType(xo);

        return new ElementaryVectorDataModel(xo.getId(),
                tipIndicator, dimIndicator,
                tips, dim,
                precisionType);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return ElementaryVectorDataModel.class;
    }

    @Override
    public String getParserName() {
        return ELEMENTARY_DATA_MODEL;
    }

    private static PrecisionType parsePrecisionType(XMLObject xo) throws XMLParseException {

        String tag = xo.getAttribute(PRECISION_TYPE, PrecisionType.ELEMENTARY.getTag());
        if (tag.compareTo(PrecisionType.ELEMENTARY.getTag()) == 0) {
            return PrecisionType.ELEMENTARY;
        } else if (tag.compareTo(PrecisionType.SCALAR.getTag()) == 0) {
            return PrecisionType.SCALAR;
        } else {
            throw new XMLParseException("Invalid precision type");
        }
    }

    private static int parsePositiveInteger(XMLObject xo, String name) throws XMLParseException {

        int value = xo.getIntegerAttribute(name);
        if (value < 1) {
            throw new XMLParseException("Integer attribute '" + name  + "' must be at least 1");
        }

        return value;
    }

    private final static XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(TIP_INDICATOR, new XMLSyntaxRule[] {
                    new ElementRule(Parameter.class),
            }),
            new ElementRule(DIM_INDICATOR, new XMLSyntaxRule[] {
                    new ElementRule(Parameter.class),
            }, true),
            AttributeRule.newStringRule(PRECISION_TYPE, true),
            AttributeRule.newIntegerRule(TIPS),
            AttributeRule.newIntegerRule(DIMENSION),
    };
}

/*
 * MaskedGradientParser.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.inferencexml.hmc;

import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */
public class RotationTranslationMaskParser extends AbstractXMLObjectParser {

    private final static String MASK = "rotationalTranslationalMask";
    private final static String DIMENSION = "dimension";

    @Override
    public String getParserName() {
        return MASK;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter parameter = (Parameter) xo.getChild(Parameter.class);
        int dim = xo.getIntegerAttribute(DIMENSION);

        if (parameter.getDimension() % dim != 0) {
            throw new XMLParseException("Dimension and parameter length are not divisible");
        }

        Parameter mask = new Parameter.Default(parameter.getDimension(), 1.0);

        int offset = 0;

        // Translational invariance
        for (int i = 0; i < dim; ++i) {
            mask.setParameterValue(offset, 0.0);
            ++offset;
        }

        // Rotational invariance
        ++offset;
        for (int i = 1; i < dim; ++i) {
            mask.setParameterValue(offset, 0.0);
            ++offset;
        }

        return mask;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule(DIMENSION),
            new ElementRule(Parameter.class),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return Parameter.class;
    }
}

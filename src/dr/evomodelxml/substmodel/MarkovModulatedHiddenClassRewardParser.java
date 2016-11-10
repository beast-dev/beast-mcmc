/*
 * MarkovModulatedHiddenClassRewardParser.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.substmodel;

import dr.evomodel.substmodel.MarkovModulatedSubstitutionModel;
import dr.evolution.datatype.HiddenDataType;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */

public class MarkovModulatedHiddenClassRewardParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "hiddenClassRewardParameter";
    public static final String NAME = "name";
    public static final String CLASS_NUMBER = "class";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        MarkovModulatedSubstitutionModel substitutionModel = (MarkovModulatedSubstitutionModel) xo.getChild(MarkovModulatedSubstitutionModel.class);
        HiddenDataType hiddenDataType = (HiddenDataType) substitutionModel.getDataType();
        int classNumber = xo.getIntegerAttribute(CLASS_NUMBER);
        int hiddenClassCount = hiddenDataType.getHiddenClassCount();
        if (classNumber < 1 || classNumber > hiddenClassCount) {
            throw new XMLParseException("Invalid class number in " + xo.getId());
        }
        classNumber--; // Use zero-indexed number
        int stateCount = hiddenDataType.getStateCount() / hiddenClassCount;

        // Construct reward parameter
        Parameter parameter = new Parameter.Default(stateCount * hiddenClassCount, 0.0);
        for (int i = 0; i < stateCount; ++i) {
            parameter.setParameterValue(i + classNumber * stateCount, 1.0);
        }

        if (xo.hasAttribute(NAME)) {
            parameter.setId((String) xo.getAttribute(NAME));
        } else {
            parameter.setId(substitutionModel.getId() + "_" + Integer.toString(classNumber + 1));
        }

        return parameter;
    }

    /**
     * @return an array of syntax rules required by this element.
     *         Order is not important.
     */
    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = {
            AttributeRule.newIntegerRule(CLASS_NUMBER),
            AttributeRule.newStringRule(NAME, true),
            new ElementRule(MarkovModulatedSubstitutionModel.class),
    };

    @Override
    public String getParserDescription() {
        return "Generates a reward parameter to log hidden classes in Markov-modulated substitutionProcess";
    }

    @Override
    public Class getReturnType() {
        return Parameter.class;
    }

    /**
     * @return Parser name, which is identical to name of xml element parsed by it.
     */
    public String getParserName() {
        return PARSER_NAME;
    }
}

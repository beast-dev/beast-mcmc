/*
 * BirthDeathSubstitutionModelParser.java
 *
 * Copyright (c) 2002-2022 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.evolution.datatype.DataType;
import dr.evomodel.substmodel.*;
import dr.inference.model.Parameter;
import dr.inference.model.ParameterParser;
import dr.xml.*;

import java.util.Arrays;
import java.util.logging.Logger;

/**
 * @author Marc Suchard
 */
public class BirthDeathSubstitutionModelParser extends AbstractXMLObjectParser {

    private static final String SUBSTITUTION_MODEL = "birthDeathSubstitutionModel";
    private static final String BIRTH_PARAMETER = "birth";
    private static final String DEATH_PARAMETER = "death";
    private static final String USE_STATIONARY_DISTRIBUTION = "useStationaryDistribution";
    private static final String FREQUENCIES = "frequencies";

    public String getParserName() {
        return SUBSTITUTION_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter birth = (Parameter) xo.getElementFirstChild(BIRTH_PARAMETER);
        Parameter death = (Parameter) xo.getElementFirstChild(DEATH_PARAMETER);
        DataType  dataType = (DataType) xo.getChild(DataType.class);
        int states = dataType.getStateCount();

        boolean useStationaryDistribution = xo.getAttribute(USE_STATIONARY_DISTRIBUTION, false);

        Logger.getLogger("dr.app.beagle.evomodel").info(
                "  Birth-death Substitution Model (stateCount=" + states + ")");

        if (useStationaryDistribution) {
            Logger.getLogger("dr.app.beagle.evomodel").info(
                    "    using stationary distribution of process as root frequencies");
        }

        BirthDeathSubstitutionModel model = new BirthDeathSubstitutionModel(xo.getId(),
                Arrays.asList(birth, death), dataType, useStationaryDistribution);

        if (xo.hasChildNamed(FREQUENCIES)) {
            ParameterParser.replaceParameter(xo.getChild(FREQUENCIES),
                    model.getFrequencyModel().getFrequencyParameter());
        }

        return model;
    }

    public String getParserDescription() {
        return "A general (truncated) birth-death substitution model for counting number patterns.";
    }

    public Class getReturnType() {
        return SVSComplexSubstitutionModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(BIRTH_PARAMETER, new XMLSyntaxRule[] {
                    new ElementRule(Parameter.class),
            }),
            new ElementRule(DEATH_PARAMETER, new XMLSyntaxRule[] {
                    new ElementRule(Parameter.class),
            }),
            new ElementRule(FREQUENCIES, new XMLSyntaxRule[] {
                    new ElementRule(Parameter.class),
            }, true),
            AttributeRule.newBooleanRule(USE_STATIONARY_DISTRIBUTION, true),
            new ElementRule(DataType.class),
    };
}

/*
 * BirthDeathCompoundParameterParser.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodelxml.speciation;

import dr.evomodel.speciation.NewBirthDeathSerialSamplingModel;
import dr.evomodel.speciation.BirthDeathCompoundParameterLogger;
import dr.xml.*;

public class BirthDeathCompoundParameterLoggerParser extends AbstractXMLObjectParser {

    public static final String BDSS_COMPOUND_PARAMETER = "birthDeathCompoundParameterLogger";
    public static final String TYPE = "compoundParameterType";


    public String getParserName() {
        return BDSS_COMPOUND_PARAMETER;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        NewBirthDeathSerialSamplingModel bdss = (NewBirthDeathSerialSamplingModel) xo.getChild(NewBirthDeathSerialSamplingModel.class);

        BirthDeathCompoundParameterLogger.BDPCompoundParameterType type = parseFromString(xo.getStringAttribute(TYPE));

        return new BirthDeathCompoundParameterLogger(bdss, type);
    }

    public  BirthDeathCompoundParameterLogger.BDPCompoundParameterType parseFromString(String text) throws XMLParseException {
        for (BirthDeathCompoundParameterLogger.BDPCompoundParameterType type : BirthDeathCompoundParameterLogger.BDPCompoundParameterType.values()) {
            if (type.getName().toLowerCase().compareToIgnoreCase(text) == 0) {
                return type;
            }
        }
        throw new XMLParseException("Unknown type '" + text + "'");
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "Computes compound parameters from birth-death serial sampling models automatically (e.g. the effective reproductive number).";
    }

    public Class getReturnType() {
        return BirthDeathCompoundParameterLogger.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(NewBirthDeathSerialSamplingModel.class, "The birth-death model containing the natural parameters (lambda, mu, psi, r)."),
            new StringAttributeRule(TYPE,"What compound parameter should be computed? Allowed: effectiveReproductiveNumber.")
    };
}

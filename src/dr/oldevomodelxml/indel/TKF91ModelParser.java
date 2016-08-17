/*
 * TKF91ModelParser.java
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

package dr.oldevomodelxml.indel;

import dr.evolution.util.Units;
import dr.oldevomodel.indel.TKF91Model;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a ConstantPopulation.
 */
public class TKF91ModelParser extends AbstractXMLObjectParser {

    public static final String TKF91_MODEL = "tkf91Model";
    public static final String TKF91_LENGTH_DIST = "lengthDistribution";
    public static final String BIRTH_RATE = "birthRate";
    public static final String DEATH_RATE = "deathRate";

    public String getParserName() {
        return TKF91_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Parameter lengthDistParameter = (Parameter) xo.getElementFirstChild(TKF91_LENGTH_DIST);
        Parameter deathParameter = (Parameter) xo.getElementFirstChild(DEATH_RATE);
        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        return new TKF91Model(lengthDistParameter, deathParameter, units);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "The TKF91 (Thorne, Kishino & Felsenstein 1991) model of insertion-deletion.";
    }

    public Class getReturnType() {
        return TKF91Model.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(TKF91_LENGTH_DIST,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            new ElementRule(DEATH_RATE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
            XMLUnits.SYNTAX_RULES[0]
    };
}

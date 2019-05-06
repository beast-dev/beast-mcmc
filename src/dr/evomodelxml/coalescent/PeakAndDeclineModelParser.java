/*
 * PeakAndDeclineModelParser.java
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

package dr.evomodelxml.coalescent;

import dr.evolution.util.Units;
import dr.evomodel.coalescent.demographicmodels.LogisticGrowthModel;
import dr.evomodel.coalescent.demographicmodels.PeakAndDeclineModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an XMLObject into LogisticGrowthModel.
 */
public class PeakAndDeclineModelParser extends AbstractXMLObjectParser {

    public static String PEAK_VALUE = "peakValue";
    public static String PEAK_AND_DECLINE_MODEL = "peakAndDecline";

    public static String SHAPE = "shape";
    public static String PEAK_TIME = "peakTime";

    public String getParserName() {
        return PEAK_AND_DECLINE_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        XMLObject cxo = xo.getChild(PEAK_VALUE);
        Parameter peakValueParam = (Parameter) cxo.getChild(Parameter.class);


        Parameter rParam;

        cxo = xo.getChild(SHAPE);
        rParam = (Parameter) cxo.getChild(Parameter.class);


        cxo = xo.getChild(PEAK_TIME);
        Parameter peakTimeParam = (Parameter) cxo.getChild(Parameter.class);

        return new PeakAndDeclineModel(peakValueParam, rParam, peakTimeParam, units);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "Logistic growth demographic model.";
    }

    public Class getReturnType() {
        return LogisticGrowthModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            XMLUnits.SYNTAX_RULES[0],
            new ElementRule(PEAK_VALUE,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),


                    new ElementRule(SHAPE,
                            new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),


            new ElementRule(PEAK_TIME,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)})
    };
}

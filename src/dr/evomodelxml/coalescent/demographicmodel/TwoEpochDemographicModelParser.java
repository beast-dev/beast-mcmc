/*
 * TwoEpochDemographicModelParser.java
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

package dr.evomodelxml.coalescent.demographicmodel;

import dr.evolution.util.Units;
import dr.evomodel.coalescent.demographicmodel.DemographicModel;
import dr.evomodel.coalescent.demographicmodel.TwoEpochDemographicModel;
import dr.evoxml.util.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Parses an element from an DOM document into a ExponentialGrowth.
 */
public class TwoEpochDemographicModelParser extends AbstractXMLObjectParser {

    public static final String EPOCH_1 = "modernEpoch";
    public static final String EPOCH_2 = "ancientEpoch";
    public static final String TRANSITION_TIME = "transitionTime";

    public static final String TWO_EPOCH_MODEL = "twoEpoch";

    public String getParserName() {
        return TWO_EPOCH_MODEL;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        Units.Type units = XMLUnits.Utils.getUnitsAttr(xo);

        XMLObject cxo = xo.getChild(EPOCH_1);
        DemographicModel demo1 = (DemographicModel) cxo.getChild(DemographicModel.class);

        cxo = xo.getChild(EPOCH_2);
        DemographicModel demo2 = (DemographicModel) cxo.getChild(DemographicModel.class);

        cxo = xo.getChild(TRANSITION_TIME);
        Parameter timeParameter = (Parameter) cxo.getChild(Parameter.class);

        return new TwoEpochDemographicModel(demo1, demo2, timeParameter, units);
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "A demographic model of two epochs.";
    }

    public Class getReturnType() {
        return TwoEpochDemographicModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(EPOCH_1,
                    new XMLSyntaxRule[]{new ElementRule(DemographicModel.class)},
                    "The demographic model for the recent epoch."),
            new ElementRule(EPOCH_2,
                    new XMLSyntaxRule[]{new ElementRule(DemographicModel.class)},
                    "The demographic model for the ancient epoch."),
            new ElementRule(TRANSITION_TIME,
                    new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                    "The time that splits the two epochs."),
            XMLUnits.SYNTAX_RULES[0]
    };
}

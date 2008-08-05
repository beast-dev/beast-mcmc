/*
 * TwoEpochDemographicModel.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.evomodel.coalescent;

import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.coalescent.TwoEpochDemographic;
import dr.evoxml.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * This class models an exponentially growing model that suffers a
 * cataclysmic event and goes into exponential decline
 * This model is nested with the constant-population size model (r=0).
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: TwoEpochDemographicModel.java,v 1.6 2005/05/24 20:25:57 rambaut Exp $
 */
public class TwoEpochDemographicModel extends DemographicModel {
    //
    // Public stuff
    //

    public static String EPOCH_1 = "modernEpoch";
    public static String EPOCH_2 = "ancientEpoch";
    public static String TRANSITION_TIME = "transitionTime";

    public static String TWO_EPOCH_MODEL = "twoEpoch";

    /**
     * Construct demographic model with default settings.
     */
    public TwoEpochDemographicModel(DemographicModel demo1, DemographicModel demo2, Parameter transitionTimeParameter, Type units) {

        this(TWO_EPOCH_MODEL, demo1, demo2, transitionTimeParameter, units);
    }

    /**
     * Construct demographic model with default settings.
     */
    public TwoEpochDemographicModel(String name, DemographicModel demo1, DemographicModel demo2, Parameter transitionTimeParameter, Type units) {

        super(name);

        this.demo1 = demo1;
        addModel(demo1);
        for (int i = 0; i < demo1.getParameterCount(); i++) {
            addParameter(demo1.getParameter(i));
        }

        this.demo2 = demo2;
        addModel(demo2);
        for (int i = 0; i < demo2.getParameterCount(); i++) {
            addParameter(demo2.getParameter(i));
        }

        this.transitionTimeParameter = transitionTimeParameter;
        addParameter(transitionTimeParameter);

        setUnits(units);
    }

    // general functions

    public DemographicFunction getDemographicFunction() {

        TwoEpochDemographic twoEpoch = new TwoEpochDemographic(demo1.getDemographicFunction(), demo2.getDemographicFunction(), getUnits());
        twoEpoch.setTransitionTime(transitionTimeParameter.getParameterValue(0));

        return twoEpoch;
    }

    /**
     * Parses an element from an DOM document into a ExponentialGrowth.
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return TWO_EPOCH_MODEL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Type units = XMLParser.Utils.getUnitsAttr(xo);

            XMLObject cxo = (XMLObject) xo.getChild(EPOCH_1);
            DemographicModel demo1 = (DemographicModel) cxo.getChild(DemographicModel.class);

            cxo = (XMLObject) xo.getChild(EPOCH_2);
            DemographicModel demo2 = (DemographicModel) cxo.getChild(DemographicModel.class);

            cxo = (XMLObject) xo.getChild(TRANSITION_TIME);
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
    };

    private Parameter transitionTimeParameter = null;
    private DemographicModel demo1 = null;
    private DemographicModel demo2 = null;
}

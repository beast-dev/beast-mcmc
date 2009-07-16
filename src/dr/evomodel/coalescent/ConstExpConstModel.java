/*
 * ConstExpConstModel.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.coalescent;

import dr.evolution.coalescent.ConstExpConst;
import dr.evolution.coalescent.DemographicFunction;
import dr.evoxml.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Exponential growth from a constant ancestral population size.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: ConstantExponentialModel.java,v 1.8 2005/10/28 02:49:17 alexei Exp $
 */
public class ConstExpConstModel extends DemographicModel {

    //
    // Public stuff
    //

    public static String CONST_EXP_CONST_MODEL = "constExpConst";
    public static String POPULATION_SIZE = "populationSize";
    public static String ANCESTRAL_POPULATION_SIZE = "ancestralPopulationSize";
    public static String FINAL_PHASE_START_TIME = "finalPhaseStartTime";

    public static String GROWTH_RATE = "growthRate";
    public static String DOUBLING_TIME = "doublingTime";

    /**
     * Construct demographic model with default settings
     */
    public ConstExpConstModel(Parameter N0Parameter, Parameter N1Parameter, Parameter timeParameter,
                              Parameter growthRateParameter, Type units, boolean usingGrowthRate) {

        this(CONST_EXP_CONST_MODEL, N0Parameter, N1Parameter, timeParameter, growthRateParameter, units, usingGrowthRate);
    }

    /**
     * Construct demographic model with default settings
     */
    public ConstExpConstModel(String name, Parameter N0Parameter, Parameter N1Parameter, Parameter timeParameter,
                              Parameter growthRateParameter, Type units, boolean usingGrowthRate) {

        super(name);

        constExpConst = new ConstExpConst(units);

        this.N0Parameter = N0Parameter;
        addParameter(N0Parameter);
        N0Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.N1Parameter = N0Parameter;
        addParameter(N1Parameter);
        N1Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.timeParameter = timeParameter;
        addParameter(timeParameter);
        timeParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.growthRateParameter = growthRateParameter;
        addParameter(growthRateParameter);
        growthRateParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.usingGrowthRate = usingGrowthRate;

        setUnits(units);
    }


    // general functions

    public DemographicFunction getDemographicFunction() {

        double time = timeParameter.getParameterValue(0);
        double N0 = N0Parameter.getParameterValue(0);
        double N1 = N1Parameter.getParameterValue(0);
        double growthRate = growthRateParameter.getParameterValue(0);

        if (!usingGrowthRate) {
            double doublingTime = growthRate;
            growthRate = Math.log(2) / doublingTime;
        }

        constExpConst.setGrowthRate(growthRate);
        constExpConst.setN0(N0);
        constExpConst.setN1(N1);
        constExpConst.setTime1(time);

        return constExpConst;
    }

    /**
     * Parses an element from an DOM document into a ConstantExponentialModel.
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return CONST_EXP_CONST_MODEL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Type units = XMLUnits.Utils.getUnitsAttr(xo);

            XMLObject cxo = xo.getChild(POPULATION_SIZE);
            Parameter N0Param = (Parameter) cxo.getChild(Parameter.class);

            cxo = xo.getChild(ANCESTRAL_POPULATION_SIZE);
            Parameter N1Param = (Parameter) cxo.getChild(Parameter.class);

            cxo = xo.getChild(FINAL_PHASE_START_TIME);
            Parameter timeParam = (Parameter) cxo.getChild(Parameter.class);

            Parameter rParam;
            boolean usingGrowthRate = true;

            if (xo.getChild(GROWTH_RATE) != null) {
                cxo = xo.getChild(GROWTH_RATE);
                rParam = (Parameter) cxo.getChild(Parameter.class);
            } else {
                cxo = xo.getChild(DOUBLING_TIME);
                rParam = (Parameter) cxo.getChild(Parameter.class);
                usingGrowthRate = false;
            }

            return new ConstExpConstModel(N0Param, N1Param, timeParam, rParam, units, usingGrowthRate);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A demographic model of constant population size followed by exponential growth.";
        }

        public Class getReturnType() {
            return ConstantExponentialModel.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                XMLUnits.SYNTAX_RULES[0],
                new ElementRule(POPULATION_SIZE,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                new ElementRule(ANCESTRAL_POPULATION_SIZE,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                new ElementRule(FINAL_PHASE_START_TIME,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                new XORRule(

                        new ElementRule(GROWTH_RATE,
                                new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                                "A value of zero represents a constant population size, negative values represent decline towards the present, " +
                                        "positive numbers represents exponential growth towards the present. " +
                                        "A random walk operator is recommended for this parameter with a starting value of 0.0 and no upper or lower limits."),
                        new ElementRule(DOUBLING_TIME,
                                new XMLSyntaxRule[]{new ElementRule(Parameter.class)},
                                "This parameter determines the doubling time.")
                )
        };
    };
    //
    // protected stuff
    //

    private final Parameter N0Parameter;
    private final Parameter N1Parameter;
    private final Parameter timeParameter;
    private final Parameter growthRateParameter;
    private final ConstExpConst constExpConst;
    private final boolean usingGrowthRate;
}
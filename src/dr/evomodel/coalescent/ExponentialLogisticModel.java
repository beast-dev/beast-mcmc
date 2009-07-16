/*
 * ExponentialLogisticModel.java
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

import dr.evolution.coalescent.DemographicFunction;
import dr.evolution.coalescent.ExponentialLogistic;
import dr.evoxml.XMLUnits;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * Exponential growth followed by Logistic growth.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id$
 */
public class ExponentialLogisticModel extends DemographicModel {

    //
    // Public stuff
    //

    public static String EXPONENTIAL_LOGISTIC_MODEL = "exponentialLogistic";
    public static String POPULATION_SIZE = "populationSize";
    public static String TRANSITION_TIME = "transitionTime";

    public static String LOGISTIC_GROWTH_RATE = "logisticGrowthRate";
    public static String LOGISTIC_SHAPE = "logisticShape";
    public static String EXPONENTIAL_GROWTH_RATE = "exponentialGrowthRate";

    public static String ALPHA = "alpha";

    /**
     * Construct demographic model with default settings
     */
    public ExponentialLogisticModel(Parameter N0Parameter,
                                    Parameter logisticGrowthParameter,
                                    Parameter logisticShapeParameter,
                                    Parameter exponentialGrowthParameter,
                                    Parameter transitionTimeParameter,
                                    double alpha, Type units) {

        this(EXPONENTIAL_LOGISTIC_MODEL,
                N0Parameter,
                logisticGrowthParameter,
                logisticShapeParameter,
                exponentialGrowthParameter,
                transitionTimeParameter,
                alpha, units);
    }

    /**
     * Construct demographic model with default settings
     */
    public ExponentialLogisticModel(String name, Parameter N0Parameter,
                                    Parameter logisticGrowthParameter,
                                    Parameter logisticShapeParameter,
                                    Parameter exponentialGrowthParameter,
                                    Parameter transistionTimeParameter,
                                    double alpha, Type units) {

        super(name);

        exponentialLogistic = new ExponentialLogistic(units);

        this.N0Parameter = N0Parameter;
        addParameter(N0Parameter);
        N0Parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.logisticGrowthParameter = logisticGrowthParameter;
        addParameter(logisticGrowthParameter);
        logisticGrowthParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.logisticShapeParameter = logisticShapeParameter;
        addParameter(logisticShapeParameter);
        logisticShapeParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.exponentialGrowthParameter = exponentialGrowthParameter;
        addParameter(exponentialGrowthParameter);
        exponentialGrowthParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.transistionTimeParameter = transistionTimeParameter;
        addParameter(transistionTimeParameter);
        transistionTimeParameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));

        this.alpha = alpha;

        setUnits(units);
    }


    // general functions

    public DemographicFunction getDemographicFunction() {
        exponentialLogistic.setN0(N0Parameter.getParameterValue(0));

        double r = logisticGrowthParameter.getParameterValue(0);
        exponentialLogistic.setGrowthRate(r);

        double r1 = exponentialGrowthParameter.getParameterValue(0);
        exponentialLogistic.setR1(r1);

        double t = transistionTimeParameter.getParameterValue(0);
        exponentialLogistic.setTime(t);

        // logisticGrowth.setShape(Math.exp(shapeParameter.getParameterValue(0)));
        exponentialLogistic.setTime50(logisticShapeParameter.getParameterValue(0));
        //exponentialLogistic.setShapeFromTimeAtAlpha(logisticShapeParameter.getParameterValue(0), alpha);

        return exponentialLogistic;
    }

    /**
     * Parses an element from an DOM document into a ExponentialGrowth.
     */
    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return EXPONENTIAL_LOGISTIC_MODEL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Type units = XMLUnits.Utils.getUnitsAttr(xo);

            XMLObject cxo = xo.getChild(POPULATION_SIZE);
            Parameter N0Param = (Parameter) cxo.getChild(Parameter.class);

            cxo = xo.getChild(LOGISTIC_GROWTH_RATE);
            Parameter logGrowthParam = (Parameter) cxo.getChild(Parameter.class);

            cxo = xo.getChild(LOGISTIC_SHAPE);
            Parameter shapeParam = (Parameter) cxo.getChild(Parameter.class);

            cxo = xo.getChild(EXPONENTIAL_GROWTH_RATE);
            Parameter expGrowthParam = (Parameter) cxo.getChild(Parameter.class);

            cxo = xo.getChild(TRANSITION_TIME);
            Parameter timeParam = (Parameter) cxo.getChild(Parameter.class);

            double alpha = xo.getDoubleAttribute(ALPHA);

            return new ExponentialLogisticModel(N0Param, logGrowthParam, shapeParam, expGrowthParam, timeParam, alpha, units);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A demographic model of constant population size followed by logistic growth.";
        }

        public Class getReturnType() {
            return ConstantLogisticModel.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                XMLUnits.SYNTAX_RULES[0],
                new ElementRule(POPULATION_SIZE,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                new ElementRule(LOGISTIC_GROWTH_RATE,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                new ElementRule(LOGISTIC_SHAPE,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                new ElementRule(EXPONENTIAL_GROWTH_RATE,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                new ElementRule(TRANSITION_TIME,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                AttributeRule.newDoubleRule(ALPHA)
        };
    };
    //
    // protected stuff
    //

    Parameter N0Parameter = null;
    Parameter logisticGrowthParameter = null;
    Parameter logisticShapeParameter = null;
    Parameter exponentialGrowthParameter = null;
    Parameter transistionTimeParameter = null;
    double alpha = 0.5;
    ExponentialLogistic exponentialLogistic = null;
}
/*
 * TemperOperator.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.inference.operators.shrinkage;

import dr.inference.model.Parameter;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.xml.*;

import static dr.inferencexml.operators.shrinkage.BayesianBridgeShrinkageOperatorParser.BAYESIAN_BRIDGE_PARSER;

/**
 * @author Marc A. Suchard
 * @author Alexander Fisher
 */

public class TemperOperator extends SimpleMCMCOperator implements GibbsOperator {

    public static final String TEMPER_OPERATOR = "temperOperator";
    private static final String TARGET_PARAMETER = "target";
    private static final String RATE = "rate";

    private final Parameter parameter;
    private final Parameter target;
    private final double rate;
    private double currentCount;
    private final double[] startValues;

    public TemperOperator(Parameter parameter,
                          Parameter target,
                          double rate,
                          double weight) {

        this.parameter = parameter;
        this.target = target;
        this.rate = rate;

        startValues = parameter.getParameterValues();
        currentCount = 0;
        setWeight(weight);
    }

    @Override
    public String getOperatorName() {
        return BAYESIAN_BRIDGE_PARSER;
    }

    @Override
    public double doOperation() {
        // plan: start + (end - start) * percentComplete
        // where percentComplete = [ exp(-Nsteps * exp(-rate)) - 1 ]  is in (0, 1), and importantly goes from *0* to *1*
        // Note we achieve percentComplete = 1 only when Nsteps * exp(-rate) := 'objective' is large

        // factor out a -1 to read code below:
        // start + (start - end) * -1*percentComplete

        currentCount = currentCount + 1;

        // seems unnecessary, could just input exp(-rate) directly as "rate" in XML if desired.
        double objective = currentCount * Math.exp(-rate);

        for (int i = 0; i < parameter.getDimension(); ++i) {

            double currentValue = parameter.getParameterValue(i);
            double targetValue = target.getParameterValue(i % target.getDimension());

            // if less than (approx) 90% of the way there, keep going: exp(-2.3) ~ 0.1
            if (objective < 2.3) {
                double startValue = startValues[i];
                double x = startValue + (startValue - targetValue) * (Math.exp(-objective) - 1);
                parameter.setParameterValue(i, x);

            } else if (currentValue != targetValue) {
                parameter.setParameterValue(i, targetValue);
            }
        }

//        System.out.println("Exponent value: " + parameter.getValue(0));
        return 0;
    }

    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return TEMPER_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double weight = xo.getDoubleAttribute(SimpleMCMCOperator.WEIGHT);

            Parameter hotParameter = (Parameter) xo.getChild(Parameter.class);

            XMLObject cxo = xo.getChild(TARGET_PARAMETER);
            Parameter targetParameter = (cxo.getChild(0) instanceof Parameter) ?
                (Parameter) cxo.getChild(Parameter.class) :
                new Parameter.Default(cxo.getDoubleChild(0));

            if (targetParameter.getDimension() > hotParameter.getDimension()) {
                throw new XMLParseException("Target parameter cannot have more dimensions than the tempered parameter.");
            }

            double rate = xo.getAttribute(RATE, 0.0);
            if (rate < 0.0) {
                throw new XMLParseException("Rate cannot be negative");
            }

            return new TemperOperator(hotParameter, targetParameter, rate, weight);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        @Override
        public String getParserDescription() {
            return "Tempers (anneals) a parameter from starting value to target value.";
        }

        @SuppressWarnings("rawtypes")
        @Override
        public Class getReturnType() {
            return TemperOperator.class;
        }

        private final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(Parameter.class),
                new ElementRule(TARGET_PARAMETER, new XMLSyntaxRule[]{
                        new XORRule(
                                new ElementRule(Parameter.class),
                                new ElementRule(Double.class)
                        )}),
                AttributeRule.newDoubleRule(RATE),
        };
    };
}

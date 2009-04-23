/*
 * ScaleOperator.java
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

package dr.inference.operators;

import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.xml.*;

/**
 * A generic scale operator for use with a multi-dimensional parameters.
 * Either scale all dimentions at once or scale one dimention at a time.
 *
 * @author Joseph Heled
 * @version $Id$
 */
public class LogRandomWalkOperator extends AbstractCoercableOperator {

    public static final String LOGRANDOMWALK_OPERATOR = "logRandomWalkOperator";
    // Use same attributes as scale operator to help users
    public static final String SCALE_ALL = ScaleOperator.SCALE_ALL;
    public static final String SCALE_ALL_IND = ScaleOperator.SCALE_ALL_IND;
    public static final String WINDOW_SIZE = "window";

    private double size;
    private Parameter parameter = null;
    private final boolean scaleAll;
    private final boolean scaleAllInd;


    public LogRandomWalkOperator(Parameter parameter, double size,
                                 CoercionMode mode, double weight, boolean scaleAll, boolean scaleAllInd) {

        super(mode);
        setWeight(weight);

        this.parameter = parameter;
        assert parameter.getDimension() == 1;
        this.size = size;

        this.scaleAll = scaleAll;
        this.scaleAllInd = scaleAllInd;
    }


    private double scaleOne(int j) {
        final double w = size * (MathUtils.nextDouble() - 0.5);
        final double newValue = parameter.getParameterValue(j) * Math.exp(w);
        parameter.setParameterValue(j, newValue);

        return w;
    }

    /**
     * change the parameter and return the hastings ratio.
     */
    public final double doOperation() throws OperatorFailedException {
        final int dim = parameter.getDimension();
        final Bounds bounds = parameter.getBounds();
        double hastingsRatio = 0;

        int checkStart = 0, checkEnd = dim;

        // Must first set all parameters first and check for boundaries later for the operator to work
        // correctly with dependent parameters such as tree node heights.

        if( scaleAllInd ) {
            for(int i = 0; i < dim; i++) {
                hastingsRatio += scaleOne(i);
            }
        } else if( scaleAll ) {
            final double w = size * (MathUtils.nextDouble() - 0.5);
            final double f = Math.exp(w);
            for(int i = 0; i < dim; i++) {
                parameter.setParameterValue(i, parameter.getParameterValue(i) * f);
            }
            hastingsRatio += dim * w;
        } else {
            int j = MathUtils.nextInt(dim);
            hastingsRatio += scaleOne(j);
            checkStart = j;
            checkEnd = j + 1;
        }

        for(int i = checkStart; i < checkEnd; i++) {
            final double value = parameter.getParameterValue(i);
            if( value < bounds.getLowerLimit(i) || value > bounds.getUpperLimit(i) ) {
                throw new OperatorFailedException("proposed value outside boundaries");
            }
        }
        return hastingsRatio;
    }

    //MCMCOperator INTERFACE
    public final String getOperatorName() {
        return "logRandomWalk" + (scaleAllInd ? "-all" : "") +
                (scaleAllInd ? "-independently" : "") +
                "(" + parameter.getParameterName() + ")";
    }

    public double getCoercableParameter() {
        return size;
    }

    public void setCoercableParameter(double value) {
        if( getMode() != CoercionMode.COERCION_OFF ) {
            size = value;
        }
    }

    public double getRawParameter() {
        return size;
    }

    public double getTargetAcceptanceProbability() {
        return 0.234;
    }

    public final String getPerformanceSuggestion() {

        double prob = Utils.getAcceptanceProbability(this);
        double targetProb = getTargetAcceptanceProbability();
        dr.util.NumberFormatter formatter = new dr.util.NumberFormatter(5);
        double sf = size;
        if( prob < getMinimumGoodAcceptanceLevel() ) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else if( prob > getMaximumGoodAcceptanceLevel() ) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else {
            return "";
        }
    }


    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return LOGRANDOMWALK_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            CoercionMode mode = CoercionMode.parseMode(xo);

            final double weight = xo.getDoubleAttribute(WEIGHT);
            final double size = xo.getDoubleAttribute(WINDOW_SIZE);
            final boolean scaleAll = xo.getAttribute(SCALE_ALL, false);
            final boolean scaleAllInd = xo.getAttribute(SCALE_ALL_IND, false);

            if( size <= 0.0 ) {
                throw new XMLParseException("size must be positive");
            }

            final Parameter parameter = (Parameter) xo.getChild(Parameter.class);


            return new LogRandomWalkOperator(parameter, size, mode, weight, scaleAll, scaleAllInd);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a scale operator on a given parameter.";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(WINDOW_SIZE),
                AttributeRule.newDoubleRule(WEIGHT),
                AttributeRule.newBooleanRule(AUTO_OPTIMIZE, true),

                AttributeRule.newBooleanRule(SCALE_ALL, true),
                AttributeRule.newBooleanRule(SCALE_ALL_IND, true),

                new ElementRule(Parameter.class),
        };

    };

    public String toString() {
        return getOperatorName();
    }
}
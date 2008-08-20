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
 * An optional bit vector and a threshold is used to vary the rate of the individual dimentions according
 * to their on/off status. For example a threshold of 1 means pick only "on" dimentions.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: ScaleOperator.java,v 1.20 2005/06/14 10:40:34 rambaut Exp $
 */
public class ScaleOperator extends AbstractCoercableOperator {

    public static final String SCALE_OPERATOR = "scaleOperator";
    public static final String SCALE_ALL = "scaleAll";
    public static final String SCALE_ALL_IND = "scaleAllIndependently";
    public static final String SCALE_FACTOR = "scaleFactor";
    public static final String DEGREES_OF_FREEDOM = "df";
    public static final String INDICATORS = "indicators";
    public static final String PICKONEPROB = "pickoneprob";

    private Parameter indicator;
    private double indicatorOnProb;

    public ScaleOperator(Parameter parameter, double scale, CoercionMode mode, double weight) {

        this(parameter, false, 0, scale, mode, null, 1.0, false);
        setWeight(weight);
    }

    public ScaleOperator(Parameter parameter, boolean scaleAll, int degreesOfFreedom, double scale,
                         CoercionMode mode, Parameter indicator, double indicatorOnProb, boolean scaleAllInd) {

        super(mode);

        this.parameter = parameter;
        this.indicator = indicator;
        this.indicatorOnProb = indicatorOnProb;
        this.scaleAll = scaleAll;
        this.scaleAllIndependently = scaleAllInd;
        this.scaleFactor = scale;
        this.degreesOfFreedom = degreesOfFreedom;
    }


    /**
     * @return the parameter this operator acts on.
     */
    public Parameter getParameter() {
        return parameter;
    }

    /**
     * change the parameter and return the hastings ratio.
     */
    public final double doOperation() throws OperatorFailedException {

        final double scale = (scaleFactor + (MathUtils.nextDouble() * ((1.0 / scaleFactor) - scaleFactor)));

        double logq;

        final Bounds bounds = parameter.getBounds();
        final int dim = parameter.getDimension();

        if (scaleAllIndependently) {
            // update all dimensions independently.
            logq = 0;
            for (int i = 0; i < dim; i++) {

                final double scaleOne = (scaleFactor + (MathUtils.nextDouble() * ((1.0 / scaleFactor) - scaleFactor)));
                final double value = scaleOne * parameter.getParameterValue(i);

                logq -= Math.log(scaleOne);

                if (value < bounds.getLowerLimit(i) || value > bounds.getUpperLimit(i)) {
                    throw new OperatorFailedException("proposed value outside boundaries");
                }

                parameter.setParameterValue(i, value);

            }
        } else if (scaleAll) {
            // update all dimensions
            // hasting ratio is dim-2 times of 1dim case. would be nice to have a reference here
            // for the proof. It is supposed to be somewhere in an Alexei/Nicholes article.
            if (degreesOfFreedom > 0)
                // For parameters with non-uniform prior on only one dimension
                logq = -degreesOfFreedom * Math.log(scale);
            else
                logq = (dim - 2) * Math.log(scale);

            for (int i = 0; i < dim; i++) {
                final double value = parameter.getParameterValue(i) * scale;
                if (value < bounds.getLowerLimit(i) || value > bounds.getUpperLimit(i)) {
                    throw new OperatorFailedException("proposed value outside boundaries");
                }
                parameter.setParameterValue(i, value);
            }
        } else {
            logq = -Math.log(scale);

            // which bit to scale
            int index;
            if (indicator != null) {
                final int idim = indicator.getDimension();
                final boolean impliedOne = idim == (dim - 1);
                // available bit locations
                int[] loc = new int[idim + 1];
                int nLoc = 0;
                // choose active or non active ones?
                final boolean takeOne = indicatorOnProb >= 1.0 || MathUtils.nextDouble() < indicatorOnProb;

                if (impliedOne && takeOne) {
                    loc[nLoc] = 0;
                    ++nLoc;
                }
                for (int i = 0; i < idim; i++) {
                    final double value = indicator.getStatisticValue(i);
                    if (takeOne == (value > 0)) {
                        loc[nLoc] = i + (impliedOne ? 1 : 0);
                        ++nLoc;
                    }
                }

                if (nLoc > 0) {
                    final int rand = MathUtils.nextInt(nLoc);
                    index = loc[rand];
                } else {
                    throw new OperatorFailedException("no active indicators");
                }
            } else {
                // any is good
                index = MathUtils.nextInt(dim);
            }

            final double oldValue = parameter.getParameterValue(index);
            final double newValue = scale * oldValue;

            if (newValue < bounds.getLowerLimit(index) || newValue > bounds.getUpperLimit(index)) {
                throw new OperatorFailedException("proposed value outside boundaries");
            }

            parameter.setParameterValue(index, newValue);

            // provides a hook for subclasses
            cleanupOperation(newValue, oldValue);
        }

        return logq;
    }

    /**
     * This method should be overridden by operators that need to do something just before the return of doOperation.
     *
     * @param newValue the proposed parameter value
     * @param oldValue the old parameter value
     */
    void cleanupOperation(double newValue, double oldValue) {
        // DO NOTHING
    }

    //MCMCOperator INTERFACE
    public final String getOperatorName() {
        return "scale(" + parameter.getParameterName() + ")";
    }

    public double getCoercableParameter() {
        return Math.log(1.0 / scaleFactor - 1.0);
    }

    public void setCoercableParameter(double value) {
        scaleFactor = 1.0 / (Math.exp(value) + 1.0);
    }

    public double getRawParameter() {
        return scaleFactor;
    }

    public double getScaleFactor() {
        return scaleFactor;
    }

    public double getTargetAcceptanceProbability() {
        return 0.234;
    }

    public final String getPerformanceSuggestion() {

        double prob = MCMCOperator.Utils.getAcceptanceProbability(this);
        double targetProb = getTargetAcceptanceProbability();
        dr.util.NumberFormatter formatter = new dr.util.NumberFormatter(5);
        double sf = OperatorUtils.optimizeScaleFactor(scaleFactor, prob, targetProb);
        if (prob < getMinimumGoodAcceptanceLevel()) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else if (prob > getMaximumGoodAcceptanceLevel()) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else return "";
    }


    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        public String getParserName() {
            return SCALE_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            boolean scaleAll = xo.getAttribute(SCALE_ALL, false);
            boolean scaleAllInd = xo.getAttribute(SCALE_ALL_IND, false);
            int degreesOfFreedom = xo.getAttribute(DEGREES_OF_FREEDOM, 0);

            CoercionMode mode = CoercionMode.parseMode(xo);

            final double weight = xo.getDoubleAttribute(WEIGHT);
            final double scaleFactor = xo.getDoubleAttribute(SCALE_FACTOR);

            if (scaleFactor <= 0.0 || scaleFactor >= 1.0) {
                throw new XMLParseException("scaleFactor must be between 0.0 and 1.0");
            }

            final Parameter parameter = (Parameter) xo.getChild(Parameter.class);

            Parameter indicator = null;
            double indicatorOnProb = 1.0;
            final XMLObject cxo = (XMLObject) xo.getChild(INDICATORS);

            if (cxo != null) {
                indicator = (Parameter) cxo.getChild(Parameter.class);
                if (cxo.hasAttribute(PICKONEPROB)) {
                    indicatorOnProb = cxo.getDoubleAttribute(PICKONEPROB);
                    if (!(0 <= indicatorOnProb && indicatorOnProb <= 1)) {
                        throw new XMLParseException("pickoneprob must be between 0.0 and 1.0");
                    }
                }
            }
            ScaleOperator operator = new ScaleOperator(parameter, scaleAll,
                    degreesOfFreedom, scaleFactor,
                    mode, indicator, indicatorOnProb,
                    scaleAllInd);
            operator.setWeight(weight);
            return operator;
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

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                AttributeRule.newDoubleRule(SCALE_FACTOR),
                AttributeRule.newBooleanRule(SCALE_ALL, true),
                AttributeRule.newBooleanRule(SCALE_ALL_IND, true),
                AttributeRule.newDoubleRule(WEIGHT),
                AttributeRule.newBooleanRule(AUTO_OPTIMIZE, true),

                new ElementRule(Parameter.class),
                new ElementRule(INDICATORS,
                        new XMLSyntaxRule[]{
                                AttributeRule.newDoubleRule(PICKONEPROB, true),
                                new ElementRule(Parameter.class)}, true),
        };

    };

    public String toString() {
        return "scaleOperator(" + parameter.getParameterName() + " [" + scaleFactor + ", " + (1.0 / scaleFactor) + "]";
    }

    //PRIVATE STUFF

    private Parameter parameter = null;
    private boolean scaleAll = false;
    private boolean scaleAllIndependently = false;
    private int degreesOfFreedom = 0;
    private double scaleFactor = 0.5;
}

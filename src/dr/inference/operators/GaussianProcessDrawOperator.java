/*
 * GaussianProcessDrawOperator.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.inference.model.Parameter;
import dr.inferencexml.operators.EllipticalSliceOperatorParser;
import dr.math.distributions.GaussianProcessRandomGenerator;
import dr.xml.*;


/**
 * @author Marc Suchard
 */
public class GaussianProcessDrawOperator extends AbstractCoercableOperator {

    public static final String GAUSSIAN_PROCESS_OPERATOR = "gaussianProcessOperator";
    public static final String SCALE_FACTOR = "scaleFactor";
    public static final String TRANSLATION_INVARIANT = EllipticalSliceOperatorParser.TRANSLATION_INVARIANT;
    public static final String ROTATION_INVARIANT = EllipticalSliceOperatorParser.ROTATION_INVARIANT;

    private double scaleFactor;
    private final Parameter parameter;
    private final GaussianProcessRandomGenerator gaussianProcess;
    private final int dim;
    private final boolean translationInvariant;
    private final boolean rotationInvariant;

    public GaussianProcessDrawOperator(Parameter parameter, double scaleFactor, double weight,
                                       CoercionMode mode, GaussianProcessRandomGenerator gaussianProcess,
                                       boolean translationInvariant, boolean rotationInvariant) {

        super(mode);
        this.scaleFactor = scaleFactor;
        this.parameter = parameter;
        this.gaussianProcess = gaussianProcess;
        setWeight(weight);

        this.translationInvariant = translationInvariant;
        this.rotationInvariant = rotationInvariant;

        if (gaussianProcess.getDimension() != parameter.getDimension()) {
            throw new IllegalArgumentException("Dimension mismatch");
        }

//        if (!gaussianProcess.isTranslationInvariant()) {
//            throw new IllegalArgumentException("Must be translationally invariant");
//        }

        this.dim = parameter.getDimension();

    }

    public double doOperation() {

        double[] x = parameter.getParameterValues();
        double[] epsilon = (double[]) gaussianProcess.nextRandom();

        for (int i = 0; i < dim; ++i) {
            x[i] += scaleFactor * epsilon[i];
        }

        EllipticalSliceOperator.transformPoint(x, translationInvariant, rotationInvariant,
                2); // TODO Make generic for dim != 2

        for (int i = 0; i < dim; i++) {
            parameter.setParameterValueQuietly(i, x[i] + scaleFactor * epsilon[i]);
        }
        parameter.fireParameterChangedEvent();

//        System.err.println("DONE");

        return 0; // TODO Not 0 if gaussianProcess mean != 0
    }

    //MCMCOperator INTERFACE
    public final String getOperatorName() {
        return parameter.getParameterName();
    }

    public double getCoercableParameter() {
        return Math.log(scaleFactor);
    }

    public void setCoercableParameter(double value) {
        scaleFactor = Math.exp(value);
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

    public double getMinimumAcceptanceLevel() {
        return 0.1;
    }

    public double getMaximumAcceptanceLevel() {
        return 0.4;
    }

    public double getMinimumGoodAcceptanceLevel() {
        return 0.20;
    }

    public double getMaximumGoodAcceptanceLevel() {
        return 0.30;
    }

    public final String getPerformanceSuggestion() {

        double prob = Utils.getAcceptanceProbability(this);
        double targetProb = getTargetAcceptanceProbability();
        dr.util.NumberFormatter formatter = new dr.util.NumberFormatter(5);
        double sf = OperatorUtils.optimizeWindowSize(scaleFactor, prob, targetProb);
        if (prob < getMinimumGoodAcceptanceLevel()) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else if (prob > getMaximumGoodAcceptanceLevel()) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else return "";
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return GAUSSIAN_PROCESS_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            CoercionMode mode = CoercionMode.parseMode(xo);

            double weight = xo.getDoubleAttribute(WEIGHT);
            double scaleFactor = xo.getDoubleAttribute(SCALE_FACTOR);

            if (scaleFactor <= 0.0) {
                throw new XMLParseException("scaleFactor must be greater than 0.0");
            }

            Parameter parameter = (Parameter) xo.getChild(Parameter.class);
            GaussianProcessRandomGenerator generator =
                    (GaussianProcessRandomGenerator) xo.getChild(GaussianProcessRandomGenerator.class);

            boolean translationInvariant = xo.getAttribute(TRANSLATION_INVARIANT, false);
            boolean rotationInvariant = xo.getAttribute(ROTATION_INVARIANT, false);

            return new GaussianProcessDrawOperator(parameter, scaleFactor, weight, mode, generator,
                    translationInvariant, rotationInvariant);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a multivariate normal random walk operator on a given parameter.";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(SCALE_FACTOR),
                AttributeRule.newDoubleRule(WEIGHT),
                AttributeRule.newBooleanRule(AUTO_OPTIMIZE, true),
                AttributeRule.newBooleanRule(TRANSLATION_INVARIANT, true),
                AttributeRule.newBooleanRule(ROTATION_INVARIANT, true),
                new ElementRule(Parameter.class),
                new ElementRule(GaussianProcessRandomGenerator.class),
        };

    };
}

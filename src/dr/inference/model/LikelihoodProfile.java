/*
 * LikelihoodProfile.java
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

package dr.inference.model;

import dr.math.NumericalDerivative;
import dr.math.UnivariateFunction;
import dr.math.distributions.RandomGenerator;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 */

public class LikelihoodProfile implements Reportable {

    public static final String SEP = ",";
    public static final String ENDL = "\n";
    public static final String LIKELIHOOD_PROFILE = "likelihoodProfile";

    private RandomGenerator generator;
    private int mcSamples = 1;
    private Parameter integrated;

    public LikelihoodProfile(Likelihood likelihood_, Parameter parameter_, final int dim_,
                             double lowerBound_, double upperBound_, int numPoints_) {
        this.likelihood = likelihood_;
        this.parameter = parameter_;
        this.dim = dim_;
        this.lowerBound = lowerBound_;
        this.upperBound = upperBound_;
        this.numPoints = numPoints_;

        function = new UnivariateFunction() {

            public double evaluate(double argument) {
                parameter.setParameterValue(dim, argument);
                return likelihood.getLogLikelihood();
            }

            public double getLowerBound() {
                return lowerBound;
            }

            public double getUpperBound() {
                return upperBound;
            }
        };
    }

    public String getReport() {
        if (profilePoints.size() == 0) {
            generateProfile();
        }
        StringBuilder sb = new StringBuilder();
        sb.append("x").append(SEP)
                .append("f").append(SEP)
                .append("df").append(SEP)
                .append("ddf").append(ENDL);
        for (ProfilePoint pt : profilePoints) {
            sb.append(pt.toString());
        }
        return sb.toString();
    }

    private int count = 0;

    private void generateProfile() {
        double delta = (upperBound - lowerBound) / (numPoints - 1);
        double x = lowerBound;
        for (int i = 0; i < numPoints; ++i) {
//            System.err.println("i = " + i);
            double fx = 0, dx = 0, ddx = 0;
            for (int j = 0; j < mcSamples; ++j) {
//                System.err.println("j = " + j);
                if (generator != null) {

//                    System.err.println("Original: ");
//                    System.err.println(new Vector(integrated.getParameterValues()));

                    integrateGenerator(x);
//                    System.err.println(new Vector(integrated.getParameterValues()));
//                    System.err.println("");
//
//                    if (count > 0) {
//                        System.exit(-1);
//                    }
//                    count++;

                }
                fx += function.evaluate(x);
                dx += NumericalDerivative.firstDerivative(function, x);
                ddx += NumericalDerivative.secondDerivative(function, x);
            }
            fx /= mcSamples;
            dx /= mcSamples;
            ddx /= mcSamples;
            profilePoints.add(new ProfilePoint(x, fx, dx, ddx));
            x += delta;
        }
    }

    private void integrateGenerator(double argument) {
        parameter.setParameterValue(dim, argument);
        double[] draw = (double[]) generator.nextRandom();
        int dim = draw.length;

        if (dim != integrated.getDimension()) {
            throw new RuntimeException("Invalid integrated parameter and generator");
        }
        for (int i = 0; i < dim; ++i) {
            draw[i] = 0;
            integrated.setParameterValue(i, draw[i]);
        }
    }

    private class ProfilePoint {
        public double x;
        public double f;
        public double df;
        public double ddf;

        ProfilePoint(double x, double f, double df, double ddf) {
            this.x = x;
            this.f = f;
            this.df = df;
            this.ddf = ddf;
        }

        public String toString() {
            return new StringBuilder()
                    .append(x).append(SEP)
                    .append(f).append(SEP)
                    .append(df).append(SEP)
                    .append(ddf).append(ENDL)
                    .toString();
        }
    }

    private final Likelihood likelihood;
    private final Parameter parameter;
    private final UnivariateFunction function;
    private int dim;
    private double lowerBound;
    private double upperBound;
    private int numPoints;
    private List<ProfilePoint> profilePoints = new ArrayList<ProfilePoint>();

    public static final String DIM = "dim";
    public static final String LOWER_BOUND = "lower";
    public static final String UPPER_BOUND = "upper";
    public static final String GRID_POINTS = "points";
    public static final String EXPECTATION = "expectation";
    public static final String MC_SAMPLES = "samples";

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return LIKELIHOOD_PROFILE;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Likelihood likelihood = (Likelihood) xo.getChild(Likelihood.class);
            Parameter parameter = (Parameter) xo.getChild(Parameter.class);
            int dim = xo.getAttribute(DIM, 0);
            double lowerBound = xo.getDoubleAttribute(LOWER_BOUND);
            double upperBound = xo.getDoubleAttribute(UPPER_BOUND);
            int numPoints = xo.getAttribute(GRID_POINTS, 100);

            LikelihoodProfile profile = new LikelihoodProfile(likelihood, parameter, dim, lowerBound, upperBound, numPoints);


            if (xo.hasChildNamed(EXPECTATION)) {
//                System.err.println("Here");
                XMLObject cxo = xo.getChild(EXPECTATION);
                int mcSamples = cxo.getAttribute(MC_SAMPLES, 1);
                RandomGenerator trait =
                        (RandomGenerator) cxo.getChild(RandomGenerator.class);

                Parameter integrated = (Parameter) cxo.getChild(Parameter.class);
                profile.addGeneratorForExpectation(trait, integrated, mcSamples);
//                System.err.println("trait null ? " + (trait == null ? "yes" : "no"));
//                System.exit(-1);
            }

            return profile;
        }

        public String getParserDescription() {
            return "This element represents a tool to profile a likelihood surface";
        }

        public Class getReturnType() {
            return LikelihoodProfile.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleArrayRule(LOWER_BOUND),
                AttributeRule.newDoubleArrayRule(UPPER_BOUND),
                AttributeRule.newIntegerRule(DIM, true),
                AttributeRule.newIntegerRule(GRID_POINTS, true),
                new ElementRule(Likelihood.class),
                new ElementRule(Parameter.class),
                new ElementRule(EXPECTATION,
                        new XMLSyntaxRule[]{
                                AttributeRule.newIntegerRule(MC_SAMPLES),
                                new ElementRule(RandomGenerator.class),
                                new ElementRule(Parameter.class),
                        }, true),
        };
    };

    private void addGeneratorForExpectation(RandomGenerator generator, Parameter integrated, int mcSamples) {
        this.generator = generator;
        this.integrated = integrated;
        this.mcSamples = mcSamples;
    }
}

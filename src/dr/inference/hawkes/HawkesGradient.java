/*
 * HawkesGradient.java
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

package dr.inference.hawkes;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.HessianWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inferencexml.operators.hmc.HamiltonianMonteCarloOperatorParser;
import dr.xml.*;

/**
 * @author Andrew Holbrook
 * @author Xiang Ji
 * @author Marc Suchard
 */
public class HawkesGradient implements GradientWrtParameterProvider, HessianWrtParameterProvider, Reportable {

    final private WrtParameter wrtParameter;
    final private HawkesLikelihood likelihood;
    final private boolean withHessian;
    final private Double tolerance;

    public HawkesGradient(WrtParameter wrtParameter,
                          HawkesLikelihood likelihood,
                          Double tolerance,
                          boolean withHessian) {
        this.wrtParameter = wrtParameter;
        this.likelihood = likelihood;
        this.tolerance = tolerance;
        this.withHessian = withHessian;
    }

    @Override
    public Likelihood getLikelihood() {
        return likelihood;
    }

    @Override
    public Parameter getParameter() {
        return wrtParameter.getParameter(likelihood);
    }

    @Override
    public int getDimension() {
        return wrtParameter.getParameter(likelihood).getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        return wrtParameter.getGradientLogDensity(likelihood);
    }

    @Override
    public String getReport() {
        return GradientWrtParameterProvider.getReportAndCheckForError(this, 0, Double.POSITIVE_INFINITY, tolerance)
                + (withHessian ? HessianWrtParameterProvider.getReportAndCheckForError(this, tolerance) : "");
    }

    @Override
    public double[] getDiagonalHessianLogDensity() {
        return wrtParameter.getDiagonalHessianLogDensity(likelihood);
    }

    @Override
    public double[][] getHessianLogDensity() {
        throw new RuntimeException("Not yet implemented!");
    }

    enum WrtParameter {

        LOCATIONS("locations") {
            @Override
            Parameter getParameter(HawkesLikelihood likelihood) {
                return likelihood.getHawkesModel().getLocationsParameter();
            }

            @Override
            double[] getGradientLogDensity(HawkesLikelihood likelihood) {
                return likelihood.getGradientLogDensity();
            }

            @Override
            double[] getDiagonalHessianLogDensity(HawkesLikelihood likelihood) {
                throw new RuntimeException("Not yet implemented!");
            }
        },

        RANDOM_RATES("randomRates") {
            @Override
            Parameter getParameter(HawkesLikelihood likelihood) {
                return likelihood.getHawkesModel().getRateProvider().getParameter();
            }

            @Override
            double[] getGradientLogDensity(HawkesLikelihood likelihood) {
                return likelihood.getHawkesModel().getRateProvider().orderByNodeIndex(likelihood.getRandomRateGradient());
            }

            @Override
            double[] getDiagonalHessianLogDensity(HawkesLikelihood likelihood) {
                return likelihood.getHawkesModel().getRateProvider().orderByNodeIndex(likelihood.getRandomRateHessian());
            }
        };
        WrtParameter(String name) {
            this.name = name;
        }

        abstract Parameter getParameter(HawkesLikelihood likelihood);

        abstract double[] getGradientLogDensity(HawkesLikelihood likelihood);

        abstract double[] getDiagonalHessianLogDensity(HawkesLikelihood likelihood);

        private final String name;

        public static HawkesGradient.WrtParameter factory(String match) {
            for (HawkesGradient.WrtParameter type : HawkesGradient.WrtParameter.values()) {
                if (match.equalsIgnoreCase(type.name)) {
                    return type;
                }
            }
            return null;
        }
    }

    // **************************************************************
    // XMLObjectParser
    // **************************************************************

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        final static String HAWKES_GRADIENT = "hawkesGradient";
        final static String WRT_PARAMETER = "wrt";
        final static String HESSIAN = "hessian";
        private static final String TOLERANCE = HamiltonianMonteCarloOperatorParser.GRADIENT_CHECK_TOLERANCE;

        public String getParserName() {
            return HAWKES_GRADIENT;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            HawkesLikelihood likelihood = (HawkesLikelihood) xo.getChild(HawkesLikelihood.class);
            String wrtParameter = (String) xo.getAttribute(WRT_PARAMETER);

            WrtParameter wrt = WrtParameter.factory(wrtParameter);

            boolean withHessian = xo.getAttribute(HESSIAN, false);
            double tolerance = xo.getAttribute(TOLERANCE, 1E-4);

            return new HawkesGradient(wrt, likelihood, tolerance, withHessian);

        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "Gradient w.r.t. parameters of HawkesLikelihood.";
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newStringRule(WRT_PARAMETER),
                new ElementRule(HawkesLikelihood.class),
                AttributeRule.newBooleanRule(HESSIAN, true),
                AttributeRule.newDoubleRule(TOLERANCE, true),
        };

        public Class getReturnType() {
            return HawkesGradient.class;
        }
    };
}

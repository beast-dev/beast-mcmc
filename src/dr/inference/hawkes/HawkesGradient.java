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
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * @author Andrew Holbrook
 * @author Xiang Ji
 * @author Marc Suchard
 */
public class HawkesGradient implements GradientWrtParameterProvider, Reportable {

    final private WrtParameter wrtParameter;
    final private HawkesLikelihood likelihood;

    public HawkesGradient(WrtParameter wrtParameter,
                          HawkesLikelihood likelihood) {
        this.wrtParameter = wrtParameter;
        this.likelihood = likelihood;
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
        return GradientWrtParameterProvider.getReportAndCheckForError(this, 0, Double.POSITIVE_INFINITY, null);
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
        },

        RANDOM_RATES("randomRates") {
            @Override
            Parameter getParameter(HawkesLikelihood likelihood) {
                return likelihood.getHawkesModel().getRateProvider().getParameter();
            }

            @Override
            double[] getGradientLogDensity(HawkesLikelihood likelihood) {
                return likelihood.getRandomRateGradient();
            }
        };
        WrtParameter(String name) {
            this.name = name;
        }

        abstract Parameter getParameter(HawkesLikelihood likelihood);

        abstract double[] getGradientLogDensity(HawkesLikelihood likelihood);

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
        final static String WRT_PARAMETER = "wrtParameter";

        public String getParserName() {
            return HAWKES_GRADIENT;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            HawkesLikelihood likelihood = (HawkesLikelihood) xo.getChild(HawkesLikelihood.class);
            String wrtParameter = (String) xo.getAttribute(WRT_PARAMETER);

            WrtParameter wrt = WrtParameter.factory(wrtParameter);

            return new HawkesGradient(wrt, likelihood);

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
                new ElementRule(HawkesLikelihood.class)
        };

        public Class getReturnType() {
            return HawkesGradient.class;
        }
    };
}

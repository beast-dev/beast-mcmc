/*
 * ConstantLikelihood.java
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

import dr.xml.*;

/**
 * A class that returns a constant log likelihood.
 * Can be used for debugging purposes when arbitrary normalizing constants are required
 *
 * @author Marc Suchard
 */
public class ConstantLikelihood extends Likelihood.Abstract {

    public static final String CONSTANT_LIKELIHOOD = "constantLikelihood";
    public static final String LOG_VALUE = "logValue";

    public ConstantLikelihood(double logValue) {
        super(null);
        this.logValue = logValue;
    }

    protected boolean getLikelihoodKnown() {
        return false;
    }

    public double calculateLogLikelihood() {
        return logValue;
    }

    public boolean evaluateEarly() {
        return true;
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public static final String DATA = "data";

        public String getParserName() {
            return CONSTANT_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double logValue = xo.getAttribute(LOG_VALUE, 0.0);
            return new ConstantLikelihood(logValue);
        }

        public String getParserDescription() {
            return "A function that returns a constant value as a likelihood.";
        }

        public Class getReturnType() {
            return ConstantLikelihood.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(LOG_VALUE),
        };
    };

    private final double logValue;
}


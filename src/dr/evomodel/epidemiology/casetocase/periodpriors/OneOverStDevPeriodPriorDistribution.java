/*
 * OneOverStDevPeriodPriorDistribution.java
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

package dr.evomodel.epidemiology.casetocase.periodpriors;

import dr.xml.*;
import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;

/**
 * This is effectively the Jeffreys prior on the set of periods assuming that they are normally (or lognormally)
 * distributed; the probability is proportional to the reciprocal of their standard deviation (or the standard deviation
 * of their logarithms).
 */
public class OneOverStDevPeriodPriorDistribution extends AbstractPeriodPriorDistribution {

    public static final String ONE_OVER_STDEV = "oneOverStDevPeriodPriorDistribution";
    public static final String LOG = "log";
    public static final String ID = "id";

    public OneOverStDevPeriodPriorDistribution(String name, boolean log){
        super(name, log);
    }

    @Override
    public void reset() {

    }

    @Override
    public double calculateLogPosteriorProbability(double newValue, double minValue) {
        return 0;
    }

    @Override
    public double calculateLogPosteriorCDF(double limit, boolean upper) {
        return 0;
    }

    public double calculateLogLikelihood(double[] values){

        DescriptiveStatistics stats = new DescriptiveStatistics(values);

        logL = -Math.log(stats.getStandardDeviation());

        return logL;

    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return ONE_OVER_STDEV;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            String id = (String) xo.getAttribute(ID);

            boolean log;
            log = xo.hasAttribute(LOG) ? xo.getBooleanAttribute(LOG) : false;

            return new OneOverStDevPeriodPriorDistribution(id, log);

        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newBooleanRule(LOG, true),
                AttributeRule.newStringRule(ID, false)
        };

        public String getParserDescription() {
            return "Calculates the probability of observing a list of doubles with probability proportional to" +
                    "1 over their standard deviation (the Jeffreys prior for normally distributed data)";
        }

        public Class getReturnType() {
            return OneOverStDevPeriodPriorDistribution.class;
        }
    };

}

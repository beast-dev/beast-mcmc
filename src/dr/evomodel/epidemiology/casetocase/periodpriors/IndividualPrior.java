/*
 * IndividualPrior.java
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

import dr.evomodel.epidemiology.casetocase.AbstractCase;
import dr.inference.distribution.LogNormalDistributionModel;
import dr.inference.distribution.NormalDistributionModel;
import dr.inference.distribution.ParametricDistributionModel;
import dr.math.distributions.LogNormalDistribution;
import dr.math.distributions.NormalDistribution;
import dr.xml.*;

import java.util.HashMap;

/**
 * Created by mhall on 03/06/2014.
 */
public class IndividualPrior extends AbstractPeriodPriorDistribution {

    ParametricDistributionModel distribution;

    public static final String INDIVIDUAL_PRIOR = "individualPrior";
    public static final String ID = "id";
    public static final String DISTRIBUTION = "distribution";

    public IndividualPrior(String name, ParametricDistributionModel distribution){
        super(name, false);
        this.distribution = distribution;
    }

    public void reset() {
        // nothing to do.
    }

    public double calculateLogPosteriorProbability(double newValue, double minValue) {
        double logNumerator = distribution.logPdf(newValue);
        double logDenominator = Math.log(1-distribution.cdf(minValue));

        if(logDenominator == Double.NEGATIVE_INFINITY){
            if(distribution instanceof LogNormalDistributionModel){
                double mean = ((LogNormalDistributionModel)distribution).getMu();
                double stdev = ((LogNormalDistributionModel)distribution).getSigma();

                double scaledValue = (Math.log(minValue)-mean)/stdev;

                logDenominator = NormalDistribution.standardCDF(-scaledValue, true);

            } else if (distribution instanceof NormalDistributionModel){
                double mean = ((NormalDistributionModel)distribution).getMean().getValue(0);
                double stdev = ((NormalDistributionModel)distribution).getStdev();

                double scaledValue = (minValue-mean)/stdev;

                logDenominator = NormalDistribution.standardCDF(-scaledValue, true);
            }
        }
        return logNumerator - logDenominator;
    }

    public double calculateLogPosteriorCDF(double limit, boolean upper) {
        double out;
        if(upper) {
            out = Math.log(1-distribution.cdf(limit));
        } else {
            out = Math.log(distribution.cdf(limit));
        }
        if(out == Double.NEGATIVE_INFINITY){
            if(distribution instanceof LogNormalDistributionModel){
                double mean = ((LogNormalDistributionModel)distribution).getMu();
                double stdev = ((LogNormalDistributionModel)distribution).getSigma();

                double scaledValue = (Math.log(limit)-mean)/stdev;

                out = NormalDistribution.standardCDF(-scaledValue, true);

            } else if (distribution instanceof NormalDistributionModel){
                double mean = ((NormalDistributionModel)distribution).getMean().getValue(0);
                double stdev = ((NormalDistributionModel)distribution).getStdev();

                double scaledValue = (limit-mean)/stdev;

                out = NormalDistribution.standardCDF(-scaledValue, true);
            }
        }
        return out;
    }

    public double calculateLogLikelihood(double[] values) {
        double out = 0;

        for (double value : values) {
            out += distribution.logPdf(value);
        }

        return out;
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return INDIVIDUAL_PRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            String id = (String) xo.getAttribute(ID);
            ParametricDistributionModel distribution =
                    (ParametricDistributionModel)xo.getElementFirstChild(DISTRIBUTION);

            return new IndividualPrior(id, distribution);

        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newStringRule(ID, false),
                new ElementRule(DISTRIBUTION, ParametricDistributionModel.class)
        };

        public String getParserDescription() {
            return "Calculates the probability of a set of doubles all being drawn from the specified prior " +
                    "distribution";
        }

        public Class getReturnType() {
            return IndividualPrior.class;
        }
    };
}

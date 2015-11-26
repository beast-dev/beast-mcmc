/*
 * TruncatedWorkingDistribution.java
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

package dr.evomodel.continuous;

import dr.inference.distribution.AbstractDistributionLikelihood;
import dr.inference.distribution.DistributionLikelihood;
import dr.inference.model.Statistic;
import dr.math.distributions.NormalDistribution;
import dr.xml.*;

/**
 * @author Marc A. Suchard
 */
public class TruncatedWorkingDistribution extends AbstractDistributionLikelihood {

    public static final String MEAN = "mean";
    public static final String STDEV = "stdev";
    public static final String WORKING_PRIOR = "truncatedWorkingPrior";

    private DistributionLikelihood baseDistribution;
    private LatentTruncation truncation;

    public TruncatedWorkingDistribution(DistributionLikelihood baseDistribution,
                                        LatentTruncation truncation) {
        super(baseDistribution.getModel());
        this.baseDistribution = baseDistribution;
        this.truncation = truncation;
    }

    @Override
    public double calculateLogLikelihood() {
        double logLikelihood = truncation.getLogLikelihood(); // Returns log(0) or log(1)
        if (logLikelihood != Double.NEGATIVE_INFINITY) {
            logLikelihood += baseDistribution.getLogLikelihood() +
                    truncation.getNormalizationConstant(baseDistribution.getDistribution());
        }
        return logLikelihood;
    }

    public static XMLObjectParser WORKING_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return WORKING_PRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double mean = xo.getDoubleAttribute(MEAN);
            double stdev = xo.getDoubleAttribute(STDEV);

            DistributionLikelihood likelihood = new DistributionLikelihood(new NormalDistribution(mean, stdev));
            for (int j = 0; j < xo.getChildCount(); j++) {
                if (xo.getChild(j) instanceof Statistic) {
                    likelihood.addData((Statistic) xo.getChild(j));
                } else if (!(xo.getChild(j) instanceof LatentTruncation)) {
                    throw new XMLParseException("Illegal element in " + xo.getName() + " element");
                }
            }

            LatentTruncation truncation = (LatentTruncation) xo.getChild(LatentTruncation.class);

            return new TruncatedWorkingDistribution(likelihood, truncation);
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(MEAN),
                AttributeRule.newDoubleRule(STDEV),
                new ElementRule(Statistic.class, 1, Integer.MAX_VALUE),
                new ElementRule(LatentTruncation.class),
        };

        public String getParserDescription() {
            return "Calculates the prior probability of some data under a latent truncated normal distribution.";
        }

        public Class getReturnType() {
            return TruncatedWorkingDistribution.class;
        }
    };
}

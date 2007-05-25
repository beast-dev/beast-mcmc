/*
 * DistributionLikelihood.java
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

package dr.inference.distribution;

import dr.inference.model.Likelihood;
import dr.inference.model.Statistic;
import dr.math.*;
import dr.xml.*;
import org.w3c.dom.Document;
import org.w3c.dom.Element;

/**
 * A class that returns the log likelihood of a set of data (statistics)
 * being distributed according to the given parametric distribution.
 * @author Alexei Drummond
 * @version $Id: DistributionLikelihood.java,v 1.11 2005/05/25 09:35:28 rambaut Exp $
 */

public class DistributionLikelihood extends AbstractDistributionLikelihood {

	public static final String DISTRIBUTION_LIKELIHOOD = "distributionLikelihood";

	public static final String DISTRIBUTION = "distribution";
	public static final String DATA = "data";

	public DistributionLikelihood(Distribution distribution) {
		super(null);
        this.distribution = distribution;
        this.offset = 0.0;
	}

    public DistributionLikelihood(Distribution distribution, double offset) {
        super(null);
        this.distribution = distribution;
        this.offset = offset;
    }

    public DistributionLikelihood(ParametricDistributionModel distributionModel) {
        super(distributionModel);
        this.distribution = distributionModel;
        this.offset = 0.0;
	}

	// **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

	/**
     * Calculate the log likelihood of the current state.
     * @return the log likelihood.
     */
	public double calculateLogLikelihood() {

		double logL = 0.0;

		for (int i = 0; i < dataList.size(); i++) {
			Statistic statistic = (Statistic)dataList.get(i);
			for (int j = 0; j < statistic.getDimension(); j++) {

                double value = statistic.getStatisticValue(j) - offset;
                if (offset > 0.0 && value < 0.0) {
                    // fixes a problem with the offset on exponential distributions not
                    // actually bounding the distribution. This only performs this check
                    // if a non-zero offset is actually given otherwise it assumes the
                    // parameter is either legitimately allowed to go negative or is bounded
                    // at zero anyway.
                    return Double.NEGATIVE_INFINITY;
                }

                logL += distribution.logPdf(value);
			}
		}
		return logL;
	}

	// **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

	public Element createElement(Document d) {
		throw new RuntimeException("Not implemented yet!");
	}


	/**
	 * Reads a distribution likelihood from a DOM Document element.
	 */
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() { return DISTRIBUTION_LIKELIHOOD; }

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			XMLObject cxo = (XMLObject)xo.getChild(DISTRIBUTION);
			ParametricDistributionModel model = (ParametricDistributionModel)cxo.getChild(ParametricDistributionModel.class);

			DistributionLikelihood likelihood = new DistributionLikelihood(model);

			cxo = (XMLObject)xo.getChild(DATA);

				for (int j = 0; j < cxo.getChildCount(); j++) {
					if (cxo.getChild(j) instanceof Statistic) {

                        likelihood.addData( (Statistic)cxo.getChild(j));
					} else {
						throw new XMLParseException("illegal element in " + cxo.getName() + " element");
					}
				}

			return likelihood;
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			new ElementRule(DISTRIBUTION,
				new XMLSyntaxRule[] { new ElementRule(ParametricDistributionModel.class) }),
			new ElementRule(DATA, new XMLSyntaxRule[] {new ElementRule(Statistic.class, 1, Integer.MAX_VALUE )})
		};

		public String getParserDescription() {
			return "Calculates the likelihood of some data given some parametric or empirical distribution.";
		}

		public Class getReturnType() { return Likelihood.class; }
	};

    public static final String UNIFORM_PRIOR = "uniformPrior";
    public static final String EXPONENTIAL_PRIOR = "exponentialPrior";
    public static final String NORMAL_PRIOR = "normalPrior";
    public static final String LOG_NORMAL_PRIOR = "logNormalPrior";
    public static final String GAMMA_PRIOR = "gammaPrior";
    public static final String UPPER = "upper";
    public static final String LOWER = "lower";
    public static final String MEAN = "mean";
	public static final String MEAN_IN_REAL_SPACE = "meanInRealSpace";
    public static final String STDEV = "stdev";
    public static final String SHAPE = "shape";
    public static final String SCALE = "scale";
    public static final String OFFSET = "offset";

    /**
     * A special parser that reads a convenient short form of priors on parameters.
     */
    public static XMLObjectParser UNIFORM_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() { return UNIFORM_PRIOR; }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double lower = xo.getDoubleAttribute(LOWER);
            double upper = xo.getDoubleAttribute(UPPER);

            DistributionLikelihood likelihood = new DistributionLikelihood(new UniformDistribution(lower, upper));
            for (int j = 0; j < xo.getChildCount(); j++) {
                if (xo.getChild(j) instanceof Statistic) {
                    likelihood.addData( (Statistic)xo.getChild(j));
                } else {
                    throw new XMLParseException("illegal element in " + xo.getName() + " element");
                }
            }

            return likelihood;
        }

        public XMLSyntaxRule[] getSyntaxRules() { return rules; }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            AttributeRule.newDoubleRule(LOWER),
            AttributeRule.newDoubleRule(UPPER),
            new ElementRule(Statistic.class, 1, Integer.MAX_VALUE )
        };

        public String getParserDescription() {
            return "Calculates the prior probability of some data under a given uniform distribution.";
		}

		public Class getReturnType() { return Likelihood.class; }
    };

    /**
     * A special parser that reads a convenient short form of priors on parameters.
     */
    public static XMLObjectParser EXPONENTIAL_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() { return EXPONENTIAL_PRIOR; }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double mean = xo.getDoubleAttribute(MEAN);
            double offset = xo.getDoubleAttribute(OFFSET);

            DistributionLikelihood likelihood = new DistributionLikelihood(new ExponentialDistribution(1.0/mean), offset);
            for (int j = 0; j < xo.getChildCount(); j++) {
                if (xo.getChild(j) instanceof Statistic) {
                    likelihood.addData( (Statistic)xo.getChild(j));
                } else {
                    throw new XMLParseException("illegal element in " + xo.getName() + " element");
                }
            }

            return likelihood;
        }

        public XMLSyntaxRule[] getSyntaxRules() { return rules; }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            AttributeRule.newDoubleRule(MEAN),
            AttributeRule.newDoubleRule(OFFSET),
            new ElementRule(Statistic.class, 1, Integer.MAX_VALUE )
	};

        public String getParserDescription() {
            return "Calculates the prior probability of some data under a given exponential distribution.";
        }

        public Class getReturnType() { return Likelihood.class; }
    };

    /**
     * A special parser that reads a convenient short form of priors on parameters.
     */
    public static XMLObjectParser NORMAL_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() { return NORMAL_PRIOR; }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double mean = xo.getDoubleAttribute(MEAN);
            double stdev = xo.getDoubleAttribute(STDEV);

            DistributionLikelihood likelihood = new DistributionLikelihood(new NormalDistribution(mean, stdev));
            for (int j = 0; j < xo.getChildCount(); j++) {
                if (xo.getChild(j) instanceof Statistic) {
                    likelihood.addData( (Statistic)xo.getChild(j));
                } else {
                    throw new XMLParseException("illegal element in " + xo.getName() + " element");
                }
            }

            return likelihood;
        }

        public XMLSyntaxRule[] getSyntaxRules() { return rules; }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            AttributeRule.newDoubleRule(MEAN),
            AttributeRule.newDoubleRule(STDEV),
            new ElementRule(Statistic.class, 1, Integer.MAX_VALUE )
        };

        public String getParserDescription() {
            return "Calculates the prior probability of some data under a given normal distribution.";
        }

        public Class getReturnType() { return Likelihood.class; }
    };

    /**
     * A special parser that reads a convenient short form of priors on parameters.
     */
    public static XMLObjectParser LOG_NORMAL_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() { return LOG_NORMAL_PRIOR; }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double mean = xo.getDoubleAttribute(MEAN);
            double stdev = xo.getDoubleAttribute(STDEV);
	        double offset = xo.getDoubleAttribute(OFFSET);
	        boolean meanInRealSpace = xo.getBooleanAttribute(MEAN_IN_REAL_SPACE);

	        if (meanInRealSpace) {
		        throw new UnsupportedOperationException("meanInRealSpace is not supported yet");
	        }

            DistributionLikelihood likelihood = new DistributionLikelihood(new LogNormalDistribution(mean, stdev), offset);
            for (int j = 0; j < xo.getChildCount(); j++) {
                if (xo.getChild(j) instanceof Statistic) {
                    likelihood.addData( (Statistic)xo.getChild(j));
                } else {
                    throw new XMLParseException("illegal element in " + xo.getName() + " element");
                }
            }

            return likelihood;
        }

        public XMLSyntaxRule[] getSyntaxRules() { return rules; }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            AttributeRule.newDoubleRule(MEAN),
            AttributeRule.newDoubleRule(STDEV),
            AttributeRule.newDoubleRule(OFFSET),
		        AttributeRule.newBooleanRule(MEAN_IN_REAL_SPACE),
            new ElementRule(Statistic.class, 1, Integer.MAX_VALUE )
        };

        public String getParserDescription() {
            return "Calculates the prior probability of some data under a given lognormal distribution.";
        }

        public Class getReturnType() { return Likelihood.class; }
    };


    /**
     * A special parser that reads a convenient short form of priors on parameters.
     */
    public static XMLObjectParser GAMMA_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() { return GAMMA_PRIOR; }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double shape = xo.getDoubleAttribute(SHAPE);
            double scale = xo.getDoubleAttribute(SCALE);
            double offset = xo.getDoubleAttribute(OFFSET);

            DistributionLikelihood likelihood = new DistributionLikelihood(new GammaDistribution(shape, scale), offset);
            for (int j = 0; j < xo.getChildCount(); j++) {
                if (xo.getChild(j) instanceof Statistic) {
                    likelihood.addData( (Statistic)xo.getChild(j));
                } else {
                    throw new XMLParseException("illegal element in " + xo.getName() + " element");
                }
            }

            return likelihood;
        }

        public XMLSyntaxRule[] getSyntaxRules() { return rules; }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            AttributeRule.newDoubleRule(SHAPE),
            AttributeRule.newDoubleRule(SCALE),
            AttributeRule.newDoubleRule(OFFSET),
            new ElementRule(Statistic.class, 1, Integer.MAX_VALUE )
        };

        public String getParserDescription() {
            return "Calculates the prior probability of some data under a given gamma distribution.";
        }

        public Class getReturnType() { return Likelihood.class; }
    };

    private final Distribution distribution;
    private final double offset;
}


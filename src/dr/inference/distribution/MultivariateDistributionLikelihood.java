/*
 * MultivariateDistributionLikelihood.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.inference.distribution;

import dr.inference.model.*;
import dr.math.distributions.*;
import dr.xml.*;
import dr.util.Attribute;


/**
 * @author Marc Suchard
 */
public class MultivariateDistributionLikelihood extends AbstractDistributionLikelihood {

    public static final String MVN_PRIOR = "multivariateNormalPrior";
    public static final String MVN_MEAN = "meanParameter";
    public static final String MVN_PRECISION = "precisionParameter";
    public static final String MVN_CV = "coefficientOfVariation";
    public static final String WISHART_PRIOR = "multivariateWishartPrior";
    public static final String INV_WISHART_PRIOR = "multivariateInverseWishartPrior";
    public static final String DIRICHLET_PRIOR = "dirichletPrior";
    public static final String DF = "df";
    public static final String SCALE_MATRIX = "scaleMatrix";
    public static final String MVGAMMA_PRIOR = "multivariateGammaPrior";
    public static final String MVGAMMA_SHAPE = "shapeParameter";
    public static final String MVGAMMA_SCALE = "scaleParameter";
    public static final String COUNTS = "countsParameter";
    public static final String NON_INFORMATIVE = "nonInformative";

    public static final String DATA = "data";

    private final MultivariateDistribution distribution;

    public MultivariateDistributionLikelihood(MultivariateDistribution distribution) {
        super(new DefaultModel());
        this.distribution = distribution;
    }

    public double calculateLogLikelihood() {
        double logL = 0.0;

        for( Attribute<double[]> data : dataList ) {
            logL += distribution.logPdf(data.getAttributeValue());
        }
        return logL;
    }

    public MultivariateDistribution getDistribution() {
        return distribution;
    }

    public static XMLObjectParser DIRICHLET_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return DIRICHLET_PRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {


            XMLObject cxo = xo.getChild(COUNTS);
            Parameter counts = (Parameter) cxo.getChild(Parameter.class);

            DirichletDistribution dirichlet = new DirichletDistribution(counts.getParameterValues());

            MultivariateDistributionLikelihood likelihood = new MultivariateDistributionLikelihood(
                    dirichlet);

            cxo = xo.getChild(DATA);
            for (int j = 0; j < cxo.getChildCount(); j++) {
                if (cxo.getChild(j) instanceof Parameter) {
                    likelihood.addData((Parameter) cxo.getChild(j));
                } else {
                    throw new XMLParseException("illegal element in " + xo.getName() + " element " + cxo.getName());
                }
            }


            return likelihood;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(COUNTS,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
        };

        public String getParserDescription() {
            return "Calculates the likelihood of some data under an Inverse-Wishart distribution.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }
    };

    public static XMLObjectParser INV_WISHART_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return INV_WISHART_PRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            int df = xo.getIntegerAttribute(DF);

            XMLObject cxo = xo.getChild(SCALE_MATRIX);
            MatrixParameter scaleMatrix = (MatrixParameter) cxo.getChild(MatrixParameter.class);
            InverseWishartDistribution invWishart = new InverseWishartDistribution(df, scaleMatrix.getParameterAsMatrix());

            MultivariateDistributionLikelihood likelihood = new MultivariateDistributionLikelihood(
                    invWishart);

            cxo = xo.getChild(DATA);
            for (int j = 0; j < cxo.getChildCount(); j++) {
                if (cxo.getChild(j) instanceof MatrixParameter) {
                    likelihood.addData((MatrixParameter) cxo.getChild(j));
                } else {
                    throw new XMLParseException("illegal element in " + xo.getName() + " element " + cxo.getName());
                }
            }


            return likelihood;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(DF),
                new ElementRule(SCALE_MATRIX,
                        new XMLSyntaxRule[]{new ElementRule(MatrixParameter.class)}),
        };

        public String getParserDescription() {
            return "Calculates the likelihood of some data under an Inverse-Wishart distribution.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }
    };

    public static XMLObjectParser WISHART_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return WISHART_PRIOR;
        }


        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            MultivariateDistributionLikelihood likelihood;

            if (xo.hasAttribute(NON_INFORMATIVE) && xo.getBooleanAttribute(NON_INFORMATIVE)) {
                // Make non-informative settings
                XMLObject cxo = xo.getChild(DATA);
                int dim = ((MatrixParameter) cxo.getChild(0)).getColumnDimension();
                likelihood = new MultivariateDistributionLikelihood(new WishartDistribution(dim));
            } else {
                if (!xo.hasAttribute(DF) || !xo.hasChildNamed(SCALE_MATRIX)) {
                    throw new XMLParseException("Must specify both a df and scaleMatrix");
                }

                int df = xo.getIntegerAttribute(DF);

                XMLObject cxo = xo.getChild(SCALE_MATRIX);
                MatrixParameter scaleMatrix = (MatrixParameter) cxo.getChild(MatrixParameter.class);

                likelihood = new MultivariateDistributionLikelihood(
                        new WishartDistribution(df, scaleMatrix.getParameterAsMatrix())
                );

            }

            XMLObject cxo = xo.getChild(DATA);
            for (int j = 0; j < cxo.getChildCount(); j++) {
                if (cxo.getChild(j) instanceof MatrixParameter) {
                    likelihood.addData((MatrixParameter) cxo.getChild(j));
                    System.err.println("added ");
                } else {
                    throw new XMLParseException("illegal element in " + xo.getName() + " element " + cxo.getName());
                }
            }

            return likelihood;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules;{
            rules = new XMLSyntaxRule[]{
                    AttributeRule.newBooleanRule(NON_INFORMATIVE, true),
                    AttributeRule.newIntegerRule(DF, true),
                    new ElementRule(SCALE_MATRIX,
                            new XMLSyntaxRule[]{new ElementRule(MatrixParameter.class)}, true),
                    new ElementRule(DATA,
                            new XMLSyntaxRule[]{new ElementRule(MatrixParameter.class, 1, Integer.MAX_VALUE)}
                    )
            };
        }

        public String getParserDescription() {
            return "Calculates the likelihood of some data under a Wishart distribution.";
        }

        public Class getReturnType() {
            return Likelihood.class;
        }
    };

    public static XMLObjectParser MVN_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return MVN_PRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            XMLObject cxo = xo.getChild(MVN_MEAN);
            Parameter mean = (Parameter) cxo.getChild(Parameter.class);

            cxo = xo.getChild(MVN_PRECISION);
            MatrixParameter precision = (MatrixParameter) cxo.getChild(MatrixParameter.class);

            if (mean.getDimension() != precision.getRowDimension() ||
                    mean.getDimension() != precision.getColumnDimension())
                throw new XMLParseException("Mean and precision have wrong dimensions in " + xo.getName() + " element");

            MultivariateDistributionLikelihood likelihood =
                    new MultivariateDistributionLikelihood(
                            new MultivariateNormalDistribution(mean.getParameterValues(),
                                    precision.getParameterAsMatrix())
                    );
            cxo = xo.getChild(DATA);
            if (cxo != null) {
                for (int j = 0; j < cxo.getChildCount(); j++) {
                    if (cxo.getChild(j) instanceof Parameter) {
                        Parameter data = (Parameter) cxo.getChild(j);
                        likelihood.addData(data);
                        if (data.getDimension() != mean.getDimension())
                            throw new XMLParseException("dim(" + data.getStatisticName() + ") != dim(" + mean.getStatisticName() + ") in " + xo.getName() + "element");
                    } else {
                        throw new XMLParseException("illegal element in " + xo.getName() + " element");
                    }
                }
            }

            return likelihood;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                new ElementRule(MVN_MEAN,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                new ElementRule(MVN_PRECISION,
                        new XMLSyntaxRule[]{new ElementRule(MatrixParameter.class)}),
                new ElementRule(DATA,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class, 1, Integer.MAX_VALUE)}, true)
        };

        public String getParserDescription() {
            return "Calculates the likelihood of some data under a given multivariate-normal distribution.";
        }

        public Class getReturnType() {
            return MultivariateDistributionLikelihood.class;
        }
    };

    public static XMLObjectParser MVGAMMA_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return MVGAMMA_PRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {


            double[] shape;
            double[] scale;


            if (xo.hasChildNamed(MVGAMMA_SHAPE)) {

                XMLObject cxo = xo.getChild(MVGAMMA_SHAPE);
                shape = ((Parameter) cxo.getChild(Parameter.class)).getParameterValues();

                cxo = xo.getChild(MVGAMMA_SCALE);
                scale = ((Parameter) cxo.getChild(Parameter.class)).getParameterValues();

                if (shape.length != scale.length)
                    throw new XMLParseException("Shape and scale have wrong dimensions in " + xo.getName() + " element");

            } else {

                XMLObject cxo = xo.getChild(MVN_MEAN);
                double[] mean = ((Parameter) cxo.getChild(Parameter.class)).getParameterValues();

                cxo = xo.getChild(MVN_CV);
                double[] cv = ((Parameter) cxo.getChild(Parameter.class)).getParameterValues();

                if (mean.length != cv.length)
                    throw new XMLParseException("Mean and CV have wrong dimensions in " + xo.getName() + " element");

                final int dim = mean.length;
                shape = new double[dim];
                scale = new double[dim];

                for (int i = 0; i < dim; i++) {
                    double c2 = cv[i] * cv[i];
                    shape[i] = 1.0 / c2;
                    scale[i] = c2 * mean[i];
                }
            }

            MultivariateDistributionLikelihood likelihood =
                    new MultivariateDistributionLikelihood(
                            new MultivariateGammaDistribution(shape, scale)
                    );
            XMLObject cxo = xo.getChild(DATA);
            for (int j = 0; j < cxo.getChildCount(); j++) {
                if (cxo.getChild(j) instanceof Parameter) {
                    Parameter data = (Parameter) cxo.getChild(j);
                    likelihood.addData(data);
                    if (data.getDimension() != shape.length)
                        throw new XMLParseException("dim(" + data.getStatisticName() + ") != " + shape.length + " in " + xo.getName() + "element");
                } else {
                    throw new XMLParseException("illegal element in " + xo.getName() + " element");
                }
            }
            return likelihood;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {

                new XORRule(
                        new ElementRule(MVGAMMA_SHAPE,
                                new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                        new ElementRule(MVN_MEAN,
                                new XMLSyntaxRule[]{new ElementRule(Parameter.class)})),
                new ElementRule(DATA,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class, 1, Integer.MAX_VALUE)})
        };

        public String getParserDescription() {
            return "Calculates the likelihood of some data under a given multivariate-gamma distribution.";
        }

        public Class getReturnType() {
            return MultivariateDistributionLikelihood.class;
        }
    };
}

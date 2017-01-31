/*
 * MultivariateDistributionLikelihood.java
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

package dr.inference.distribution;

import dr.evomodel.continuous.FullyConjugateMultivariateTraitLikelihood;
import dr.evomodel.continuous.TreeTraitNormalDistributionModel;
import dr.inference.model.*;
import dr.inferencexml.distribution.DistributionLikelihoodParser;
import dr.math.distributions.*;
import dr.util.Attribute;
import dr.util.Citable;
import dr.util.Transform;
import dr.xml.*;

import java.util.logging.Logger;


/**
 * @author Marc Suchard
 * @author Guy Baele
 */
public class MultivariateDistributionLikelihood extends AbstractDistributionLikelihood {

    public static final String MVN_PRIOR = "multivariateNormalPrior";
    public static final String MVN_MEAN = "meanParameter";
    public static final String MVN_PRECISION = "precisionParameter";
    public static final String MVN_CV = "coefficientOfVariation";
    public static final String WISHART_PRIOR = "multivariateWishartPrior";
    public static final String INV_WISHART_PRIOR = "multivariateInverseWishartPrior";
    public static final String DIRICHLET_PRIOR = "dirichletParameterPrior";
    public static final String SUM_TO_NUMBER_OF_ELEMENTS = "sumToNumberOfElements";
    public static final String DF = "df";
    public static final String SCALE_MATRIX = "scaleMatrix";
    public static final String MVGAMMA_PRIOR = "multivariateGammaPrior";
    public static final String MVGAMMA_SHAPE = "shapeParameter";
    public static final String MVGAMMA_SCALE = "scaleParameter";
    public static final String COUNTS = "countsParameter";
    public static final String NON_INFORMATIVE = "nonInformative";
    public static final String MULTIVARIATE_LIKELIHOOD = "multivariateDistributionLikelihood";
    public static final String DATA_AS_MATRIX = "dataAsMatrix";
    // public static final String TREE_TRAIT = "treeTraitNormalDistribution";
    public static final String TREE_TRAIT = "treeTraitNormalDistributionLikelihood";
    public static final String TREE_TRAIT_NORMAL = "treeTraitNormalDistribution";
    public static final String ROOT_VALUE = "rootValue";
    public static final String CONDITION = "conditionOnRoot";

    public static final String DATA = "data";

    private final MultivariateDistribution distribution;
    private final Transform[] transforms;
    private Parameter parameter = null;

    public MultivariateDistributionLikelihood(String name, ParametricMultivariateDistributionModel model) {
        this(name, model, null);
    }

    public MultivariateDistributionLikelihood(String name, ParametricMultivariateDistributionModel model,
                                              Transform[] transforms) {
        super(model);
        this.distribution = model;
        this.transforms = transforms;
    }

    public MultivariateDistributionLikelihood(String name, MultivariateDistribution distribution) {
        this(name, distribution, null);
    }

    public MultivariateDistributionLikelihood(String name, MultivariateDistribution distribution,
                                              Transform[] transforms) {
        super(new DefaultModel(name));
        this.distribution = distribution;
        this.transforms = transforms;
    }

    public MultivariateDistributionLikelihood(MultivariateDistribution distribution) {
        this(distribution, null);
    }

    public MultivariateDistributionLikelihood(MultivariateDistribution distribution, Transform[] transforms) {
        this(distribution.getType(), distribution, transforms);
    }

    public String toString() {
        return getClass().getName() + "(" + getLogLikelihood() + ")";
    }

    public double calculateLogLikelihood() {
        double logL = 0.0;

        for (Attribute<double[]> data : dataList) {
            double[] x =
//                    (data instanceof Parameter) ?
//                            ((Parameter) data).getParameterValues() :
                            data.getAttributeValue();
            if (transforms != null) {
                double[] y = new double[x.length];
                for (int i = 0; i < x.length; ++i) {
                    logL += transforms[i].getLogJacobian(x[i]);
                    y[i] = transforms[i].transform(x[i]);
                }
                logL += distribution.logPdf(y);
            } else {
                logL += distribution.logPdf(x);
            }
        }
        return logL;
    }

    public void addData(Parameter parameter) {
        this.parameter = parameter;

        addData((Attribute<double[]>)parameter);
    }

    public Parameter getDataParameter() {
        return parameter;
    }

    @Override
    public void addData(Attribute<double[]> data) {
        super.addData(data);

        if (data instanceof Variable && getModel() instanceof DefaultModel) {
            ((DefaultModel) getModel()).addVariable((Variable) data);
        }
    }

    public MultivariateDistribution getDistribution() {
        return distribution;
    }

    public static Transform[] parseListOfTransforms(XMLObject xo, int maxDim) throws XMLParseException {
        Transform[] transforms = null;

        boolean anyTransforms = false;
        for (int i = 0; i < xo.getChildCount(); ++i) {
            if (xo.getChild(i) instanceof Transform.ParsedTransform) {
                Transform.ParsedTransform t = (Transform.ParsedTransform) xo.getChild(i);
                if (transforms == null) {
                    transforms = Transform.Util.getListOfNoTransforms(maxDim);
                }

                t.end = Math.max(t.end, maxDim);
                if (t.start < 0 || t.end < 0 || t.start > t.end) {
                    throw new XMLParseException("Invalid bounds for transform in " + xo.getId());
                }
                for (int j = t.start; j < t.end; j += t.every) {
                    transforms[j] = t.transform;
                    anyTransforms = true;
                }
            }
        }
        if (anyTransforms) {
            StringBuilder sb = new StringBuilder("Using distributional transforms in " + xo.getId() + "\n");
            for (int i = 0; i < transforms.length; ++i) {
                if (transforms[i] != Transform.NONE) {
                    sb.append("\t").append(transforms[i].getTransformName()).append(" on index ")
                            .append(i + 1).append("\n");
                }
            }
            Logger.getLogger("dr.utils.Transform").info(sb.toString());
        }
        return transforms;
    }

    public static XMLObjectParser DIRICHLET_PRIOR_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return DIRICHLET_PRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            boolean sumConstraint = false;
            if (xo.hasAttribute(SUM_TO_NUMBER_OF_ELEMENTS)) {
                sumConstraint = xo.getBooleanAttribute(SUM_TO_NUMBER_OF_ELEMENTS);
            }

            XMLObject cxo = xo.getChild(COUNTS);
            Parameter counts = (Parameter) cxo.getChild(Parameter.class);

            DirichletDistribution dirichlet = new DirichletDistribution(counts.getParameterValues(), sumConstraint);

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
                AttributeRule.newBooleanRule(SUM_TO_NUMBER_OF_ELEMENTS, true),
                new ElementRule(COUNTS,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                new ElementRule(DATA,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class)}, 1, Integer.MAX_VALUE),
        };

        public String getParserDescription() {
            return "Calculates the likelihood of some data under a Dirichlet distribution.";
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

                double df = xo.getDoubleAttribute(DF);

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
                    AttributeRule.newDoubleRule(DF, true),
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

    public static XMLObjectParser MULTIVARIATE_LIKELIHOOD_PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return MULTIVARIATE_LIKELIHOOD;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            XMLObject cxo = xo.getChild(DistributionLikelihoodParser.DISTRIBUTION);
            ParametricMultivariateDistributionModel distribution = (ParametricMultivariateDistributionModel)
                    cxo.getChild(ParametricMultivariateDistributionModel.class);

            // Parse transforms here
            int maxDim = distribution.getMean().length;
            Transform[] transforms = parseListOfTransforms(xo, maxDim);

            MultivariateDistributionLikelihood likelihood = new MultivariateDistributionLikelihood(xo.getId(),
                    distribution, transforms);

            boolean dataAsMatrix = xo.getAttribute(DATA_AS_MATRIX, false);

            cxo = xo.getChild(DATA);
            if (cxo != null) {
                for (int j = 0; j < cxo.getChildCount(); j++) {
                    if (cxo.getChild(j) instanceof Parameter) {
                        Parameter data = (Parameter) cxo.getChild(j);
                        if (data instanceof MatrixParameter) {
                            MatrixParameter matrix = (MatrixParameter) data;
                            if (dataAsMatrix) {
                                likelihood.addData(matrix);
                            } else {
                                if (matrix.getParameter(0).getDimension() != distribution.getMean().length)
                                    throw new XMLParseException("dim(" + data.getStatisticName() + ") = " + matrix.getParameter(0).getDimension()
                                            + " is not equal to dim(" + distribution.getType() + ") = " + distribution.getMean().length
                                            + " in " + xo.getName() + "element");

                                for (int i = 0; i < matrix.getParameterCount(); i++) {
                                    likelihood.addData(matrix.getParameter(i));
                                }
                            }
                        } else {
                            if (data.getDimension() != distribution.getMean().length)
                                throw new XMLParseException("dim(" + data.getStatisticName() + ") = " + data.getDimension()
                                        + " is not equal to dim(" + distribution.getType() + ") = " + distribution.getMean().length
                                        + " in " + xo.getName() + "element");
                            likelihood.addData(data);
                        }
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
                new ElementRule(DistributionLikelihoodParser.DISTRIBUTION,
                        new XMLSyntaxRule[]{new ElementRule(ParametricMultivariateDistributionModel.class)}
                ),
                AttributeRule.newBooleanRule(DATA_AS_MATRIX, true),
                new ElementRule(Transform.ParsedTransform.class, 0, Integer.MAX_VALUE),
                new ElementRule(DATA,
                        new XMLSyntaxRule[]{new ElementRule(Parameter.class, 1, Integer.MAX_VALUE)}, true)
        };

        public String getParserDescription() {
            return "Calculates the likelihood of some data under a given multivariate distribution.";
        }

        public Class getReturnType() {
            return MultivariateDistributionLikelihood.class;
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

            Transform[] transforms = parseListOfTransforms(xo, mean.getDimension());

            MultivariateDistributionLikelihood likelihood =
                    new MultivariateDistributionLikelihood(
                            new MultivariateNormalDistribution(mean.getParameterValues(),
                                    precision.getParameterAsMatrix()), transforms
                    );
            cxo = xo.getChild(DATA);
            if (cxo != null) {
                for (int j = 0; j < cxo.getChildCount(); j++) {
                    if (cxo.getChild(j) instanceof Parameter) {
                        Parameter data = (Parameter) cxo.getChild(j);
                        if (data instanceof MatrixParameter) {
                            MatrixParameter matrix = (MatrixParameter) data;
                            if (matrix.getParameter(0).getDimension() != mean.getDimension())
                                throw new XMLParseException("dim(" + data.getStatisticName() + ") = " + matrix.getParameter(0).getDimension()
                                        + " is not equal to dim(" + mean.getStatisticName() + ") = " + mean.getDimension()
                                        + " in " + xo.getName() + "element");

                            for (int i = 0; i < matrix.getParameterCount(); i++) {
                                likelihood.addData(matrix.getParameter(i));
                            }
                        } else {
                            if (data.getDimension() != mean.getDimension())
                                throw new XMLParseException("dim(" + data.getStatisticName() + ") = " + data.getDimension()
                                        + " is not equal to dim(" + mean.getStatisticName() + ") = " + mean.getDimension()
                                        + " in " + xo.getName() + "element");
                            likelihood.addData(data);
                        }
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
                new ElementRule(Transform.ParsedTransform.class, 0, Integer.MAX_VALUE),
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
                new XORRule(
                        new ElementRule(MVGAMMA_SCALE,
                                new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
                        new ElementRule(MVN_CV,
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

    public static XMLObjectParser TREE_TRAIT_MODEL = new AbstractXMLObjectParser() {

        public String getParserName() {
            return TREE_TRAIT_NORMAL;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            boolean conditionOnRoot = xo.getAttribute(CONDITION, false);

            FullyConjugateMultivariateTraitLikelihood traitModel = (FullyConjugateMultivariateTraitLikelihood)
                    xo.getChild(FullyConjugateMultivariateTraitLikelihood.class);

            TreeTraitNormalDistributionModel treeTraitModel;

            if (xo.getChild(ROOT_VALUE) != null) {
                XMLObject cxo = xo.getChild(ROOT_VALUE);
                Parameter rootValue = (Parameter) cxo.getChild(Parameter.class);
                treeTraitModel = new TreeTraitNormalDistributionModel(traitModel, rootValue, conditionOnRoot);
            } else {
                treeTraitModel = new TreeTraitNormalDistributionModel(traitModel, conditionOnRoot);
            }
            return treeTraitModel;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newBooleanRule(CONDITION, true),
                new ElementRule(FullyConjugateMultivariateTraitLikelihood.class)
        };

        public String getParserDescription() {
            return "Parses TreeTraitNormalDistributionModel";
        }

        public Class getReturnType() {
            return TreeTraitNormalDistributionModel.class;
        }
    };


    public static XMLObjectParser TREE_TRAIT_DISTRIBUTION = new AbstractXMLObjectParser() {

        public String getParserName() {
            return TREE_TRAIT;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            /*
            boolean conditionOnRoot = xo.getAttribute(CONDITION, false);

            FullyConjugateMultivariateTraitLikelihood traitModel = (FullyConjugateMultivariateTraitLikelihood)
                    xo.getChild(FullyConjugateMultivariateTraitLikelihood.class);
            */

            TreeTraitNormalDistributionModel treeTraitModel = (TreeTraitNormalDistributionModel)
                    xo.getChild(TreeTraitNormalDistributionModel.class);

            MultivariateDistributionLikelihood likelihood =
                    new MultivariateDistributionLikelihood(
                            //  new TreeTraitNormalDistributionModel(traitModel, conditionOnRoot)
                            treeTraitModel
                    );

            XMLObject cxo = xo.getChild(DATA);
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
                //  AttributeRule.newBooleanRule(CONDITION, true),
                //   new ElementRule(FullyConjugateMultivariateTraitLikelihood.class),
                new ElementRule(TreeTraitNormalDistributionModel.class),
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

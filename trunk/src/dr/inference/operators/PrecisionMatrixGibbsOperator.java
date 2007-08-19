/*
 * PrecisionMatrixGibbsOperator.java
 *
 * Copyright (C) 2002-2007 Alexei Drummond and Andrew Rambaut
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

package dr.inference.operators;

import dr.evolution.tree.NodeRef;
import dr.evomodel.tree.TreeModel;
import dr.inference.distribution.MultivariateDistributionLikelihood;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.math.WishartDistribution;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.SymmetricMatrix;
import dr.xml.*;

/**
 * @author Marc Suchard
 */
public class PrecisionMatrixGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {

    public static final String VARIANCE_OPERATOR = "precisionGibbsOperator";
    public static final String PRECISION_MATRIX = "precisionMatrix";
    public static final String TREE_MODEL = "treeModel";
    public static final String OUTCOME = "outcome";
    public static final String MEAN = "mean";
    public static final String PRIOR = "prior";

    //	private Parameter outcomeParam;
    //	private Parameter meanParam;
    private MatrixParameter precisionParam;
    private WishartDistribution priorDistribution;
    private int priorDf;
    //	private double[][] priorScaleMatrix;
    private SymmetricMatrix priorScaleMatrix;
    private TreeModel treeModel;
    private int dim;
    private int numberObservations;
//	private int weight;

    public PrecisionMatrixGibbsOperator(//Parameter outcomeParam, Parameter meanParam,
                                        MatrixParameter precisionParam,
                                        WishartDistribution priorDistribution,
                                        TreeModel treeModel,
                                        double weight) {
        super();
//		this.outcomeParam = outcomeParam;
//		this.meanParam = meanParam;
        this.precisionParam = precisionParam;
        this.priorDistribution = priorDistribution;
        this.priorDf = priorDistribution.df();
//		this.priorScaleMatrix = (SymmetricMatrix)(new SymmetricMatrix(priorDistribution.inverseScaleMatrix())).inverse();
        this.priorScaleMatrix = new SymmetricMatrix(priorDistribution.inverseScaleMatrix());
        setWeight(weight);
        this.treeModel = treeModel;
        dim = precisionParam.getRowDimension(); // assumed to be square
        numberObservations = treeModel.getNodeCount() - 1; // do not count the root


    }


    public int getStepCount() {
        return 1;
    }

    private void incrementsOuterProduct(double[][] S, NodeRef node) {

//		System.err.println("CALLED");
        if (!treeModel.isRoot(node)) {
            NodeRef parent = treeModel.getParent(node);
            double[] parentTrait = treeModel.getMultivariateNodeTrait(parent, "trait");
            // todo fix trait name
            double[] childTrait = treeModel.getMultivariateNodeTrait(node, "trait");
            double time = treeModel.getBranchLength(node);

            double sqrtTime = Math.sqrt(time);

            double[] delta = new double[dim];

            for (int i = 0; i < dim; i++)
                delta[i] = (childTrait[i] - parentTrait[i]) / sqrtTime; // todo check is * time?

            for (int i = 0; i < dim; i++) {            // symmetric matrix,
                for (int j = i; j < dim; j++)
                    S[j][i] = S[i][j] += delta[i] * delta[j];
//				for(int j=0; j<dim; j++)
//					S[i][j] += delta[i] * delta[j];
            }
        }
        // recurse down tree
        for (int i = 0; i < treeModel.getChildCount(node); i++)
            incrementsOuterProduct(S, treeModel.getChild(node, i));
    }


    public double doOperation() throws OperatorFailedException {
        // calculate sum-of-the-weighted-squares matrix over tree

//		System.err.println("DRAWING NEW PRECISION");
        double[][] S = new double[dim][dim];
        SymmetricMatrix S2 = null;
        SymmetricMatrix inverseS2 = null;
//		Matrix S = new Matrix(dim,dim)l
        incrementsOuterProduct(S, treeModel.getRoot());

//		System.err.println("S = \n"+(new SymmetricMatrix(S)).toString());
        try {
            S2 = priorScaleMatrix.add(new SymmetricMatrix(S));
            inverseS2 = (SymmetricMatrix) S2.inverse();
//			System.err.println(priorScaleMatrix.toString());
//			System.err.println(S2.toString());
//			System.err.println(inverseS2.toString());

        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }

        int df = priorDf + numberObservations;

        double[][] draw = WishartDistribution.nextWishart(df, inverseS2.toComponents());
//		System.err.println("SETTING VALUES");

        for (int i = 0; i < dim; i++) {
            Parameter column = precisionParam.getParameter(i);
            for (int j = 0; j < dim; j++)
                column.setParameterValueQuietly(j, draw[j][i]);
        }
//		System.err.println("FIRING CHANGE");
        precisionParam.fireParameterChangedEvent();
//		System.err.println("CHANGE DONE");

//		System.exit(0);


        return 0;  //To change body of implemented methods use File | Settings | File Templates.

    }

    public String getPerformanceSuggestion() {
        return null;  //To change body of implemented methods use File | Settings | File Templates.
    }

    public String getOperatorName() {
        return VARIANCE_OPERATOR;
    }

    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        public String getParserName() {
            return VARIANCE_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

//			int mode = CoercableMCMCOperator.DEFAULT;
//
//			if (xo.hasAttribute(AUTO_OPTIMIZE)) {
//				if (xo.getBooleanAttribute(AUTO_OPTIMIZE)) {
//					mode = CoercableMCMCOperator.COERCION_ON;
//				} else {
//					mode = CoercableMCMCOperator.COERCION_OFF;
//				}
//			}

            double weight = xo.getDoubleAttribute(WEIGHT);
//			double scaleFactor = xo.getDoubleAttribute(SCALE_FACTOR);
//
//			if (scaleFactor <= 0.0) {
//				throw new XMLParseException("scaleFactor must be greater than 0.0");
//			}

//			XMLObject cxo = (XMLObject) xo.getChild(OUTCOME);
//			Parameter outcomeParam = (Parameter) cxo.getChild(Parameter.class);

//			cxo = (XMLObject) xo.getChild(MEAN);
//			Parameter meanParam = (Parameter) cxo.getChild(Parameter.class);

//			XMLObject cxo = (XMLObject) xo.getChild(TREE_MODEL);
            TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

//			cxo = (XMLObject) xo.getChild(PRECISION_MATRIX);
            MatrixParameter precMatrix = (MatrixParameter) xo.getChild(MatrixParameter.class);
            MultivariateDistributionLikelihood prior = (MultivariateDistributionLikelihood) xo.getChild(MultivariateDistributionLikelihood.class);
            if (!(prior.getDistribution() instanceof WishartDistribution))
                throw new RuntimeException("Only a Wishart distribution is conjugate for Gibbs sampling");

//			WishartDistribution prior = (WishartDistribution) xo.getChild(WishartDistribution.class);

            // Make sure precMatrix is square and dim(precMatrix) = dim(parameter)

            if (precMatrix.getColumnDimension() != precMatrix.getRowDimension())
                throw new XMLParseException("The variance matrix is not square");

//			if (precMatrix.getColumnDimension() != outcomeParam.getDimension())
//				throw new XMLParseException("The parameter and variance matrix have differing dimensions");

//			if (outcomeParam.getDimension() != meanParam.getDimension())
//				throw new XMLParseException("The outcome and mean parameters have differing dimensions.");

//			return new MultivariateNormalOperator(parameter, scaleFactor, precMatrix, weight, mode);
            return new PrecisionMatrixGibbsOperator(
                    precMatrix, (WishartDistribution) prior.getDistribution(), treeModel, weight);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "This element returns a multivariate normal random walk operator on a given parameter.";
        }

        public Class getReturnType() {
            return MCMCOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
//				AttributeRule.newDoubleRule(SCALE_FACTOR),
                AttributeRule.newDoubleRule(WEIGHT),
//				AttributeRule.newBooleanRule(AUTO_OPTIMIZE, true),
//				new ElementRule(OUTCOME,
//						new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
//				new ElementRule(MEAN,
//								new XMLSyntaxRule[]{new ElementRule(Parameter.class)}),
//				new ElementRule(TREE_MODEL,
//								new XMLSyntaxRule[]{
                new ElementRule(TreeModel.class),
//										)}),
//				new ElementRule(PRIOR,
//								new XMLSyntaxRule[]{new ElementRule(WishartDistribution.class)}),
//
//				new ElementRule(Parameter.class),
                new ElementRule(MultivariateDistributionLikelihood.class),
//				new ElementRule(PRECISION_MATRIX,
//						new XMLSyntaxRule[]{
                new ElementRule(MatrixParameter.class)
        };

    };

}

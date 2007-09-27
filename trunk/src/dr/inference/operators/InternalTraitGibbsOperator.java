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
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.MatrixParameter;
import dr.math.MathUtils;
import dr.math.MultivariateNormalDistribution;
import dr.xml.*;

/**
 * @author Marc Suchard
 */
public class InternalTraitGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {

	public static final String GIBBS_OPERATOR = "internalTraitGibbsOperator";
//	public static final String PRECISION_MATRIX = "precisionMatrix";
//	public static final String TREE_MODEL = "treeModel";
//	public static final String OUTCOME  = "outcome";
//	public static final String MEAN = "mean";
//	public static final String PRIOR = "prior";

//	private Parameter outcomeParam;
//	private Parameter meanParam;
	//	private MatrixParameter precisionParam;
	//	private WishartDistribution priorDistribution;
	//	private int priorDf;
	//	private double[][] priorScaleMatrix;
	//	private SymmetricMatrix priorScaleMatrix;
	private TreeModel treeModel;
	private MatrixParameter precisionMatrixParameter;
	private int dim;
//	private int numberObservations;
//	private int weight;

	public InternalTraitGibbsOperator(TreeModel treeModel, MultivariateDiffusionModel diffusionModel) {
		super();
		this.treeModel = treeModel;
//		System.err.println("D: "+diffusionModel.toString());
//		System.exit(0);
		precisionMatrixParameter = diffusionModel.getPrecisionMatrixParameter();
		dim = treeModel.getMultivariateNodeTrait(treeModel.getRoot(), "trait").length;
		// todo need to fix trait name
	}


	public int getStepCount() {
		return 1;
	}


	public double doOperation() throws OperatorFailedException {

		double[][] precision = precisionMatrixParameter.getParameterAsMatrix();

		NodeRef node = treeModel.getRoot();
//		System.err.println("CNT: "+treeModel.getInternalNodeCount());
		while (node == treeModel.getRoot()) {
			node = treeModel.getInternalNode(MathUtils.nextInt(
					treeModel.getInternalNodeCount()));
		} // select any internal node but the root.

		NodeRef parent = treeModel.getParent(node);
//		int childCount = treeModel.getChildCount(node);

		double[] weightedAverage = new double[dim];
//		double weightTotal = 0;

//		double[] weights = new double[childCount+1];
		double weight = 1.0 / treeModel.getBranchLength(node);
		double[] trait = treeModel.getMultivariateNodeTrait(parent, "trait");

		for (int i = 0; i < dim; i++)
			weightedAverage[i] = trait[i] * weight;
//		weightedAverage[i] = i;

		double weightTotal = weight;
		for (int j = 0; j < treeModel.getChildCount(node); j++) {
			NodeRef child = treeModel.getChild(node, j);
			trait = treeModel.getMultivariateNodeTrait(child, "trait");
			weight = 1.0 / treeModel.getBranchLength(child);
			for (int i = 0; i < dim; i++)
				weightedAverage[i] += trait[i] * weight;
			weightTotal += weight;
		}

		for (int i = 0; i < dim; i++) {
			weightedAverage[i] /= weightTotal;
			for (int j = i; j < dim; j++)
				precision[j][i] = precision[i][j] *= weightTotal;
		}

		double[] draw = MultivariateNormalDistribution.nextMultivariateNormal(
				weightedAverage, precision);

		treeModel.setMultivariateTrait(node, "trait", draw);
//		for(int i=0; i<dim; i++)

		return 0;  //To change body of implemented methods use File | Settings | File Templates.

	}

	public String getPerformanceSuggestion() {
		return null;
	}

	public String getOperatorName() {
		return GIBBS_OPERATOR;
	}

	public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

		public String getParserName() {
			return GIBBS_OPERATOR;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			double weight = xo.getIntegerAttribute(WEIGHT);

			TreeModel treeModel = (TreeModel) xo.getChild(TreeModel.class);

			MultivariateDiffusionModel diffusionModel = (MultivariateDiffusionModel) xo.getChild(MultivariateDiffusionModel.class);

			InternalTraitGibbsOperator operator = new InternalTraitGibbsOperator(treeModel, diffusionModel);
			operator.setWeight(weight);
			return operator;
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "This element returns a multivariate Gibbs operator on an internal node trait.";
		}

		public Class getReturnType() {
			return MCMCOperator.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
				AttributeRule.newDoubleRule(WEIGHT),
				new ElementRule(TreeModel.class),
				new ElementRule(MultivariateDiffusionModel.class)
		};

	};

}

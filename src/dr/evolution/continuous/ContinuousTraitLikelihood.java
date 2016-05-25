/*
 * ContinuousTraitLikelihood.java
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

package dr.evolution.continuous;

import dr.evolution.io.NewickImporter;
import dr.evolution.tree.*;
import dr.matrix.Matrix;
import dr.matrix.MutableMatrix;
import dr.geo.math.SphericalPolarCoordinates;

import java.io.StringReader;


/**
 * Calculates the likelihood of a trait and tree.
 *
 * @version $Id: ContinuousTraitLikelihood.java,v 1.8 2006/06/18 16:20:58 alexei Exp $
 *
 * @author Alexei Drummond
 */
public class ContinuousTraitLikelihood {
		
	/**
	 * Calculates the likelihood of the traits on the given tree.
	 * @param tree this is the tree on which to calculate the continuous attributes likelihood
	 */	
	public double calculateLikelihood(MutableTree tree, String[] attributes, Contrastable[] mles, double kappa) {
		
		ContrastedTraitNode contrastNode = new ContrastedTraitNode(tree, tree.getRoot(), attributes);
		contrastNode.calculateContrasts(kappa);
		
		for (int i =0; i < mles.length; i++) {
			mles[i] = contrastNode.getTraitValue(i);
		}

		return calculateTraitsLikelihood(contrastNode);
	}
	
	/**
	 * Calculates the likelihood of the contrasted node
	 * @param contrastNode this is the node that the likelihood is calculated for
	 */	
	private double calculateTraitsLikelihood(ContrastedTraitNode contrastNode) {
	
		int count = contrastNode.getTraitCount();
		
		if (count == 1) return calculateSingleTraitLikelihood(contrastNode);
		
		return calculateMultipleTraitsLikelihood(contrastNode, count);
	}
		
	
	/**
	 * Calculates the likelihood of the contrasted node
	 * @param contrastNode this is the node that the likelihood is calculated for
	 */	
	private double calculateMultipleTraitsLikelihood(ContrastedTraitNode contrastNode, int traitCount) {
		
		SimpleTree contrastTree = new SimpleTree(contrastNode);
		
		double[][] w = new double[traitCount][traitCount];
		for (int j =0; j < traitCount; j++) {
			for (int k = j; k < traitCount; k++) {
				double wjk = 0.0;
				for (int i = 0; i < contrastTree.getInternalNodeCount(); i++) {
					ContrastedTraitNode ctNode = (ContrastedTraitNode)contrastTree.getInternalNode(i);
					
					wjk += (ctNode.contrast[j] * ctNode.contrast[k]) / ctNode.contrastVariance;
				}
				
				
				//System.out.println("w["+j+"]["+k+"]="+wjk);
				wjk /= (double)contrastTree.getInternalNodeCount();
				w[j][k] = wjk;
				w[k][j] = wjk;
			}
		}
		
		MutableMatrix answer = Matrix.Util.createMutableMatrix(new double[1][1]);
		MutableMatrix temp = Matrix.Util.createMutableMatrix(w);
		double detW = 0.0;
		try {
			detW = Matrix.Util.det(temp);
		} catch (Matrix.NotSquareException nse) {  nse.printStackTrace(System.out); }	
	
		//System.out.println("W matrix");
		//System.out.println(temp);
		
		//System.out.println("|W|=" + detW);
		
		MutableMatrix invW = Matrix.Util.createMutableMatrix(w);
		try {
			Matrix.Util.invert(invW);
		} catch (Matrix.NotSquareException nse) { nse.printStackTrace(System.out); }	
	
		//System.out.println("inverse of W matrix");
		//System.out.println(invW);
	
		double logL = 0.0;
		int n = contrastTree.getInternalNodeCount() + 1;
		
		for (int i =0; i < contrastTree.getInternalNodeCount(); i++) {
			ContrastedTraitNode ctNode = (ContrastedTraitNode)contrastTree.getInternalNode(i);
			double[] contrasts = ctNode.getTraitContrasts();
			
			Matrix uT = Matrix.Util.createRowVector(contrasts);
			Matrix u = Matrix.Util.createColumnVector(contrasts);
			
			try {
				Matrix.Util.product(invW, u, temp);
				Matrix.Util.product(uT, temp, answer);
			} catch (Matrix.WrongDimensionException wde) { wde.printStackTrace(System.out); }	
			
			logL += answer.getElement(0,0) / ctNode.getContrastVariance();
			logL += traitCount * Math.log(ctNode.getContrastVariance());
		}	
		
		// root variance
		logL += traitCount * Math.log(contrastNode.getNodeVariance());
		logL += n * Math.log(detW);
		logL += n * traitCount * Math.log(2*Math.PI);
		logL = -logL / 2.0;

        //System.out.println("root node variance = " + contrastNode.getNodeVariance());

        return logL;
	}
	
	/**
	 * Calculate the likelihood of a single continuous trait on the given tree.
	 */
	private double calculateSingleTraitLikelihood(ContrastedTraitNode contrastNode) {
		
		SimpleTree contrastTree = new SimpleTree(contrastNode);
		
		double s2 = 0.0;
		double sssContrast = 0.0;
		double slogCV = 0.0;
        for (int i = 0; i < contrastTree.getInternalNodeCount(); i++) {
			ContrastedTraitNode ctNode = (ContrastedTraitNode)contrastTree.getInternalNode(i);
			double contrast = ctNode.getTraitContrasts()[0];
			double cv =  ctNode.getContrastVariance();
			sssContrast += (contrast * contrast) / cv;
			slogCV += Math.log(cv);
			if (ctNode.isRoot()) {
                slogCV += Math.log(ctNode.getNodeVariance());
            }
        }

        double tl = 0.0;
        for (int i = 0; i < contrastTree.getNodeCount(); i++) {
            NodeRef node = contrastTree.getNode(i);
            if (!contrastTree.isRoot(node)) {
                tl += contrastTree.getBranchLength(node);
            }
        }

        s2 = sssContrast / contrastTree.getInternalNodeCount();
		
		int n = contrastTree.getInternalNodeCount() + 1;
		
		double logL = n * Math.log(2.0*Math.PI*s2);
		logL += slogCV;

        logL += sssContrast / s2;
		logL = -logL / 2.0;
		
		return logL;
	} 
	
	class ContrastedTraitNode extends SimpleNode {
	
		public ContrastedTraitNode(MutableTree tree, NodeRef node, String[] attributeNames) {
			
			init(tree, node, attributeNames.length);

            if (!tree.isExternal(node)) {
				
				if (tree.getChildCount(node) != 2) { throw new IllegalArgumentException("Tree must be strictly bifurcating!"); }
			
				addChild(new ContrastedTraitNode(tree, tree.getChild(node, 0), attributeNames));
				addChild(new ContrastedTraitNode(tree, tree.getChild(node, 1), attributeNames));
			} else {
				for (int i =0; i < attributeNames.length; i++) {
					Object obj = tree.getNodeTaxon(node).getAttribute(attributeNames[i]);
					
					if (obj == null) throw new IllegalArgumentException("attribute " + attributeNames[i] + " does not exist in " + tree.getTaxonId(node.getNumber()));
					
					if (obj instanceof Number) {
						traitValue[i] = new Continuous(((Number)obj).doubleValue());
					} else if (obj instanceof String) {
						traitValue[i] = new Continuous(Double.parseDouble((String)obj));
					} else if (obj instanceof Continuous) {
						traitValue[i] = (Continuous)obj;
					} else if (obj instanceof SphericalPolarCoordinates) {
                        traitValue[i] = (SphericalPolarCoordinates)obj;
                    }

                    tree.setNodeAttribute(node,attributeNames[i], traitValue[i]);
                }
			}

            this.traitNames = attributeNames;
        }
		
		private void init(MutableTree tree, NodeRef node, int traitCount) {
			setHeight(tree.getNodeHeight(node));
			setRate(tree.getNodeRate(node));
			setId(tree.getTaxonId(node.getNumber()));
			setNumber(node.getNumber());
			setTaxon(tree.getNodeTaxon(node));

			contrast = new double[traitCount];
			contrastVariance = 0.0;
			traitValue = new Contrastable[traitCount];
			nodeVariance = 0.0;
            this.tree = tree;
            this.node = node;
        }
		
		public double[] getTraitContrasts() {
			return contrast;
		}
		
		public double getContrastVariance() { return contrastVariance; }
		
		public double getNodeVariance() { return nodeVariance; }
		
		public Contrastable getTraitValue(int traitIndex) {
			
		    return traitValue[traitIndex];
		}
		
		
		public int getTraitCount() { return traitValue.length; }
		
		/**
		 * Recursively calculate the contrast information for the continuous trait nodes.
		 */
		private void calculateContrasts(double kappa) {
			if (!isExternal()) {
	
				ContrastedTraitNode left = (ContrastedTraitNode)getChild(0);
				ContrastedTraitNode right = (ContrastedTraitNode)getChild(1);
			
				left.calculateContrasts(kappa);
				right.calculateContrasts(kappa);
				
				double leftNodeBranchVariance = left.nodeVariance + Math.pow(getHeight() - left.getHeight(),kappa);
				double rightNodeBranchVariance = right.nodeVariance + Math.pow(getHeight() - right.getHeight(), kappa);
				
				// calculate the contrast variances
				contrastVariance = leftNodeBranchVariance + rightNodeBranchVariance;
				
				// estimate the variance of the tree value at this node
				nodeVariance = (leftNodeBranchVariance * rightNodeBranchVariance) / (leftNodeBranchVariance + rightNodeBranchVariance);
				
				// calculate the weights for each child
				double invVarLeft = 1.0 / leftNodeBranchVariance;
				double invVarRight = 1.0 / rightNodeBranchVariance;
				
				// estimate the contrasts for each trait
				for (int i = 0; i < getTraitCount(); i++) {
					// calculate the contrast
					contrast[i] = left.traitValue[i].getDifference(right.traitValue[i]);
					//left.traitValue[i] - right.traitValue[i];
					
					// estimate the variance weighted mean of the two child observations of the ith trait
					traitValue[i] = left.traitValue[i].getWeightedMean(invVarLeft, left.traitValue[i], invVarRight, right.traitValue[i]);
					//traitValue[i] = (invVarLeft * left.traitValue[i] + invVarRight * right.traitValue[i]) / (invVarLeft + invVarRight);	

                    tree.setNodeAttribute(node, traitNames[i], traitValue[i]);
                }
			}
		}
	
		// the contrast for each trait at this node.
		private double[] contrast;
		
		// the contrast variance at this node.
		private double contrastVariance;
		
		// the trait value at this node.
		private Contrastable[] traitValue;
		
		// the tree variance at this node.
		private double nodeVariance;

        // the original tree
        private MutableTree tree;
        private NodeRef node;

        private String[] traitNames;
		
        // @todo find out what should be done with this variable
		// the sum of the log of the contrast variances below this node
        // private double slContrastVariance = 0;	
	}
	
	public static void main(String[] args) throws Exception {
		
		String testTree = "((A:1, B:1):1,(C:1, D:1):1);";
	
		NewickImporter newickImporter = new NewickImporter(new StringReader(testTree));
		
		MutableTree tree = (MutableTree)newickImporter.importTree(null);
	
		tree.setTaxonAttribute(0, "U1", new Continuous(1.10));
		tree.setTaxonAttribute(1, "U1", new Continuous(1.95));
		tree.setTaxonAttribute(2, "U1", new Continuous(3.15));
		tree.setTaxonAttribute(3, "U1", new Continuous(4.39));
		
		tree.setTaxonAttribute(0, "U2", new Continuous(5.2));
		tree.setTaxonAttribute(1, "U2", new Continuous(3.8));
		tree.setTaxonAttribute(2, "U2", new Continuous(3.1));
		tree.setTaxonAttribute(3, "U2", new Continuous(1.95));
		
		ContinuousTraitLikelihood ctLikelihood = new ContinuousTraitLikelihood();

		Contrastable[] mles = new Contrastable[2];
		double logL = ctLikelihood.calculateLikelihood(tree, new String[] {"U1", "U2"}, mles, 1.0);
			
				
		System.out.println("logL = " + logL);
		System.out.println("mle(trait1) = " + mles[0]);
		System.out.println("mle(trait2) = " + mles[1]);

		Contrastable[] mle = new Contrastable[1];
		System.out.println("logL (trait1) = " + ctLikelihood.calculateLikelihood(tree, new String[] {"U1"}, mle, 1.0));
		System.out.println("mle(trait1) = " + mle[0]);
		System.out.println("logL (trait2) = " + ctLikelihood.calculateLikelihood(tree, new String[] {"U2"}, mle, 1.0));
		System.out.println("mle(trait2) = " + mle[0]);


	}
}
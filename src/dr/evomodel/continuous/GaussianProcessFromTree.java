/*
 * GaussianProcessFromTree.java
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

import dr.evolution.tree.NodeRef;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.math.distributions.GaussianProcessRandomGenerator;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.CholeskyDecomposition;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.SymmetricMatrix;

/**
 * @author Marc A. Suchard
 * @author Max R. Tolkoff
 */
public class GaussianProcessFromTree implements GaussianProcessRandomGenerator {

    private final FullyConjugateMultivariateTraitLikelihood traitModel;

    public GaussianProcessFromTree(FullyConjugateMultivariateTraitLikelihood traitModel) {
        this.traitModel = traitModel;
    }

    public Likelihood getLikelihood() {
        return traitModel;
    }

    //    boolean firstTime=true;
    public double[] nextRandomFast() {

        double[] random = new double[traitModel.getTreeModel().getExternalNodeCount()*traitModel.getDimTrait()];
        NodeRef root = traitModel.getTreeModel().getRoot();
        double[] traitStart=traitModel.getPriorMean();
        double[][] varianceCholesky=null;
        double[][] temp= new SymmetricMatrix(traitModel.getDiffusionModel().getPrecisionmatrix()).inverse().toComponents();
        try {
            varianceCholesky = (new CholeskyDecomposition(temp).getL());
        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }
//        if(traitModel.getTreeModel().isExternal(root)) {
//            random[0] = traitModel.getTreeModel().getMultivariateNodeTrait(root, traitModel.getTraitName())[i];
//        }
//        else{


//        double[][] var = MultivariateTraitUtils.computeTreeVariance(traitModel, true);
//        if(firstTime) {
//            for (int j = 0; j < var[0].length; j++) {
//                for (int k = 0; k < var[0].length; k++) {
//                    if(j!=k)
//                        var[j][k] = var[j][k] / Math.sqrt(var[k][k] * var[j][j]);
//                }
//            }
//
//
//
//            for (int j = 0; j < var[0].length; j++) {
//                String empty = "";
//                for (int k = 0; k < var[0].length; k++) {
//                    empty += Double.toString(var[j][k]) + "\t";
//                }
//                System.out.println(empty);
//            }
//            firstTime=false;
//        }
        nextRandomFast(traitStart, root, random, varianceCholesky);
//        }
        return random;
    }

    private void nextRandomFast(double[] currentValue, NodeRef currentNode, double[] random, double[][] varianceCholesky) {

        double rescaledLength = (traitModel.getTreeModel().isRoot(currentNode)) ?
            1.0 / traitModel.getPriorSampleSize() :
                traitModel.getRescaledBranchLengthForPrecision(currentNode);

        double scale = Math.sqrt(rescaledLength);

        // draw ~ MNV(mean = currentVale, variance = scale * scale * L^t L)
        double[] draw = MultivariateNormalDistribution.nextMultivariateNormalCholesky(currentValue, varianceCholesky, scale);

        if (traitModel.getTreeModel().isExternal(currentNode)) {
            System.arraycopy(draw, 0, random, currentNode.getNumber() * draw.length, draw.length);
        } else {
            int childCount = traitModel.getTreeModel().getChildCount(currentNode);
            for (int i = 0; i < childCount; i++) {
                nextRandomFast(draw, traitModel.getTreeModel().getChild(currentNode, i), random, varianceCholesky);
            }
        }
    }

    @Override
    public Object nextRandom() {
        return nextRandomFast();
    }

    @Override
    public double logPdf(Object x) {

        double[] v = (double[]) x;
        Parameter variable = traitModel.getTraitParameter();
        for (int i = 0; i < v.length; ++i) {
            variable.setParameterValueQuietly(i, v[i]);
        }
        variable.fireParameterChangedEvent();

        return traitModel.getLogLikelihood();
    }
}

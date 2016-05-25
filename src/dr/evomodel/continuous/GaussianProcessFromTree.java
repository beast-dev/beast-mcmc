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
import dr.math.KroneckerOperation;
import dr.math.distributions.GaussianProcessRandomGenerator;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.CholeskyDecomposition;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
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

    @Override
    public Likelihood getLikelihood() {
        return traitModel;
    }

    @Override
    public int getDimension() {
        return traitModel.getTreeModel().getExternalNodeCount() * traitModel.getDimTrait();
    }

    @Override
    public double[][] getPrecisionMatrix() {
        final boolean includeRoot = false; // TODO make an option

        double[][] treeVariance;
//        long startTime1 = System.nanoTime();
        treeVariance = traitModel.computeTreeVariance2(includeRoot);
//        long estimatedTime1 = System.nanoTime() - startTime1;

//        long startTime2 = System.nanoTime();
//        treeVariance = traitModel.computeTreeVariance(includeRoot);
//        long estimatedTime2 = System.nanoTime() - startTime2;

        double[][] traitPrecision = traitModel.getDiffusionModel().getPrecisionmatrix();


//        for (int i = 0; i < treeVariance2.length; ++i) {
//            for (int j = 0; j < treeVariance2[i].length; ++j) {
//                if (treeVariance2[i][j] != treeVariance[i][j]) {
//                    System.err.println(i + " " + j);
//                    System.err.println(treeVariance2[i][j] + " " + treeVariance[i][j]);
//                    System.exit(-1);
//                }
//            }
//        }

//        System.err.println("T1: " + estimatedTime1);
//        System.err.println("T2: " + estimatedTime2);

//        System.err.println("\t\tSTART prec");
        Matrix treePrecision = new Matrix(treeVariance).inverse();

//        System.err.println("\t\tSTART kron");

        double[][] jointPrecision = KroneckerOperation.product(treePrecision.toComponents(), traitPrecision); // TODO Double-check order

        return jointPrecision;
    }

    private static void scale(double[][] matrix, double scale) {
        for (int i = 0; i < matrix.length; ++i) {
            for (int j = 0; j < matrix[i].length; ++j) {
                matrix[i][j] *= scale;
            }
        }
    }

    public double getLogLikelihood() { return traitModel.getLogLikelihood(); }

    //    boolean firstTime=true;
    public double[] nextRandomFast() {

        double[] random = new double[traitModel.getTreeModel().getExternalNodeCount() * traitModel.getDimTrait()];
        NodeRef root = traitModel.getTreeModel().getRoot();
        double[] traitStart = traitModel.getPriorMean();
        double[][] varianceCholesky = null;
        double[][] temp = new SymmetricMatrix(traitModel.getDiffusionModel().getPrecisionmatrix()).inverse().toComponents();
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
        if (USE_BUFFER) {
            final int length = traitModel.getDimTrait();
            final int nodeCount = traitModel.getTreeModel().getNodeCount();
            double[] currentValue = new double[(nodeCount + 1) * length];
            double[] epsilon = new double[length];
            final int priorOffset = nodeCount * length;
            System.arraycopy(traitStart, 0, currentValue, priorOffset, length);
            nextRandomFast2(currentValue, priorOffset, root, random, varianceCholesky, epsilon);
        } else {
            nextRandomFast(traitStart, root, random, varianceCholesky);
        }
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

    private void nextRandomFast2(double[] currentValue, int parentOffset, NodeRef currentNode, double[] random,
                                 double[][] varianceCholesky, double[] epsilon) {

        final int length = varianceCholesky.length;

        double rescaledLength = (traitModel.getTreeModel().isRoot(currentNode)) ?
                1.0 / traitModel.getPriorSampleSize() :
                traitModel.getRescaledBranchLengthForPrecision(currentNode);

        double scale = Math.sqrt(rescaledLength);

        final int currentOffset = currentNode.getNumber() * length;

        // draw ~ MNV(mean = currentValue at parent, variance = scale * scale * L^t L)
        MultivariateNormalDistribution.nextMultivariateNormalCholesky(
                currentValue, parentOffset, // mean at parent
                varianceCholesky, scale,
                currentValue, currentOffset, // result at current
                epsilon);

        if (traitModel.getTreeModel().isExternal(currentNode)) {
            System.arraycopy(
                    currentValue, currentOffset, // result at tip
                    random, currentOffset, // into final results buffer
                    length);
        } else {
            int childCount = traitModel.getTreeModel().getChildCount(currentNode);
            for (int i = 0; i < childCount; i++) {
                nextRandomFast2(
                        currentValue, currentOffset,
                        traitModel.getTreeModel().getChild(currentNode, i),
                        random, varianceCholesky, epsilon);
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

    private static final boolean USE_BUFFER = true;
}

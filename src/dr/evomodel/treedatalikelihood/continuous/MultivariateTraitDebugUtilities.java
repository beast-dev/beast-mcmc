/*
 * MultivariateTraitDebugUtilities.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.BranchRates;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.evomodel.treedatalikelihood.continuous.cdi.SafeMultivariateActualizedWithDriftIntegrator;
import dr.math.KroneckerOperation;
import dr.math.matrixAlgebra.Matrix;
import org.ejml.data.Complex64F;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.interfaces.decomposition.EigenDecomposition;
import org.ejml.ops.CommonOps;
import org.ejml.ops.EigenOps;

import java.util.HashSet;
import java.util.Set;

/**
 * @author Marc A. Suchard
 */
public class MultivariateTraitDebugUtilities {

    public static double getLengthToRoot(final Tree tree, final BranchRates branchRates, final NodeRef nodeRef) {

        double length = 0;

        if (!tree.isRoot(nodeRef)) {
            NodeRef parent = tree.getParent(nodeRef);
            double rateScalar = 1.0;
            if (branchRates != null) {
                rateScalar = branchRates.getBranchRate(tree, nodeRef);
            }
            length += rateScalar * tree.getBranchLength(nodeRef) + getLengthToRoot(tree, branchRates, parent);
        }

        return length;
    }

    private static NodeRef findMRCA(final Tree tree, final int iTip, final int jTip) {
        Set<String> leafNames = new HashSet<String>();
        leafNames.add(tree.getTaxonId(iTip));
        leafNames.add(tree.getTaxonId(jTip));
        return TreeUtils.getCommonAncestorNode(tree, leafNames);
    }

    public static void insertPrecision(Tree tree, NodeRef parent, NodeRef child, double[][] matrix, double normalization) {
        matrix[parent.getNumber()][child.getNumber()] =
                           matrix[child.getNumber()][parent.getNumber()] =-1.0 / (tree.getBranchLength(child) * normalization);
        recurseGraph(tree, child, matrix, normalization);

    }

    public static void recurseGraph(Tree tree, NodeRef node, double[][] matrix, double normalization) {
        if (!tree.isExternal(node)) {
            insertPrecision(tree, node, tree.getChild(node, 0), matrix, normalization);
            insertPrecision(tree, node, tree.getChild(node, 1), matrix, normalization);
        }
    }

    public static double[][] getGraphVariance(final Tree tree, final BranchRates branchRates, final double normalization, final double priorSampleSize) {

        final int nodeCount = tree.getNodeCount();
//        final double[][] variance = new double[nodeCount + 1][nodeCount + 1];
        final double[][] variance = new double[nodeCount][nodeCount];

        for (int i = 0; i < nodeCount; i++) {

            // Fill in diagonal
            double marginalTime = getLengthToRoot(tree, branchRates, tree.getNode(i)) * normalization;
            variance[i][i] = marginalTime;

            // Fill in upper right triangle,
            for (int j = i + 1; j < nodeCount; j++) {
                NodeRef mrca = TreeUtils.getCommonAncestorSafely(tree, tree.getNode(i), tree.getNode(j));
                variance[i][j] = getLengthToRoot(tree, branchRates, mrca) * normalization;
            }
        }

        makeSymmetric(variance);

        addPrior(variance, priorSampleSize);

        return variance;
    }

    public static double[][] getTreeVariance(final Tree tree, final BranchRates branchRates, final double normalization, final double priorSampleSize) {

        final int tipCount = tree.getExternalNodeCount();

        double[][] variance = new double[tipCount][tipCount];

        for (int i = 0; i < tipCount; i++) {

            // Fill in diagonal
            double marginalTime = getLengthToRoot(tree, branchRates, tree.getExternalNode(i)) * normalization;
            variance[i][i] = marginalTime;

            // Fill in upper right triangle,

            for (int j = i + 1; j < tipCount; j++) {
                NodeRef mrca = findMRCA(tree, i, j);
                variance[i][j] = getLengthToRoot(tree, branchRates, mrca) * normalization;
            }
        }

        makeSymmetric(variance);

        addPrior(variance, priorSampleSize);

        return variance;
    }

    private static void makeSymmetric(final double[][] variance) {
        for (int i = 0; i < variance.length; i++) {
            for (int j = i + 1; j < variance[i].length; j++) {
                variance[j][i] = variance[i][j];
            }
        }
    }

    public static double[][] getJointVariance(final Tree tree, final double normalization,
                                              final double priorSampleSize, final TreeDataLikelihood callbackLikelihood,
                                              Matrix traitVariance, StringBuilder sb, DiffusionProcessDelegate diffusionProcessDelegate) {

        if (diffusionProcessDelegate instanceof DiagonalOrnsteinUhlenbeckDiffusionModelDelegate) {
            // Eigen of strength of selection matrix
            double[] eigVals = ((DiagonalOrnsteinUhlenbeckDiffusionModelDelegate) diffusionProcessDelegate).getStrengthOfSelection();
            int n = eigVals.length;

            // Computation of matrix
            DenseMatrix64F transTraitVariance = new DenseMatrix64F(traitVariance.toComponents());

            double[][] treeSharedLengths = getTreeVariance(tree, callbackLikelihood.getBranchRateModel(), normalization, Double.POSITIVE_INFINITY);
            int ntaxa = tree.getExternalNodeCount();
            double ti; double tj; double tij; double ep; double eq;
            DenseMatrix64F varTemp = new DenseMatrix64F(n,n);
            double[][] jointVariance = new double[n * ntaxa][n * ntaxa];
            for (int i = 0; i < ntaxa; ++i) {
                for (int j = 0; j < ntaxa; ++j) {
                    ti = treeSharedLengths[i][i];
                    tj = treeSharedLengths[j][j];
                    tij = treeSharedLengths[i][j];
                    for (int p = 0; p < n; ++p) {
                        for (int q = 0; q < n; ++q) {
                            ep = eigVals[p];
                            eq = eigVals[q];
                            varTemp.set(p, q, Math.exp(-ep * ti) * Math.exp(-eq * tj) * ((Math.exp((ep + eq) * tij) - 1)/(ep + eq) + 1/priorSampleSize) * transTraitVariance.get(p, q));
                        }
                    }
                    for (int p = 0; p < n; ++p) {
                        for (int q = 0; q < n; ++q) {
                            jointVariance[i * n + p][j * n + q] = varTemp.get(p, q);
                        }
                    }
                }
            }

            sb.append("Trait variance:\n");
            sb.append(traitVariance);
            sb.append("\n\n");

            sb.append("Joint variance:\n");
            sb.append(new Matrix(jointVariance));
            sb.append("\n\n");

            return jointVariance;

        } else {
            if (diffusionProcessDelegate instanceof OrnsteinUhlenbeckDiffusionModelDelegate) {
                // Eigen of strength of selection matrix
                double[][] alphaMat = ((OrnsteinUhlenbeckDiffusionModelDelegate) diffusionProcessDelegate).getStrengthOfSelection();
                DenseMatrix64F A = new DenseMatrix64F(alphaMat);
                int n = A.numCols;
                if (n != A.numRows) throw new RuntimeException("Selection strength A matrix must be square.");
                EigenDecomposition eigA = DecompositionFactory.eig(n, true);
                if (!eigA.decompose(A)) throw new RuntimeException("Eigen decomposition failed.");

                // Transformed variance
                DenseMatrix64F V = EigenOps.createMatrixV(eigA);
                DenseMatrix64F Vinv = new DenseMatrix64F(n, n);
                CommonOps.invert(V, Vinv);

                DenseMatrix64F transTraitVariance = new DenseMatrix64F(traitVariance.toComponents());

                DenseMatrix64F tmp = new DenseMatrix64F(n, n);
                CommonOps.mult(Vinv, transTraitVariance, tmp);
                CommonOps.multTransB(tmp, Vinv, transTraitVariance);

                // Eigenvalues
                double[] eigVals = new double[n];
                for (int p = 0; p < n; ++p) {
                    Complex64F ev = eigA.getEigenvalue(p);
                    if (!ev.isReal())
                        throw new RuntimeException("Selection strength A should only have real eigenvalues.");
                    eigVals[p] = ev.real;
                }

                // inverse of eigenvalues
                double[][] invEigVals = new double[n][n];
                for (int p = 0; p < n; ++p) {
                    for (int q = 0; q < n; ++q) {
                        invEigVals[p][q] = 1 / (eigVals[p] + eigVals[q]);
                    }
                }

                // Computation of matrix
                double[][] treeSharedLengths = getTreeVariance(tree, callbackLikelihood.getBranchRateModel(), normalization, Double.POSITIVE_INFINITY);
                int ntaxa = tree.getExternalNodeCount();
                double ti;
                double tj;
                double tij;
                double ep;
                double eq;
                DenseMatrix64F varTemp = new DenseMatrix64F(n, n);
                double[][] jointVariance = new double[n * ntaxa][n * ntaxa];
                for (int i = 0; i < ntaxa; ++i) {
                    for (int j = 0; j < ntaxa; ++j) {
                        ti = treeSharedLengths[i][i];
                        tj = treeSharedLengths[j][j];
                        tij = treeSharedLengths[i][j];
                        for (int p = 0; p < n; ++p) {
                            for (int q = 0; q < n; ++q) {
                                ep = eigVals[p];
                                eq = eigVals[q];
                                varTemp.set(p, q, Math.exp(-ep * ti) * Math.exp(-eq * tj) * (invEigVals[p][q] * (Math.exp((ep + eq) * tij) - 1) + 1 / priorSampleSize) * transTraitVariance.get(p, q));
                            }
                        }
                        CommonOps.mult(V, varTemp, tmp);
                        CommonOps.multTransB(tmp, V, varTemp);
                        for (int p = 0; p < n; ++p) {
                            for (int q = 0; q < n; ++q) {
                                jointVariance[i * n + p][j * n + q] = varTemp.get(p, q);
                            }
                        }
                    }
                }

                sb.append("Trait variance:\n");
                sb.append(traitVariance);
                sb.append("\n\n");

                sb.append("Joint variance:\n");
                sb.append(new Matrix(jointVariance));
                sb.append("\n\n");

                return jointVariance;

            } else {
                double[][] treeVariance = MultivariateTraitDebugUtilities.getTreeVariance(tree, callbackLikelihood.getBranchRateModel(), normalization, priorSampleSize);

                double[][] jointVariance = KroneckerOperation.product(treeVariance, traitVariance.toComponents());

                Matrix treeV = new Matrix(treeVariance);
                Matrix treeP = treeV.inverse();

                sb.append("Tree variance:\n");
                sb.append(treeV);
                sb.append("Tree precision:\n");
                sb.append(treeP);
                sb.append("\n\n");

                sb.append("Trait variance:\n");
                sb.append(traitVariance);
                sb.append("\n\n");

                sb.append("Joint variance:\n");
                sb.append(new Matrix(jointVariance));
                sb.append("\n\n");

                return jointVariance;
            }
        }
    }

    private static void addPrior(final double[][] variance, double priorSampleSize) {
        if (!Double.isInfinite(priorSampleSize)) {
            for (int i = 0; i < variance.length; ++i) {
                for (int j = 0; j < variance[i].length; ++j) {
                    variance[i][j] += 1.0 / priorSampleSize;
                }
            }
        }
    }

    public static double[][] getTreeDrift(Tree tree, DiffusionProcessDelegate diffusion, double[] priorMean, StringBuilder sb) {

        final int dim = diffusion.getDiffusionModel(0).getPrecisionParameter().getColumnDimension();
        final double[][] drift = new double[tree.getExternalNodeCount()][dim];

        if (diffusion instanceof DriftDiffusionModelDelegate) {
            for (int tip = 0; tip < tree.getExternalNodeCount(); ++tip) {
                drift[tip] = ((DriftDiffusionModelDelegate) diffusion).getAccumulativeDrift(
                        tree.getExternalNode(tip), priorMean);
            }
        }

        if (diffusion instanceof OrnsteinUhlenbeckDiffusionModelDelegate) {
            for (int tip = 0; tip < tree.getExternalNodeCount(); ++tip) {
                drift[tip] = ((OrnsteinUhlenbeckDiffusionModelDelegate) diffusion).getAccumulativeDrift(
                        tree.getExternalNode(tip), priorMean);
            }
        }

        if (diffusion instanceof DiagonalOrnsteinUhlenbeckDiffusionModelDelegate) {
            for (int tip = 0; tip < tree.getExternalNodeCount(); ++tip) {
                drift[tip] = ((DiagonalOrnsteinUhlenbeckDiffusionModelDelegate) diffusion).getAccumulativeDrift(
                        tree.getExternalNode(tip), priorMean);
            }
        }

        if (diffusion.hasDrift()) {
            sb.append("Tree drift (including root mean):\n");
            sb.append(new Matrix(drift));
            sb.append("\n\n");
        }

        return drift;
    }

    public static double[][] getGraphDrift(Tree tree, DiffusionProcessDelegate diffusion) {

        final int dim = diffusion.getDiffusionModel(0).getPrecisionParameter().getColumnDimension();
        final double[][] drift = new double[tree.getNodeCount()][dim];

        final double[] temp = new double[dim];

        if (diffusion instanceof DriftDiffusionModelDelegate) {
            for (int node = 0; node < tree.getNodeCount(); ++node) {
                drift[node] = ((DriftDiffusionModelDelegate) diffusion).getAccumulativeDrift(
                        tree.getNode(node), temp);
            }
        }

        return drift;
    }
}

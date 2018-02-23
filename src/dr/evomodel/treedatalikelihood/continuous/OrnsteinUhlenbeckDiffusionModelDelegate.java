/*
 * OrnsteinUhlenbeckDiffusionModelDelegate.java
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

package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.inference.model.MatrixParameterInterface;
import dr.math.matrixAlgebra.missingData.MissingOps;
import org.ejml.data.Complex64F;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.interfaces.decomposition.EigenDecomposition;
import org.ejml.ops.CommonOps;
import org.ejml.ops.EigenOps;
import org.ejml.ops.MatrixFeatures;

import java.util.List;

import static dr.math.matrixAlgebra.missingData.MissingOps.wrap;

/**
 * A simple OU diffusion model delegate with branch-specific drift and constant diffusion
 *
 * @author Marc A. Suchard
 * @author Paul Bastide
 * @version $Id$
 */
public final class OrnsteinUhlenbeckDiffusionModelDelegate extends AbstractOUDiffusionModelDelegate {

    // NOTE TO PB: Need to merge DriftDiffusion*, OrnsteinUhlenbeck* into single class and then delegate specialized work


    private MatrixParameterInterface strengthOfSelectionMatrixParameter;
    private EigenDecomposition eigenDecompositionStrengthOfSelection;


    public OrnsteinUhlenbeckDiffusionModelDelegate(Tree tree,
                                                   MultivariateDiffusionModel diffusionModel,
                                                   List<BranchRateModel> branchRateModels,
                                                   MatrixParameterInterface strengthOfSelectionMatrixParam) {
        this(tree, diffusionModel, branchRateModels, strengthOfSelectionMatrixParam, 0);
    }

    private OrnsteinUhlenbeckDiffusionModelDelegate(Tree tree,
                                                    MultivariateDiffusionModel diffusionModel,
                                                    List<BranchRateModel> branchRateModels,
                                                    MatrixParameterInterface strengthOfSelectionMatrixParam,
                                                    int partitionNumber) {
        super(tree, diffusionModel, branchRateModels, partitionNumber);

        // Strength of selection matrix
        this.strengthOfSelectionMatrixParameter = strengthOfSelectionMatrixParam;
        addVariable(strengthOfSelectionMatrixParameter);
        // Eigen decomposition
        this.eigenDecompositionStrengthOfSelection = decomposeStrenghtOfSelection(strengthOfSelectionMatrixParam);
    }

    private EigenDecomposition decomposeStrenghtOfSelection(MatrixParameterInterface Aparam) {
        DenseMatrix64F A = MissingOps.wrap(Aparam);
        int n = A.numCols;
        // Checks
        if (n != A.numRows) throw new RuntimeException("Selection strength A matrix must be square.");
        if (!MatrixFeatures.isSymmetric(A))
            throw new RuntimeException("Selection strength A matrix must be symmetric."); // TODO : this is not strictly necessary, but might be good to impose the constraint ?
        // Decomposition
        EigenDecomposition eigA = DecompositionFactory.eig(n, true, true);
        if (!eigA.decompose(A)) throw new RuntimeException("Eigen decomposition failed.");
        return eigA;
    }

    @Override
    public double[][] getStrengthOfSelection() {
        return strengthOfSelectionMatrixParameter.getParameterAsMatrix();
    }

    @Override
    public double[] getEigenValuesStrengthOfSelection() {
        double[] eigA = new double[dim];
        for (int p = 0; p < dim; ++p) {
            Complex64F ev = eigenDecompositionStrengthOfSelection.getEigenvalue(p);
            if (!ev.isReal()) throw new RuntimeException("Selection strength A should only have real eigenvalues.");
            eigA[p] = ev.real;
        }
        return eigA;
    }

    @Override
    public double[] getEigenVectorsStrengthOfSelection() {
        return EigenOps.createMatrixV(eigenDecompositionStrengthOfSelection).data;
    }


    @Override
    public double[] getAccumulativeDrift(final NodeRef node, double[] priorMean, ContinuousDiffusionIntegrator cdi) {
        final DenseMatrix64F drift = new DenseMatrix64F(dim, 1, true, priorMean);
        double[] displacement = new double[dim];
        double[] actualization = new double[dim * dim];
        recursivelyAccumulateDrift(node, drift, cdi, displacement, actualization);
        return drift.data;
    }

    private void recursivelyAccumulateDrift(final NodeRef node, final DenseMatrix64F drift, ContinuousDiffusionIntegrator cdi, double[] displacement, double[] actualization) {
        if (!tree.isRoot(node)) {

            // Compute parent
            recursivelyAccumulateDrift(tree.getParent(node), drift, cdi, displacement, actualization);

            // Node
            cdi.getBranchDisplacement(getMatrixBufferOffsetIndex(node.getNumber()), displacement);
            cdi.getBranchActualization(getMatrixBufferOffsetIndex(node.getNumber()), actualization);

            DenseMatrix64F temp = new DenseMatrix64F(dim, 1);
            CommonOps.mult(new DenseMatrix64F(dim, dim, true, actualization), drift, temp);
            CommonOps.add(temp, new DenseMatrix64F(dim, 1, true, displacement), drift);

        }
    }

    @Override
    public double[][] getJointVariance(final double priorSampleSize, final double[][] treeVariance, final double[][] treeSharedLengths, final double[][] traitVariance) {

        double[] eigVals = this.getEigenValuesStrengthOfSelection();
        DenseMatrix64F V = wrap(this.getEigenVectorsStrengthOfSelection(), 0, dim, dim);
        DenseMatrix64F Vinv = new DenseMatrix64F(dim, dim);
        CommonOps.invert(V, Vinv);

        DenseMatrix64F transTraitVariance = new DenseMatrix64F(traitVariance);

        DenseMatrix64F tmp = new DenseMatrix64F(dim, dim);
        CommonOps.mult(Vinv, transTraitVariance, tmp);
        CommonOps.multTransB(tmp, Vinv, transTraitVariance);

        // inverse of eigenvalues
        double[][] invEigVals = new double[dim][dim];
        for (int p = 0; p < dim; ++p) {
            for (int q = 0; q < dim; ++q) {
                invEigVals[p][q] = 1 / (eigVals[p] + eigVals[q]);
            }
        }

        // Computation of matrix
        int ntaxa = tree.getExternalNodeCount();
        double ti;
        double tj;
        double tij;
        double ep;
        double eq;
        DenseMatrix64F varTemp = new DenseMatrix64F(dim, dim);
        double[][] jointVariance = new double[dim * ntaxa][dim * ntaxa];
        for (int i = 0; i < ntaxa; ++i) {
            for (int j = 0; j < ntaxa; ++j) {
                ti = treeSharedLengths[i][i];
                tj = treeSharedLengths[j][j];
                tij = treeSharedLengths[i][j];
                for (int p = 0; p < dim; ++p) {
                    for (int q = 0; q < dim; ++q) {
                        ep = eigVals[p];
                        eq = eigVals[q];
                        varTemp.set(p, q, Math.exp(-ep * ti) * Math.exp(-eq * tj) * (invEigVals[p][q] * (Math.exp((ep + eq) * tij) - 1) + 1 / priorSampleSize) * transTraitVariance.get(p, q));
                    }
                }
                CommonOps.mult(V, varTemp, tmp);
                CommonOps.multTransB(tmp, V, varTemp);
                for (int p = 0; p < dim; ++p) {
                    for (int q = 0; q < dim; ++q) {
                        jointVariance[i * dim + p][j * dim + q] = varTemp.get(p, q);
                    }
                }
            }
        }
        return jointVariance;
    }
}
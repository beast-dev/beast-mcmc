/*
 * HomogenousDiffusionModelDelegate.java
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
import dr.evomodel.treedatalikelihood.continuous.cdi.SafeMultivariateWithDriftIntegrator;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Model;
import dr.math.matrixAlgebra.missingData.MissingOps;
import org.ejml.data.Complex64F;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.interfaces.decomposition.EigenDecomposition;
import org.ejml.ops.CommonOps;
import org.ejml.ops.EigenOps;
import org.ejml.ops.MatrixFeatures;

import java.util.List;

/**
 * A simple OU diffusion model delegate with branch-specific drift and constant diffusion
 * @author Marc A. Suchard
 * @author Paul Bastide
 * @version $Id$
 */
public final class OrnsteinUhlenbeckDiffusionModelDelegate extends AbstractDiffusionModelDelegate {

    // QUESTION: would it be better to nest this under DriftDiffusionModelDelegate ?
    // it was declared final, so I didn't dare to do it...

    // Here, branchRateModels represents optimal values


    private final int dim;
    private final List<BranchRateModel> branchRateModels;

    protected MatrixParameterInterface strengthOfSelectionMatrixParameter;
    protected EigenDecomposition eigenDecompositionStrengthOfSelection;

//    private double[][] strengthOfSelectionMatrix;
//    private double[][] savedStrengthOfSelectionMatrix;

//    private double[][] stationaryPrecisionMatrix;
//    private double[][] savedStationaryPrecisionMatrix;

//    private final BufferIndexHelper matrixActualizationBufferHelper;

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
        super(tree, diffusionModel, partitionNumber);
        this.branchRateModels = branchRateModels;

        dim = diffusionModel.getPrecisionParameter().getColumnDimension();

        if (branchRateModels != null) {

            for (BranchRateModel rateModel : branchRateModels) {
                addModel(rateModel);
            }

            if (branchRateModels.size() != dim) {
                throw new IllegalArgumentException("Invalid dimensions");
            }
        }

        // Strength of selection matrix
        this.strengthOfSelectionMatrixParameter = strengthOfSelectionMatrixParam;
//        calculateStrengthOfSelectionInfo(diffusionModel);
        addVariable(strengthOfSelectionMatrixParameter);
        // Eigen decomposition
        this.eigenDecompositionStrengthOfSelection = decomposeStrenghtOfSelection(strengthOfSelectionMatrixParam);
    }

    private EigenDecomposition decomposeStrenghtOfSelection(MatrixParameterInterface Aparam){
        DenseMatrix64F A = MissingOps.wrap(Aparam);
        int n = A.numCols;
        // Checks
        if (n != A.numRows) throw new RuntimeException("Selection strength A matrix must be square.");
        if (!MatrixFeatures.isSymmetric(A)) throw new RuntimeException("Selection strength A matrix must be symmetric."); // TODO : this is not strictly necessary, but might be good to impose the constraint ?
        // Decomposition
        EigenDecomposition eigA = DecompositionFactory.eig(n, true, true);
        if( !eigA.decompose(A) ) throw new RuntimeException("Eigen decomposition failed.");
        return eigA;
    }

//    protected void calculateStrengthOfSelectionInfo(MultivariateDiffusionModel diffusionModel) {
//        strengthOfSelectionMatrix = strengthOfSelectionMatrixParameter.getParameterAsMatrix();
////        stationaryPrecisionMatrix = computeStationaryVariance(strengthOfSelectionMatrix, diffusionModel.getPrecisionParameter())
//    }

    public double[][] getStrengthOfSelection() {
        return strengthOfSelectionMatrixParameter.getParameterAsMatrix();
    }

    public double[] getEigenValuesStrengthOfSelection() {
        double[] eigA = new double[dim];
        for (int p = 0; p < dim; ++p) {
            Complex64F ev = eigenDecompositionStrengthOfSelection.getEigenvalue(p);
            if (!ev.isReal()) throw new RuntimeException("Selection strength A should only have real eigenvalues.");
            eigA[p] = ev.real;
        }
        return eigA;
    }

    public double[] getEigenVectorsStrengthOfSelection() {
        return EigenOps.createMatrixV(eigenDecompositionStrengthOfSelection).data;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

        if (branchRateModels.contains(model)) {
            fireModelChanged(model);
        } else {
            super.handleModelChangedEvent(model, object, index);
        }
    }

    @Override
    public boolean hasDrift() { return true; }

    @Override
    public boolean hasActualization() { return true; }

    @Override
    protected double[] getDriftRates(int[] branchIndices, int updateCount) {

        final double[] drift = new double[updateCount * dim];  // TODO Reuse?

        if (branchRateModels != null) {

            int offset = 0;
            for (int i = 0; i < updateCount; ++i) {

                final NodeRef node = tree.getNode(branchIndices[i]); // TODO Check if correct node

                for (int model = 0; model < dim; ++model) {
                    drift[offset] = branchRateModels.get(model).getBranchRate(tree, node);
                    ++offset;
                }
            }
        }
        return drift;
    }

    @Override
    public void setDiffusionModels(ContinuousDiffusionIntegrator cdi, boolean flip) {
        super.setDiffusionModels(cdi, flip);

        cdi.setDiffusionStationaryVariance(getEigenBufferOffsetIndex(0),
                getEigenValuesStrengthOfSelection(), getEigenVectorsStrengthOfSelection());
    }

    @Override
    public void updateDiffusionMatrices(ContinuousDiffusionIntegrator cdi, int[] branchIndices, double[] edgeLengths,
                                        int updateCount, boolean flip) {

        int[] probabilityIndices = new int[updateCount];

        for (int i = 0; i < updateCount; i++) {
            if (flip) {
                flipMatrixBufferOffset(branchIndices[i]);
            }
            probabilityIndices[i] = getMatrixBufferOffsetIndex(branchIndices[i]);
        }

        cdi.updateOrnsteinUhlenbeckDiffusionMatrices(
                getEigenBufferOffsetIndex(0),
                probabilityIndices,
                edgeLengths,
                getDriftRates(branchIndices, updateCount),
                getEigenValuesStrengthOfSelection(),
                getEigenVectorsStrengthOfSelection(),
                updateCount);
    }

    double[] getAccumulativeDrift(final NodeRef node, double[] priorMean) {
        final DenseMatrix64F drift = new DenseMatrix64F(dim, 1, true, priorMean);
        recursivelyAccumulateDrift(node, drift);
        return drift.data;
    }

    private void recursivelyAccumulateDrift(final NodeRef node, final DenseMatrix64F drift) {
        if (!tree.isRoot(node)) {

            // Compute parent
            recursivelyAccumulateDrift(tree.getParent(node), drift);

            // Actualize
            int[] branchIndice = new int[1];
            branchIndice[0] = getMatrixBufferOffsetIndex(node.getNumber());

            final double length = tree.getBranchLength(node);

            DenseMatrix64F actualization = new DenseMatrix64F(dim, dim);
            computeActualizationBranch(-length, actualization);

            DenseMatrix64F temp = new DenseMatrix64F(dim, 1);
            CommonOps.mult(actualization, drift, temp);
            CommonOps.scale(1.0, temp, drift);

            // Add optimal value
            DenseMatrix64F idMinusAct = CommonOps.identity(dim);
            CommonOps.addEquals(idMinusAct, -1.0, actualization);

            DenseMatrix64F optVal = new DenseMatrix64F(dim, 1, true, getDriftRates(branchIndice, 1));

            CommonOps.multAdd(idMinusAct, optVal, drift);

        }
    }

    private void computeActualizationBranch(double lambda, DenseMatrix64F C){
        DenseMatrix64F A = new DenseMatrix64F(strengthOfSelectionMatrixParameter.getParameterAsMatrix());
        EigenDecomposition eigA = DecompositionFactory.eig(dim, true);
        if( !eigA.decompose(A) ) throw new RuntimeException("Eigen decomposition failed.");
        DenseMatrix64F expDiag = CommonOps.identity(dim);
        for (int p = 0; p < dim; ++p) {
            Complex64F ev = eigA.getEigenvalue(p);
            if (!ev.isReal()) throw new RuntimeException("Selection strength A should only have real eigenvalues.");
            expDiag.set(p, p, Math.exp(lambda * ev.real));
        }
        DenseMatrix64F V = EigenOps.createMatrixV(eigA);
        DenseMatrix64F tmp = new DenseMatrix64F(dim, dim);
        CommonOps.mult(V, expDiag, tmp);
        CommonOps.invert(V);
        CommonOps.mult(tmp, V, C);
    }
}
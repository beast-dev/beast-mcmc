package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.evomodel.treedatalikelihood.continuous.cdi.SafeMultivariateDiagonalActualizedWithDriftIntegrator;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import dr.math.matrixAlgebra.missingData.MissingOps;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

final class OUBranchTimeGradient extends ContinuousTraitGradientForBranch.Default {

    private static final double ZERO_TOLERANCE = 1.0e-12;

    private final OUDiffusionModelDelegate diffusionProcessDelegate;
    private final ContinuousDiffusionIntegrator cdi;
    private final DenseMatrix64F matrix0;
    private final DenseMatrix64F matrix1;
    private final DenseMatrix64F vector0;
    private final double[] actualization;

    OUBranchTimeGradient(int dim,
                         Tree tree,
                         ContinuousDataLikelihoodDelegate likelihoodDelegate) {
        super(dim, tree);

        if (!(likelihoodDelegate.getDiffusionProcessDelegate() instanceof OUDiffusionModelDelegate)) {
            throw new IllegalArgumentException("OUBranchTimeGradient requires an OUDiffusionModelDelegate.");
        }

        this.diffusionProcessDelegate =
                (OUDiffusionModelDelegate) likelihoodDelegate.getDiffusionProcessDelegate();
        if (!diffusionProcessDelegate.hasDiagonalActualization()) {
            throw new IllegalArgumentException("OUBranchTimeGradient currently requires diagonal OU actualization.");
        }
        if (!(likelihoodDelegate.getIntegrator() instanceof SafeMultivariateDiagonalActualizedWithDriftIntegrator)) {
            throw new IllegalArgumentException(
                    "OUBranchTimeGradient requires SafeMultivariateDiagonalActualizedWithDriftIntegrator.");
        }

        this.cdi = likelihoodDelegate.getIntegrator();
        this.matrix0 = new DenseMatrix64F(dim, dim);
        this.matrix1 = new DenseMatrix64F(dim, dim);
        this.vector0 = new DenseMatrix64F(dim, 1);
        this.actualization = new double[dim];
    }

    @Override
    public int getDimension() {
        return 1;
    }

    @Override
    double[] chainRule(BranchSufficientStatistics statistics, NodeRef node,
                       DenseMatrix64F gradQInv, DenseMatrix64F gradN) {

        final double[] gradient = new double[1];
        gradient[0] += getGradientVarianceWrtBranchTime(statistics, node, gradQInv);
        gradient[0] += getGradientBranchVarianceWrtBranchTime(node, gradQInv);
        gradient[0] += getGradientDisplacementWrtBranchTime(statistics, node, gradN);

        return gradient;
    }

    @Override
    double[] chainRuleRoot(BranchSufficientStatistics statistics, NodeRef node,
                           DenseMatrix64F gradQInv, DenseMatrix64F gradN) {
        return new double[1];
    }

    private double getGradientVarianceWrtBranchTime(BranchSufficientStatistics statistics,
                                                    NodeRef node,
                                                    DenseMatrix64F gradient) {
        final int matrixIndex = diffusionProcessDelegate.getMatrixBufferOffsetIndex(node.getNumber());

        final DenseMatrix64F wi = statistics.getAbove().getRawVarianceCopy();

        final double[] branchVariance = matrix0.getData();
        cdi.getBranchVariance(matrixIndex,
                diffusionProcessDelegate.getEigenBufferOffsetIndex(0),
                branchVariance);
        final DenseMatrix64F sigma = MissingOps.wrap(branchVariance, 0, dim, dim);

        CommonOps.addEquals(wi, -1.0, sigma);
        CommonOps.multTransB(wi, gradient, matrix1);
        CommonOps.scale(2.0, matrix1);

        CommonOps.extractDiag(matrix1, vector0);
        cdi.getBranchActualization(matrixIndex, actualization);
        return chainActualizationGradientWrtBranchTime(vector0);
    }

    private double getGradientBranchVarianceWrtBranchTime(NodeRef node,
                                                          DenseMatrix64F gradient) {

        final int matrixIndex = diffusionProcessDelegate.getMatrixBufferOffsetIndex(node.getNumber());
        final int eigenIndex = diffusionProcessDelegate.getEigenBufferOffsetIndex(0);
        final double time = cdi.getBranchLength(matrixIndex);

        cdi.getBranchActualization(matrixIndex, actualization);

        final double[] attenuation = diffusionProcessDelegate.getEigenValuesStrengthOfSelection();
        final double[] stationaryVariance =
                ((SafeMultivariateDiagonalActualizedWithDriftIntegrator) cdi).getStationaryVariance(eigenIndex);

        final double[] branchVariance = matrix0.getData();
        cdi.getBranchVariance(matrixIndex, eigenIndex, branchVariance);

        double result = 0.0;
        for (int row = 0; row < dim; ++row) {
            for (int col = 0; col < dim; ++col) {
                final double attenuationSum = attenuation[row] + attenuation[col];
                final double stationary = stationaryVariance[row * dim + col];
                final double derivative;
                if (Math.abs(attenuationSum) <= ZERO_TOLERANCE || !Double.isFinite(stationary)) {
                    derivative = Math.abs(time) <= ZERO_TOLERANCE ? 0.0 : branchVariance[row * dim + col] / time;
                } else {
                    derivative = stationary * attenuationSum * actualization[row] * actualization[col];
                }
                result += derivative * gradient.unsafe_get(row, col);
            }
        }

        return result;
    }

    private double getGradientDisplacementWrtBranchTime(BranchSufficientStatistics statistics,
                                                        NodeRef node,
                                                        DenseMatrix64F gradient) {

        final DenseMatrix64F ni = statistics.getAbove().getRawMean();
        final DenseMatrix64F beta = MissingOps.wrap(diffusionProcessDelegate.getDriftRate(node), 0, dim, 1);

        CommonOps.add(ni, -1.0, beta, vector0);
        CommonOps.multTransB(gradient, vector0, matrix0);
        CommonOps.extractDiag(matrix0, vector0);

        cdi.getBranchActualization(
                diffusionProcessDelegate.getMatrixBufferOffsetIndex(node.getNumber()),
                actualization);
        return chainActualizationGradientWrtBranchTime(vector0);
    }

    private double chainActualizationGradientWrtBranchTime(DenseMatrix64F gradient) {
        final double[] attenuation = diffusionProcessDelegate.getEigenValuesStrengthOfSelection();

        double result = 0.0;
        for (int row = 0; row < dim; ++row) {
            result -= attenuation[row] * gradient.unsafe_get(row, 0);
        }
        return result;
    }
}

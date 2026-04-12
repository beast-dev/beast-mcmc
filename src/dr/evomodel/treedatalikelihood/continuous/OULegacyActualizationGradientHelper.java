package dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.continuous.MultivariateElasticModel;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.evomodel.treedatalikelihood.continuous.cdi.MultivariateIntegrator;
import dr.evomodel.treedatalikelihood.continuous.cdi.SafeMultivariateActualizedWithDriftIntegrator;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import static dr.math.matrixAlgebra.missingData.MissingOps.wrap;
import static dr.math.matrixAlgebra.missingData.MissingOps.diagMult;

public final class OULegacyActualizationGradientHelper {

    private OULegacyActualizationGradientHelper() {
    }

    public static void actualizeGradient(final MultivariateElasticModel elasticModel,
                                         final int dim,
                                         final ContinuousDiffusionIntegrator cdi,
                                         final int matrixBufferOffsetIndex,
                                         final int nodeIndex,
                                         final DenseMatrix64F gradient) {
        if (elasticModel.hasBlockStructure()) {
            throw new UnsupportedOperationException(
                    "Legacy eigen-basis actualization gradient is not valid for block-structured selection matrices; use the canonical block bridge path.");
        }
        final double[] attenuationRotation = elasticModel.getEigenVectorsStrengthOfSelection();
        final DenseMatrix64F p = wrap(attenuationRotation, 0, dim, dim);

        SafeMultivariateActualizedWithDriftIntegrator.transformMatrix(
                gradient, p, elasticModel.hasOrthogonalActualizationBasis());
        actualizeGradientDiagonal(elasticModel, dim, cdi, matrixBufferOffsetIndex, nodeIndex, gradient);
        SafeMultivariateActualizedWithDriftIntegrator.transformMatrixBack(gradient, p);
    }

    public static void actualizeGradientDiagonal(final MultivariateElasticModel elasticModel,
                                                 final int dim,
                                                 final ContinuousDiffusionIntegrator cdi,
                                                 final int matrixBufferOffsetIndex,
                                                 final int nodeIndex,
                                                 final DenseMatrix64F gradient) {
        if (elasticModel.hasBlockStructure()) {
            throw new UnsupportedOperationException(
                    "Legacy diagonal actualization gradient is not valid for block-structured selection matrices; use the canonical block bridge path.");
        }
        final double[] attenuation = elasticModel.getEigenValuesStrengthOfSelection();
        final double edgeLength = cdi.getBranchLength(matrixBufferOffsetIndex);
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                gradient.unsafe_set(i, j,
                        factorFunction(attenuation[i] + attenuation[j], edgeLength) * gradient.unsafe_get(i, j));
            }
        }
    }

    public static double factorFunction(final double x, final double l) {
        if (x == 0) return l;
        return -Math.expm1(-x * l) / x;
    }

    public static DenseMatrix64F getGradientVarianceWrtAttenuationDiagonal(final MultivariateElasticModel elasticModel,
                                                                           final int dim,
                                                                           final ContinuousDiffusionIntegrator cdi,
                                                                           final int matrixBufferOffsetIndex,
                                                                           final int eigenBufferOffsetIndex,
                                                                           final BranchSufficientStatistics statistics,
                                                                           final int nodeIndex,
                                                                           final DenseMatrix64F gradient) {
        final DenseMatrix64F gradActualization = getGradientVarianceWrtActualizationDiagonal(
                dim, cdi, matrixBufferOffsetIndex, eigenBufferOffsetIndex, statistics, nodeIndex, gradient);
        final DenseMatrix64F gradStationary = getGradientBranchVarianceWrtAttenuationDiagonal(
                elasticModel, dim, cdi, matrixBufferOffsetIndex, eigenBufferOffsetIndex, nodeIndex, gradient);
        CommonOps.addEquals(gradActualization, gradStationary);
        return gradActualization;
    }

    public static DenseMatrix64F getGradientVarianceWrtActualizationDiagonal(final int dim,
                                                                             final ContinuousDiffusionIntegrator cdi,
                                                                             final int matrixBufferOffsetIndex,
                                                                             final int eigenBufferOffsetIndex,
                                                                             final BranchSufficientStatistics statistics,
                                                                             final int nodeIndex,
                                                                             final DenseMatrix64F gradient) {
        final DenseMatrix64F wi = statistics.getAbove().getRawVarianceCopy();
        final double[] branchVariance = new double[dim * dim];
        cdi.getBranchVariance(matrixBufferOffsetIndex, eigenBufferOffsetIndex, branchVariance);
        final DenseMatrix64F sigmaI = wrap(branchVariance, 0, dim, dim);

        final DenseMatrix64F res = new DenseMatrix64F(dim, dim);
        CommonOps.addEquals(wi, -1, sigmaI);
        CommonOps.multTransB(wi, gradient, res);
        CommonOps.scale(2.0, res);

        final DenseMatrix64F resDiag = new DenseMatrix64F(dim, 1);
        CommonOps.extractDiag(res, resDiag);
        final double ti = cdi.getBranchLength(matrixBufferOffsetIndex);
        chainRuleActualizationWrtAttenuationDiagonal(ti, resDiag);
        return resDiag;
    }

    public static void chainRuleActualizationWrtAttenuationDiagonal(final double ti,
                                                                    final DenseMatrix64F grad) {
        CommonOps.scale(-ti, grad);
    }

    public static DenseMatrix64F getGradientBranchVarianceWrtAttenuationDiagonal(final MultivariateElasticModel elasticModel,
                                                                                 final int dim,
                                                                                 final ContinuousDiffusionIntegrator cdi,
                                                                                 final int matrixBufferOffsetIndex,
                                                                                 final int eigenBufferOffsetIndex,
                                                                                 final int nodeIndex,
                                                                                 final DenseMatrix64F gradient) {
        final double[] attenuation = elasticModel.getEigenValuesStrengthOfSelection();
        final DenseMatrix64F variance = wrap(
                ((MultivariateIntegrator) cdi).getVariance(eigenBufferOffsetIndex),
                0, dim, dim);

        final double ti = cdi.getBranchLength(matrixBufferOffsetIndex);
        final DenseMatrix64F res = new DenseMatrix64F(dim, 1);
        CommonOps.elementMult(variance, gradient);

        for (int k = 0; k < dim; k++) {
            double sum = 0.0;
            for (int l = 0; l < dim; l++) {
                sum -= variance.unsafe_get(k, l) *
                        computeAttenuationFactorActualized(attenuation[k] + attenuation[l], ti);
            }
            res.unsafe_set(k, 0, sum);
        }

        return res;
    }

    public static double computeAttenuationFactorActualized(final double lambda, final double ti) {
        if (lambda == 0) return ti * ti;
        final double em1 = Math.expm1(-lambda * ti);
        return 2.0 * (em1 * em1 - (em1 + lambda * ti) * Math.exp(-lambda * ti)) / lambda / lambda;
    }

    public static DenseMatrix64F getGradientDisplacementWrtAttenuationDiagonal(final int dim,
                                                                               final ContinuousDiffusionIntegrator cdi,
                                                                               final int matrixBufferOffsetIndex,
                                                                               final BranchSufficientStatistics statistics,
                                                                               final double[] driftRate,
                                                                               final DenseMatrix64F gradient) {
        final DenseMatrix64F ni = statistics.getAbove().getRawMean();
        final DenseMatrix64F betai = wrap(driftRate, 0, dim, 1);

        final DenseMatrix64F resFull = new DenseMatrix64F(dim, dim);
        final DenseMatrix64F resDiag = new DenseMatrix64F(dim, 1);
        CommonOps.add(ni, -1, betai, resDiag);
        CommonOps.multTransB(gradient, resDiag, resFull);
        CommonOps.extractDiag(resFull, resDiag);

        final double ti = cdi.getBranchLength(matrixBufferOffsetIndex);
        chainRuleActualizationWrtAttenuationDiagonal(ti, resDiag);
        return resDiag;
    }

    public static void actualizeDisplacementGradientDiagonal(final int dim,
                                                             final ContinuousDiffusionIntegrator cdi,
                                                             final int matrixBufferOffsetIndex,
                                                             final DenseMatrix64F gradient) {
        final double[] qi = new double[dim];
        cdi.getBranch1mActualization(matrixBufferOffsetIndex, qi);
        diagMult(qi, gradient);
    }

    public static void actualizeDisplacementGradient(final int dim,
                                                     final ContinuousDiffusionIntegrator cdi,
                                                     final int matrixBufferOffsetIndex,
                                                     final DenseMatrix64F gradient) {
        final double[] qi = new double[dim * dim];
        cdi.getBranch1mActualization(matrixBufferOffsetIndex, qi);
        final DenseMatrix64F actu = wrap(qi, 0, dim, dim);
        CommonOps.scale(-1.0, actu);
        final DenseMatrix64F tmp = new DenseMatrix64F(dim, 1);
        CommonOps.multTransA(actu, gradient, tmp);
        CommonOps.scale(-1.0, tmp, gradient);
    }

    public static double[] actualizeRootGradientDiagonal(final int dim,
                                                         final ContinuousDiffusionIntegrator cdi,
                                                         final int matrixBufferOffsetIndex,
                                                         final DenseMatrix64F gradient) {
        final double[] qi = new double[dim];
        cdi.getBranchActualization(matrixBufferOffsetIndex, qi);
        final DenseMatrix64F tmp = new DenseMatrix64F(dim, 1);
        diagMult(qi, gradient, tmp);
        return tmp.getData();
    }

    public static double[] actualizeRootGradientFull(final int dim,
                                                     final ContinuousDiffusionIntegrator cdi,
                                                     final int matrixBufferOffsetIndex,
                                                     final DenseMatrix64F gradient) {
        final double[] qi = new double[dim * dim];
        cdi.getBranchActualization(matrixBufferOffsetIndex, qi);
        final DenseMatrix64F actu = wrap(qi, 0, dim, dim);
        final DenseMatrix64F tmp = new DenseMatrix64F(dim, 1);
        CommonOps.multTransA(actu, gradient, tmp);
        return tmp.getData();
    }

}

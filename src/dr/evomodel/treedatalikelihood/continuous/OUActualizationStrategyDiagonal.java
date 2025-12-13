package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.NodeRef;
import dr.evomodel.continuous.MultivariateElasticModel;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.evomodel.treedatalikelihood.continuous.cdi.MultivariateIntegrator;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import dr.math.matrixAlgebra.missingData.MissingOps;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import static dr.math.matrixAlgebra.missingData.MissingOps.wrap;

/**
 * Diagonal-S strategy: uses exact elementwise formulas and never touches
 * Fréchet / block-Schur machinery. This gathers all the diagonal-only
 * branches that previously lived inside OUDiffusionModelDelegateNewest.
 *
 * Assumes:
 *  - elasticModel.isDiagonal() is true,
 *  - CDI actualization is diagonal in the chosen basis.
 */
public final class OUActualizationStrategyDiagonal implements OUActualizationStrategy {

    private final OUDelegateContext ctx;
    private final int dim;
    private final MultivariateElasticModel elasticModel;

    public OUActualizationStrategyDiagonal(OUDelegateContext ctx) {
        this.ctx = ctx;
        this.dim = ctx.getDim();
        this.elasticModel = ctx.getElasticModel();
    }

    // ------------------------------------------------------------------
    // Variance wrt variance parameters
    // ------------------------------------------------------------------

    @Override
    public DenseMatrix64F gradientVarianceWrtVariance(NodeRef node,
                                                      ContinuousDiffusionIntegrator cdi,
                                                      ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                      DenseMatrix64F gradient) {
        // Delegate has already guarded root; we can assume non-root here.
        DenseMatrix64F result = gradient.copy();
        actualizeGradientDiagonal(cdi, node.getNumber(), result);
        return result;
    }

    // ------------------------------------------------------------------
    // Variance wrt attenuation (λ)
    // ------------------------------------------------------------------

    @Override
    public DenseMatrix64F gradientVarianceWrtAttenuation(NodeRef node,
                                                         ContinuousDiffusionIntegrator cdi,
                                                         BranchSufficientStatistics statistics,
                                                         DenseMatrix64F dSigma) {
        final int nodeIndex = node.getNumber();
        // Exact diagonal path: actualization + stationary-variance contributions.
        return getGradientVarianceWrtAttenuationDiagonalExact(
                cdi, statistics, nodeIndex, dSigma);
    }

    // ------------------------------------------------------------------
    // Displacement wrt drift and attenuation
    // ------------------------------------------------------------------

    @Override
    public DenseMatrix64F gradientDisplacementWrtDrift(NodeRef node,
                                                       ContinuousDiffusionIntegrator cdi,
                                                       DenseMatrix64F gradient) {
        DenseMatrix64F result = gradient.copy();
        actualizeDisplacementGradientDiagonal(cdi, node.getNumber(), result);
        return result;
    }

    @Override
    public DenseMatrix64F gradientDisplacementWrtAttenuation(NodeRef node,
                                                             ContinuousDiffusionIntegrator cdi,
                                                             BranchSufficientStatistics statistics,
                                                             DenseMatrix64F gradient) {
        final int nodeIndex = node.getNumber();
        final double ti = cdi.getBranchLength(ctx.getMatrixBufferOffsetIndex(nodeIndex));

        DenseMatrix64F ni    = statistics.getAbove().getRawMean().copy();      // d × 1
        DenseMatrix64F betai = wrap(ctx.getDriftRate(node), 0, dim, 1);        // d × 1

        boolean gradientIsVector = (gradient.getNumCols() == 1);

        DenseMatrix64F diff = new DenseMatrix64F(dim, 1);
        CommonOps.add(ni, -1.0, betai, diff);

        DenseMatrix64F resDiag = new DenseMatrix64F(dim, 1);
        if (gradientIsVector) {
            for (int k = 0; k < dim; k++) {
                resDiag.unsafe_set(
                        k, 0,
                        -ti * gradient.unsafe_get(k, 0) * diff.unsafe_get(k, 0)
                );
            }
        } else {
            DenseMatrix64F resFull = new DenseMatrix64F(dim, dim);
            CommonOps.multTransB(gradient, diff, resFull);
            CommonOps.extractDiag(resFull, resDiag);
            CommonOps.scale(-ti, resDiag);
        }
        return resDiag;
    }

    // ------------------------------------------------------------------
    // Heritability-style summaries
    // ------------------------------------------------------------------

    @Override
    public void meanTipVariances(double priorSampleSize,
                                 double[] treeLengths,
                                 DenseMatrix64F traitVariance,
                                 DenseMatrix64F varSum) {
        getMeanTipVariancesDiagonal(priorSampleSize, treeLengths, traitVariance, varSum);
    }

    @Override
    public double[] rootGradient(int index, ContinuousDiffusionIntegrator cdi, DenseMatrix64F gradient) {
        double[] ai = new double[dim];
//        private double[] rootGradient(ContinuousDiffusionIntegrator cdi, int nodeindex, DenseMatrix64F gradient) {
//      TODO the first entry should contain getMatrixBufferOffsetIndex(nodeIndex) --> nodeindex should be already transformed
        cdi.getBranchActualization(index, ai);
        DenseMatrix64F tmp = new DenseMatrix64F(dim, 1);
        MissingOps.diagMult(ai, gradient, tmp);
        return tmp.getData();
    }

    @Override
    public DenseMatrix64F gradientVarianceDirectBackprop(NodeRef node, ContinuousDiffusionIntegrator cdi, BranchSufficientStatistics statistics, DenseMatrix64F dL_dJ, DenseMatrix64F dL_deta, double dL_dc) {
        return null;
    }

    @Override
    public DenseMatrix64F gradientDisplacementDirectBackprop(NodeRef node, ContinuousDiffusionIntegrator cdi, BranchSufficientStatistics statistics, DenseMatrix64F dL_dJ, DenseMatrix64F dL_deta, double dL_dc) {
        return null;
    }


    // ------------------------------------------------------------------
    // Diagonal-only helpers
    // ------------------------------------------------------------------

    /**
     * Variance wrt variance parameters: apply elementwise attenuation factor
     * f(λ_i + λ_j, t) to each gradient entry.
     */
    private void actualizeGradientDiagonal(ContinuousDiffusionIntegrator cdi,
                                           int nodeIndex,
                                           DenseMatrix64F gradient) {
        double[] lam = elasticModel.getBasisEigenValues();
        double t = cdi.getBranchLength(ctx.getMatrixBufferOffsetIndex(nodeIndex));

        for (int i = 0; i < dim; i++) {
            double li = lam[i];
            for (int j = 0; j < dim; j++) {
                double lj = lam[j];
                double x  = li + lj;
                double g  = gradient.unsafe_get(i, j);
                gradient.unsafe_set(i, j, factorFunction(x, t) * g);
            }
        }
    }

    /**
     * Diagonal S, exact path:
     *   dL/dλ = dL/dλ (actualization) + dL/dλ (stationary variance).
     */
    private DenseMatrix64F getGradientVarianceWrtAttenuationDiagonalExact(
            ContinuousDiffusionIntegrator cdi,
            BranchSufficientStatistics statistics,
            int nodeIndex,
            DenseMatrix64F gradient /* ∂L/∂Σ_i in eigen/diag basis */) {

        // 1) actualization part
        DenseMatrix64F gradActualization =
                getGradientVarianceWrtActualizationDiagonalExact(
                        cdi, statistics, nodeIndex, gradient);

        // 2) stationary-variance part
        DenseMatrix64F gradStationary =
                getGradientBranchVarianceWrtAttenuationDiagonalExact(
                        cdi, nodeIndex, gradient);

        CommonOps.addEquals(gradActualization, gradStationary);
        return gradActualization;
    }

    /**
     * Diagonal S: exact actualization contribution dL/dλ.
     * Works entirely in the eigen/diag basis, no expm/Fréchet.
     */
    private DenseMatrix64F getGradientVarianceWrtActualizationDiagonalExact(
            ContinuousDiffusionIntegrator cdi,
            BranchSufficientStatistics statistics,
            int nodeIndex,
            DenseMatrix64F gradient /* ∂L/∂Σ_parent in eigen basis */) {

        final int d = dim;

        // W_i = above variance (in working/eigen basis in diagonal mode)
        DenseMatrix64F Wi = statistics.getAbove().getRawVarianceCopy();

        // Σ_i = branch/process variance (same basis)
        double[] branchVariance = new double[d * d];
        cdi.getBranchVariance(ctx.getMatrixBufferOffsetIndex(nodeIndex),
                ctx.getEigenBufferOffsetIndex(0),
                branchVariance);
        DenseMatrix64F Sigma_i = wrap(branchVariance, 0, d, d);

        // Wdiff = Wi - Σ_i  (reuse Wi buffer)
        DenseMatrix64F Wdiff = Wi;
        CommonOps.addEquals(Wdiff, -1.0, Sigma_i);

        // res = 2 * (Wi - Σ_i)^T * G
        DenseMatrix64F res = new DenseMatrix64F(d, d);
        CommonOps.multTransB(Wdiff, gradient, res);
        CommonOps.scale(2.0, res);

        // extract diagonal → components wrt λ_k
        DenseMatrix64F resDiag = new DenseMatrix64F(d, 1);
        CommonOps.extractDiag(res, resDiag);

        // chain rule ∂q/∂λ = -t * exp(-λ t); in diagonal mode q is implicit via A
        double ti = cdi.getBranchLength(ctx.getMatrixBufferOffsetIndex(nodeIndex));
        CommonOps.scale(-ti, resDiag);

        return resDiag;
    }

    /**
     * Diagonal S: exact stationary-variance contribution dL/dλ, using elementwise formulas
     * involving Γ and the attenuation eigenvalues.
     */
    private DenseMatrix64F getGradientBranchVarianceWrtAttenuationDiagonalExact(
            ContinuousDiffusionIntegrator cdi,
            int nodeIndex,
            DenseMatrix64F gradient /* ∂L/∂Σ_i in eigen/diag basis */) {

        final int d = dim;

        double[] attenuation = elasticModel.getBasisEigenValues(); // λ_k
        DenseMatrix64F variance = wrap(
                ((MultivariateIntegrator) cdi).getVariance(ctx.getEigenBufferOffsetIndex(0)),
                0, d, d);

        double ti = cdi.getBranchLength(ctx.getMatrixBufferOffsetIndex(nodeIndex));

        DenseMatrix64F res = new DenseMatrix64F(d, 1);

        // H = Γ .* gradient  (EJML mutates first arg)
        DenseMatrix64F H = variance.copy();
        CommonOps.elementMult(H, gradient);

        for (int k = 0; k < d; k++) {
            double sum = 0.0;
            for (int l = 0; l < d; l++) {
                double lamSum = attenuation[k] + attenuation[l];
                sum -= H.unsafe_get(k, l) * computeAttenuationFactorActualized(lamSum, ti);
            }
            res.unsafe_set(k, 0, sum);
        }

        return res;
    }

    /**
     * Displacement gradient wrt drift, diagonal S:
     * multiply each component by (1 - A_kk) coming from CDI.
     */
    private void actualizeDisplacementGradientDiagonal(ContinuousDiffusionIntegrator cdi,
                                                       int nodeIndex,
                                                       DenseMatrix64F gradient) {
        double[] qi = new double[dim];
        cdi.getBranch1mActualization(ctx.getMatrixBufferOffsetIndex(nodeIndex), qi);
        MissingOps.diagMult(qi, gradient);
    }

    /**
     * Heritability-style mean tip variances in diagonal basis: exact elementwise formula.
     */
    private void getMeanTipVariancesDiagonal(final double priorSampleSize,
                                             final double[] treeLengths,
                                             final DenseMatrix64F traitVariance,
                                             final DenseMatrix64F varSum) {
        double[] eigVals = elasticModel.getBasisEigenValues();
        int ntaxa = ctx.getTree().getExternalNodeCount();

        for (int i = 0; i < ntaxa; ++i) {
            double ti = treeLengths[i];
            for (int p = 0; p < dim; ++p) {
                double ep = eigVals[p];
                for (int q = 0; q < dim; ++q) {
                    double eq = eigVals[q];
                    double sum = ep + eq;
                    double var = (sum == 0.0)
                            ? (ti + 1.0 / priorSampleSize) * traitVariance.get(p, q)
                            : Math.exp(-sum * ti)
                            * (Math.expm1(sum * ti) / sum + 1.0 / priorSampleSize)
                            * traitVariance.get(p, q);
                    varSum.set(p, q, varSum.get(p, q) + var);
                }
            }
        }
        CommonOps.scale(1.0 / treeLengths.length, varSum);
    }

    // ------------------------------------------------------------------
    // Scalar helpers
    // ------------------------------------------------------------------

    /** f(x, t) = (1 - e^{-xt}) / x with good behaviour when x → 0. */
    private static double factorFunction(double x, double t) {
        if (x == 0.0) return t;
        return -Math.expm1(-x * t) / x;
    }

    /**
     * Attenuation factor used in the stationary-variance gradient.
     * Same formula as in the delegate, duplicated here for independence.
     */
    private static double computeAttenuationFactorActualized(double lambda, double ti) {
        if (lambda == 0.0) return ti * ti;
        double em1 = Math.expm1(-lambda * ti);
        double em  = Math.exp(-lambda * ti);
        return 2.0 * (em1 * em1 - (em1 + lambda * ti) * em) / (lambda * lambda);
    }
}

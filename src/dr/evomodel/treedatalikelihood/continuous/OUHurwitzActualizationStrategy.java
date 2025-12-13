package dr.evomodel.treedatalikelihood.continuous;

import cern.colt.matrix.DoubleMatrix2D;
import cern.colt.matrix.impl.DenseDoubleMatrix2D;
import cern.colt.matrix.linalg.Algebra;
import dr.evolution.tree.NodeRef;
import dr.evomodel.continuous.MultivariateElasticModel;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.evomodel.treedatalikelihood.continuous.cdi.MultivariateIntegrator;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import static dr.math.matrixAlgebra.missingData.MissingOps.wrap;

/**
 * OUActualizationStrategy for the Hurwitz (real-diagonalizable) case.
 *
 * Assumes S is real diagonalizable: S = R Λ R^{-1} with real eigenvalues Λ and
 * (possibly non-orthogonal) right eigenvectors R.
 *
 * This is the non-diagonal path extracted from
 * OUDiffusionModelDelegateWorkingWithEnum, reshaped into the
 * OUActualizationStrategy interface.
 */
public final class OUHurwitzActualizationStrategy implements OUActualizationStrategy {

    private final OUDelegateContext ctx;
    private final int dim;
    private final MultivariateElasticModel elasticModel;

    // Basis S = R Λ R^{-1}
    private double[] evals;
    private DenseMatrix64F R;
    private DenseMatrix64F Rinv;

    // Small temporaries
    private final DenseMatrix64F tmpDxD;
    private final DenseMatrix64F tmpDx1;

    public OUHurwitzActualizationStrategy(OUDelegateContext ctx) {
        this.ctx = ctx;
        this.dim = ctx.getDim();
        this.elasticModel = ctx.getElasticModel();

        this.tmpDxD = new DenseMatrix64F(dim, dim);
        this.tmpDx1 = new DenseMatrix64F(dim, 1);

        refreshBasis();
    }

    /**
     * Call this from the delegate when the elastic model changes.
     */
    public void refreshBasis() {
        this.evals = elasticModel.getBasisEigenValues();     // real spectrum
        double[] evecsFlat = elasticModel.getBasisEigenVectors(); // right eigenvectors
        this.R = wrap(evecsFlat, 0, dim, dim);
        this.Rinv = new DenseMatrix64F(dim, dim);
        CommonOps.invert(R, Rinv);
    }

    // -------------------------------------------------------------------------
    // OUActualizationStrategy interface
    // -------------------------------------------------------------------------

    @Override
    public DenseMatrix64F gradientVarianceWrtVariance(NodeRef node,
                                                      ContinuousDiffusionIntegrator cdi,
                                                      ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                      DenseMatrix64F gradient) {
        // Delegate guarantees non-root here.
        DenseMatrix64F result = gradient.copy();
        actualizeVarianceGradient(cdi, node.getNumber(), result);
        return result;
    }

    @Override
    public DenseMatrix64F gradientVarianceWrtAttenuation(NodeRef node,
                                                         ContinuousDiffusionIntegrator cdi,
                                                         BranchSufficientStatistics statistics,
                                                         DenseMatrix64F dSigma) {
        final int d = dim;
        final int nodeIndex = node.getNumber();

        // (1) Actualization contribution: ∂ℓ/∂S from A = exp(-t S)
        DenseMatrix64F GS_from_A =
                getGradientVarianceWrtActualizationGeneral(cdi, statistics, nodeIndex, dSigma);

        // (2) Push ∂ℓ/∂S to (R, Λ) with S = R Λ R^{-1}
        DenseMatrix64F H = new DenseMatrix64F(d, d); // H = R^{-1} G_S R
        {
            DenseMatrix64F tmp = new DenseMatrix64F(d, d);
            CommonOps.mult(Rinv, GS_from_A, tmp);
            CommonOps.mult(tmp, R, H);
        }

        // (2a) eigenvalue gradient from the A-path: diag(H)
        DenseMatrix64F dLambda_from_A = new DenseMatrix64F(d, 1);
        for (int i = 0; i < d; i++) {
            dLambda_from_A.set(i, 0, H.get(i, i));
        }

        // (2b) eigenvector gradient from the A-path: G_R = R^{-T} [Λ, H^T]
        double[] lam = evals;
        DenseMatrix64F HT = H.copy();
        CommonOps.transpose(HT);

        DenseMatrix64F DHT = new DenseMatrix64F(d, d);
        DenseMatrix64F HTD = new DenseMatrix64F(d, d);

        // DHT = Λ H^T
        for (int i = 0; i < d; i++) {
            for (int j = 0; j < d; j++) {
                DHT.set(i, j, lam[i] * HT.get(i, j));
            }
        }
        // HTD = H^T Λ
        for (int i = 0; i < d; i++) {
            for (int j = 0; j < d; j++) {
                HTD.set(i, j, HT.get(i, j) * lam[j]);
            }
        }

        DenseMatrix64F comm = new DenseMatrix64F(d, d);
        CommonOps.subtract(DHT, HTD, comm); // [Λ, H^T]

        DenseMatrix64F RinvT = new DenseMatrix64F(d, d);
        CommonOps.transpose(Rinv, RinvT);
        DenseMatrix64F gradR_from_A = new DenseMatrix64F(d, d);
        CommonOps.mult(RinvT, comm, gradR_from_A);

        // (3) Branch-variance contribution (only depends on eigenvalues)
        DenseMatrix64F Geig_tmp = new DenseMatrix64F(d, d);
        DenseMatrix64F Geig = new DenseMatrix64F(d, d);
        CommonOps.multTransA(dSigma, Rinv, Geig_tmp); // dSigma^T Rinv = dSigma Rinv (sym)
        CommonOps.mult(RinvT, Geig_tmp, Geig);        // R^{-T} (dSigma R^{-1})

        DenseMatrix64F dLambda_from_Var =
                getGradientBranchVarianceWrtAttenuationDiagonal(cdi, nodeIndex, Geig);

        // (4) Return concatenated [ dΛ ; first two entries of dR ]
        DenseMatrix64F dLambda_total = dLambda_from_A.copy();
        CommonOps.addEquals(dLambda_total, dLambda_from_Var);

        DenseMatrix64F result = new DenseMatrix64F(d + 2, 1);
        for (int i = 0; i < d; i++) {
            result.set(i, 0, dLambda_total.get(i, 0));
        }
        result.set(d,     0, gradR_from_A.get(0, 0));
        result.set(d + 1, 0, gradR_from_A.get(0, 1)); //TODO generalize
        return result;
    }

    @Override
    public DenseMatrix64F gradientDisplacementWrtDrift(NodeRef node,
                                                       ContinuousDiffusionIntegrator cdi,
                                                       DenseMatrix64F gradient) {
        DenseMatrix64F result = gradient.copy();
        actualizeDisplacementGradient(cdi, node.getNumber(), result);
        return result;
    }

    @Override
    public DenseMatrix64F gradientDisplacementWrtAttenuation(NodeRef node,
                                                             ContinuousDiffusionIntegrator cdi,
                                                             BranchSufficientStatistics statistics,
                                                             DenseMatrix64F gradient) {
        final int nodeIndex = node.getNumber();
        final double ti = cdi.getBranchLength(ctx.getMatrixBufferOffsetIndex(nodeIndex));
        final boolean gradientIsVector = (gradient.getNumCols() == 1);

        DenseMatrix64F ni    = statistics.getAbove().getRawMean().copy();   // d x 1
        DenseMatrix64F betai = wrap(ctx.getDriftRate(node), 0, dim, 1);     // d x 1

        // Rotate (n_i, beta_i) to eigenbasis
        DenseMatrix64F niEig    = new DenseMatrix64F(dim, 1);
        DenseMatrix64F betaiEig = new DenseMatrix64F(dim, 1);
        CommonOps.mult(Rinv, ni,    niEig);
        CommonOps.mult(Rinv, betai, betaiEig);

        DenseMatrix64F diffEig = new DenseMatrix64F(dim, 1);
        CommonOps.add(niEig, -1.0, betaiEig, diffEig);

        // First compute the λ-part: dim × 1
        DenseMatrix64F resDiag = new DenseMatrix64F(dim, 1);

        if (gradientIsVector) {
            DenseMatrix64F mEig = new DenseMatrix64F(dim, 1);
            CommonOps.mult(Rinv, gradient, mEig);
            for (int k = 0; k < dim; k++) {
                resDiag.unsafe_set(
                        k, 0,
                        -ti * mEig.unsafe_get(k, 0) * diffEig.unsafe_get(k, 0)
                );
            }
        } else {
            DenseMatrix64F G = gradient.copy();
            transformMatrixToEig(G, R, Rinv);
            DenseMatrix64F resFull = new DenseMatrix64F(dim, dim);
            CommonOps.multTransB(G, diffEig, resFull);
            CommonOps.extractDiag(resFull, resDiag);
            CommonOps.scale(-ti, resDiag);
        }

        // Now extend to match the Hurwitz parameter vector length (dim + 2):
        // [ dℓ/dλ ; 0 ; 0 ]  – displacement does not depend on eigenvector params.
        DenseMatrix64F extended = new DenseMatrix64F(dim + 2, 1);
        for (int i = 0; i < dim; i++) {
            extended.set(i, 0, resDiag.get(i, 0));
        }
        // extended[d] and extended[d+1] stay 0.0

        return extended;
    }

    @Override
    public void meanTipVariances(double priorSampleSize,
                                 double[] treeLengths,
                                 DenseMatrix64F traitVariance,
                                 DenseMatrix64F varSum) {
        getMeanTipVariancesFull(priorSampleSize, treeLengths, traitVariance, varSum);
    }

    @Override
    public double[] rootGradient(int index,
                                 ContinuousDiffusionIntegrator cdi,
                                 DenseMatrix64F gradient) {
        // Here we always apply the *full* actualization A to the gradient:
        double[] aFlat = new double[dim * dim];
        cdi.getBranchActualization(ctx.getMatrixBufferOffsetIndex(index), aFlat);
        DenseMatrix64F A = wrap(aFlat, 0, dim, dim);
        DenseMatrix64F tmp = new DenseMatrix64F(dim, 1);
        CommonOps.mult(A, gradient, tmp);
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

    // -------------------------------------------------------------------------
    // Non-diagonal helpers (Hurwitz case)
    // -------------------------------------------------------------------------

    private void actualizeVarianceGradient(ContinuousDiffusionIntegrator cdi,
                                           int nodeIndex,
                                           DenseMatrix64F grad) {
        // Rotate to eigenbasis, apply diagonal factor, rotate back
        transformMatrixToEig(grad, R, Rinv);
        actualizeGradientDiagonal(cdi, nodeIndex, grad);
        transformMatrixFromEig(grad, R);
    }

    private void actualizeGradientDiagonal(ContinuousDiffusionIntegrator cdi,
                                           int nodeIndex,
                                           DenseMatrix64F gradientEig) {
        double[] lam = evals;
        double t = cdi.getBranchLength(ctx.getMatrixBufferOffsetIndex(nodeIndex));
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                double x = lam[i] + lam[j];
                gradientEig.unsafe_set(i, j,
                        factorFunction(x, t) * gradientEig.unsafe_get(i, j));
            }
        }
    }

    private static double factorFunction(double x, double t) {
        if (x == 0.0) return t;
        return -Math.expm1(-x * t) / x; // (1 - e^{-xt})/x, numerically stable
    }

    /**
     * ∂ℓ/∂S contribution coming from A = exp(-t S), via Fréchet adjoint.
     * All matrices are in the *working* (original) basis.
     */
    private DenseMatrix64F getGradientVarianceWrtActualizationGeneral(
            ContinuousDiffusionIntegrator cdi,
            BranchSufficientStatistics statistics,
            int nodeIndex,
            DenseMatrix64F G /* = gradient in working basis */) {

        final int d = dim;

        // A = exp(-t S)
        DenseMatrix64F A = new DenseMatrix64F(d, d);
        cdi.getBranchActualization(ctx.getMatrixBufferOffsetIndex(nodeIndex), A.getData());

        // Σ_child
        DenseMatrix64F Wchild = statistics.getAbove().getRawVarianceCopy();

        // Σ_i (branch/process variance)
        double[] branchVar = new double[d * d];
        cdi.getBranchVariance(ctx.getMatrixBufferOffsetIndex(nodeIndex),
                ctx.getEigenBufferOffsetIndex(0),
                branchVar);
        DenseMatrix64F Sigma_i = wrap(branchVar, 0, d, d);

        // A^{-T}
        DenseMatrix64F Ai = A.copy();
        CommonOps.invert(Ai);
        DenseMatrix64F AiT = new DenseMatrix64F(d, d);
        CommonOps.transpose(Ai, AiT);

        DenseMatrix64F Wdiff = Wchild.copy();
        CommonOps.addEquals(Wdiff, -1.0, Sigma_i);
        DenseMatrix64F AW = new DenseMatrix64F(d, d);
        CommonOps.mult(Wdiff, AiT, AW); // AW = (Σ_child - Σ_i) A^{-T} = A W_parent

        // M = G * AW; Msym = M + M^T
        DenseMatrix64F M = new DenseMatrix64F(d, d);
        CommonOps.mult(G, AW, M);
        DenseMatrix64F Msym = M.copy();
        DenseMatrix64F Mt   = M.copy();
        CommonOps.transpose(Mt);
        CommonOps.addEquals(Msym, Mt);

        // Fréchet adjoint at −t S via block exponential
        double ti = cdi.getBranchLength(ctx.getMatrixBufferOffsetIndex(nodeIndex));
        DenseMatrix64F S = wrap(elasticModel.getStrengthOfSelectionMatrixAsVector(), 0, d, d);
        DenseMatrix64F minus_tS = S.copy();
        CommonOps.scale(-ti, minus_tS);

        DenseMatrix64F K = new DenseMatrix64F(d, d);
        frechetAdjExpBlock(minus_tS, Msym, K);  // K = L_exp(−tS)^*[Msym]

        DenseMatrix64F GS_from_A = new DenseMatrix64F(d, d);
        CommonOps.scale(-ti, K, GS_from_A);
        return GS_from_A;
    }

    private DenseMatrix64F getGradientBranchVarianceWrtAttenuationDiagonal(ContinuousDiffusionIntegrator cdi,
                                                                           int nodeIndex,
                                                                           DenseMatrix64F gradientEig) {
        double[] lam = evals;
        DenseMatrix64F variance =
                wrap(((MultivariateIntegrator) cdi).getVariance(ctx.getEigenBufferOffsetIndex(0)),
                        0, dim, dim);
        double ti = cdi.getBranchLength(ctx.getMatrixBufferOffsetIndex(nodeIndex));

        DenseMatrix64F res = new DenseMatrix64F(dim, 1);

        DenseMatrix64F H = gradientEig.copy(); // H ← Γ .* gradientEig
        CommonOps.elementMult(H, variance);

        for (int k = 0; k < dim; k++) {
            double sum = 0.0;
            for (int l = 0; l < dim; l++) {
                double ls = lam[k] + lam[l];
                sum -= H.unsafe_get(k, l) * computeAttenuationFactorActualized(ls, ti);
            }
            res.unsafe_set(k, 0, sum);
        }
        return res;
    }

    private double computeAttenuationFactorActualized(double lambda, double ti) {
        if (lambda == 0.0) return ti * ti;
        double em1 = Math.expm1(-lambda * ti);
        double em  = Math.exp(-lambda * ti);
        return 2.0 * (em1 * em1 - (em1 + lambda * ti) * em) / (lambda * lambda);
    }

    private void actualizeDisplacementGradient(ContinuousDiffusionIntegrator cdi,
                                               int nodeIndex,
                                               DenseMatrix64F gradient) {
        double[] qFlat = new double[dim * dim];
        cdi.getBranch1mActualization(ctx.getMatrixBufferOffsetIndex(nodeIndex), qFlat);
        DenseMatrix64F Actu = wrap(qFlat, 0, dim, dim);
        CommonOps.scale(-1.0, Actu); // −(I − A)
        DenseMatrix64F tmp = new DenseMatrix64F(dim, 1);
        CommonOps.mult(Actu, gradient, tmp);
        CommonOps.scale(-1.0, tmp, gradient); // (I − A) gradient
    }

    // -------------------------------------------------------------------------
    // Heritability-style summaries (full / non-diagonal)
    // -------------------------------------------------------------------------

    private void getMeanTipVariancesFull(double priorSampleSize,
                                         double[] treeLengths,
                                         DenseMatrix64F traitVariance,
                                         DenseMatrix64F varSum) {
        DenseMatrix64F transV = traitVariance.copy();
        transformMatrixToEig(transV, R, Rinv);

        getMeanTipVariancesDiagonal(priorSampleSize, treeLengths, transV, varSum);

        DenseMatrix64F tmp = new DenseMatrix64F(dim, dim);
        CommonOps.mult(R, varSum, tmp);
        CommonOps.multTransB(tmp, R, varSum);
    }

    private void getMeanTipVariancesDiagonal(double priorSampleSize,
                                             double[] treeLengths,
                                             DenseMatrix64F traitVariance,
                                             DenseMatrix64F varSum) {
        double[] eigVals = evals;
        int ntaxa = treeLengths.length;

        for (int i = 0; i < ntaxa; ++i) {
            double ti = treeLengths[i];
            for (int p = 0; p < dim; ++p) {
                double ep = eigVals[p];
                for (int q = 0; q < dim; ++q) {
                    double eq = eigVals[q];
                    double sum = ep + eq;
                    double var = (sum == 0.0)
                            ? (ti + 1.0 / priorSampleSize) * traitVariance.get(p, q)
                            : Math.exp(-sum * ti) *
                            (Math.expm1(sum * ti) / sum + 1.0 / priorSampleSize) *
                            traitVariance.get(p, q);
                    varSum.set(p, q, varSum.get(p, q) + var);
                }
            }
        }
        CommonOps.scale(1.0 / ntaxa, varSum);
    }

    // -------------------------------------------------------------------------
    // Basis helpers
    // -------------------------------------------------------------------------

    private void transformMatrixToEig(DenseMatrix64F X,
                                      DenseMatrix64F R,
                                      DenseMatrix64F Rinv) {
        // X <- R^{-1} X R^{-T}
        CommonOps.mult(Rinv, X, tmpDxD);
        CommonOps.multTransB(tmpDxD, Rinv, X);
    }

    private void transformMatrixFromEig(DenseMatrix64F Xeig,
                                        DenseMatrix64F R) {
        // Xeig <- R Xeig R^T
        CommonOps.mult(R, Xeig, tmpDxD);
        CommonOps.multTransB(tmpDxD, R, Xeig);
    }

    // -------------------------------------------------------------------------
    // Fréchet helpers (block exponential)
    // -------------------------------------------------------------------------

    static void frechetExpBlock(final DenseMatrix64F A,
                                final DenseMatrix64F E,
                                final DenseMatrix64F out) {
        final int d = A.numRows;
        if (A.numCols != d || E.numRows != d || E.numCols != d ||
                out.numRows != d || out.numCols != d) {
            throw new IllegalArgumentException("Dimension mismatch in frechetExpBlock");
        }

        DenseMatrix64F B = new DenseMatrix64F(2 * d, 2 * d);
        CommonOps.insert(A, B, 0, 0);   // TL
        CommonOps.insert(E, B, 0, d);   // TR
        CommonOps.insert(A, B, d, d);   // BR

        double[][] MM = new double[2 * d][2 * d];
        for (int i = 0; i < 2 * d; i++) {
            for (int j = 0; j < 2 * d; j++) {
                MM[i][j] = B.get(i, j);
            }
        }
        DoubleMatrix2D Mc = new DenseDoubleMatrix2D(MM);
        DoubleMatrix2D EB = MatrixExponential.expmPade13(Mc);
        DoubleMatrix2D Kc = EB.viewPart(0, d, d, d).copy();

        for (int i = 0; i < d; i++) {
            for (int j = 0; j < d; j++) {
                out.set(i, j, Kc.getQuick(i, j));
            }
        }
    }

    static void frechetAdjExpBlock(final DenseMatrix64F A,
                                   final DenseMatrix64F X,
                                   final DenseMatrix64F out) {
        final int d = A.numRows;
        if (A.numCols != d || X.numRows != d || X.numCols != d ||
                out.numRows != d || out.numCols != d) {
            throw new IllegalArgumentException("Dimension mismatch in frechetAdjExpBlock");
        }

        DenseMatrix64F At = A.copy();
        CommonOps.transpose(At);

        DenseMatrix64F Xt = X.copy();
        CommonOps.transpose(Xt);

        DenseMatrix64F tmp = new DenseMatrix64F(d, d);
        frechetExpBlock(At, Xt, tmp);

        CommonOps.transpose(tmp, out);
    }

    // -------------------------------------------------------------------------
    // Matrix exponential via Colt (Padé(13))
    // -------------------------------------------------------------------------

    static final class MatrixExponential {

        private static final Algebra ALG = Algebra.DEFAULT;
        private static final double THETA_13 = 4.25;

        private MatrixExponential() {}

        public static DoubleMatrix2D expmPade13(DoubleMatrix2D A) {
            final int n = A.rows();
            if (n != A.columns()) {
                throw new IllegalArgumentException("A must be square");
            }
            if (n == 0) return new DenseDoubleMatrix2D(0, 0);
            if (isZero(A)) return eye(n);

            double A1 = norm1(A);
            int s = 0;
            if (A1 > THETA_13) {
                s = (int) Math.max(0, Math.ceil(log2(A1 / THETA_13)));
            }

            DoubleMatrix2D As = A.copy();
            if (s > 0) scaleInPlace(As, 1.0 / Math.pow(2.0, s));

            DoubleMatrix2D A2 = ALG.mult(As, As);
            DoubleMatrix2D A4 = ALG.mult(A2, A2);
            DoubleMatrix2D A6 = ALG.mult(A2, A4);

            final double[] c = new double[] {
                    64764752532480000.0,
                    32382376266240000.0,
                    7771770303897600.0,
                    1187353796428800.0,
                    129060195264000.0,
                    10559470521600.0,
                    670442572800.0,
                    33522128640.0,
                    1323241920.0,
                    40840800.0,
                    960960.0,
                    16380.0,
                    182.0,
                    1.0
            };

            DoubleMatrix2D I = eye(n);

            // Odd chain
            DoubleMatrix2D tmp1 = linComb(A6, c[13], A4, c[11]);
            axpyInPlace(tmp1, A2, c[9]);

            DoubleMatrix2D polyU = ALG.mult(A6, tmp1);
            axpyInPlace(polyU, A6, c[7]);
            axpyInPlace(polyU, A4, c[5]);
            axpyInPlace(polyU, A2, c[3]);
            axpyInPlace(polyU, I,  c[1]);
            DoubleMatrix2D U = ALG.mult(As, polyU);

            // Even chain
            DoubleMatrix2D innerV = linComb(A6, c[12], A4, c[10]);
            axpyInPlace(innerV, A2, c[8]);

            DoubleMatrix2D V = ALG.mult(A6, innerV);
            axpyInPlace(V, A6, c[6]);
            axpyInPlace(V, A4, c[4]);
            axpyInPlace(V, A2, c[2]);
            axpyInPlace(V, I,  c[0]);

            // Solve (V - U) X = (V + U)
            DoubleMatrix2D VmU = V.copy();
            axpyInPlace(VmU, U, -1.0);
            DoubleMatrix2D VpU = V.copy();
            axpyInPlace(VpU, U,  1.0);
            DoubleMatrix2D X = ALG.solve(VmU, VpU);

            for (int k = 0; k < s; k++) {
                X = ALG.mult(X, X);
            }
            return X;
        }

        private static DoubleMatrix2D eye(int n) {
            DenseDoubleMatrix2D I = new DenseDoubleMatrix2D(n, n);
            for (int i = 0; i < n; i++) I.setQuick(i, i, 1.0);
            return I;
        }

        private static boolean isZero(DoubleMatrix2D A) {
            final int r = A.rows(), c = A.columns();
            for (int i = 0; i < r; i++) {
                for (int j = 0; j < c; j++) {
                    if (A.getQuick(i, j) != 0.0) return false;
                }
            }
            return true;
        }

        private static double norm1(DoubleMatrix2D A) {
            final int n = A.rows();
            double max = 0.0;
            for (int j = 0; j < n; j++) {
                double s = 0.0;
                for (int i = 0; i < n; i++) {
                    s += Math.abs(A.getQuick(i, j));
                }
                if (s > max) max = s;
            }
            return max;
        }

        private static double log2(double x) {
            return Math.log(x) / Math.log(2.0);
        }

        private static void scaleInPlace(DoubleMatrix2D A, double alpha) {
            final int r = A.rows(), c = A.columns();
            for (int i = 0; i < r; i++) {
                for (int j = 0; j < c; j++) {
                    A.setQuick(i, j, A.getQuick(i, j) * alpha);
                }
            }
        }

        private static DoubleMatrix2D linComb(DoubleMatrix2D X, double a,
                                              DoubleMatrix2D Y, double b) {
            DenseDoubleMatrix2D Z = new DenseDoubleMatrix2D(X.rows(), X.columns());
            for (int i = 0; i < X.rows(); i++) {
                for (int j = 0; j < X.columns(); j++) {
                    Z.setQuick(i, j, a * X.getQuick(i, j) + b * Y.getQuick(i, j));
                }
            }
            return Z;
        }

        private static void axpyInPlace(DoubleMatrix2D X, DoubleMatrix2D Y, double alpha) {
            final int r = X.rows(), c = X.columns();
            for (int i = 0; i < r; i++) {
                for (int j = 0; j < c; j++) {
                    X.setQuick(i, j,
                            X.getQuick(i, j) + alpha * Y.getQuick(i, j));
                }
            }
        }
    }
}

package dr.evomodel.treedatalikelihood.continuous.backprop;

import dr.evolution.tree.NodeRef;
import dr.inference.model.BlockDiagonalCosSinMatrixParameter;
import dr.inference.model.Parameter;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

/**
 * Maps primitive OU gradients (dL/dS, ...) to gradients of:
 *   (a, r, theta)   and   R  in  S = R D R^{-1}.
 *
 * Output layout in target[] (starting at offset):
 *   [0 .. paramDim-1]                 : gradients wrt (a,r,theta)
 *   [paramDim .. paramDim + dim*dim-1]: gradients wrt R entries, row-major
 *
 * Notes:
 * - This class assumes dL/dS is already computed elsewhere (e.g., in the backprop strategy).
 * - No matrix inversion is performed here: R and R^{-1} are provided by the parameter.
 * - Uses D only in compressed form (diag | upper | lower), so dL/dR is computed without forming dense D.
 */
public final class BlockDiagCosSinPrimitiveGradientMapper
        implements SingleParameterOUPrimitiveGradientMapper {

    private final BlockDiagonalCosSinMatrixParameter parameter;

    private final int dim;
    private final int numBlocks;
    private final int num2x2Blocks;
    private final int paramDim;   // numBlocks + 2*num2x2Blocks
    private final int rDim;       // dim*dim
    private final int totalDim;   // paramDim + rDim

    // Cached block layout once (avoid clones per call)
    private final int[] blockStarts;
    private final int[] blockSizes;

    // Reused EJML matrices
    private final DenseMatrix64F R;
    private final DenseMatrix64F Rinv;

    // Workspaces (all dim×dim)
    private final DenseMatrix64F W1;
    private final DenseMatrix64F W2;
    private final DenseMatrix64F W3;
    private final DenseMatrix64F dLdD;   // R^{-T} dLdS R^T
    private final DenseMatrix64F dLdR;   // gradient wrt R

    // Reused buffers (no allocations per call)
    // source layout expected by parameter.chainGradient:
    // [diag (dim), upper (num2x2Blocks), lower (num2x2Blocks)]
    private final double[] source;
    private final double[] localGrad;     // length = paramDim
    private final double[] blockDParams;  // length = 3*dim - 2 (diag | upper | lower)

    public BlockDiagCosSinPrimitiveGradientMapper(BlockDiagonalCosSinMatrixParameter parameter) {
        this.parameter = parameter;

        this.dim = parameter.getRowDimension();
        this.numBlocks = parameter.getNumBlocks();
        this.num2x2Blocks = parameter.getNum2x2Blocks();

        this.paramDim = numBlocks + 2 * num2x2Blocks;
        this.rDim = dim * dim;
        this.totalDim = paramDim + rDim;

        // Cache block layout ONCE (parameter returns clones)
        this.blockStarts = parameter.getBlockStarts();
        this.blockSizes  = parameter.getBlockSizes();

        this.R    = new DenseMatrix64F(dim, dim);
        this.Rinv = new DenseMatrix64F(dim, dim);

        this.W1 = new DenseMatrix64F(dim, dim);
        this.W2 = new DenseMatrix64F(dim, dim);
        this.W3 = new DenseMatrix64F(dim, dim);

        this.dLdD = new DenseMatrix64F(dim, dim);
        this.dLdR = new DenseMatrix64F(dim, dim);

        this.source       = new double[dim + 2 * num2x2Blocks];
        this.localGrad    = new double[paramDim];
        this.blockDParams = new double[3 * dim - 2];
    }

    @Override
    public Parameter getParameter() {
        return parameter;
    }

    @Override
    public int getDimension() {
        return totalDim;
    }

    @Override
    public void mapPrimitiveToParameter(NodeRef node,
                                        DenseMatrix64F dLdS,
                                        DenseMatrix64F dLdSigmaSt,
                                        DenseMatrix64F dLdMu,
                                        DenseMatrix64F dLdSigma,
                                        double[] target,
                                        int offset) {

        // --- Fetch R and R^{-1} from the parameter (NO inversion here) ---
        parameter.fillRAndRinv(R.getData(), Rinv.getData()); //TODO THIS IS A USELESS ALLOCATION, FIX LATER

        // --- Fetch D in compressed form (diag|upper|lower), NO allocations ---
        parameter.fillBlockDiagonalElements(blockDParams);

        // -----------------------------------------------------------------
        // 1) dL/dD = R^{-T} (dL/dS) R^T
        //    dLdD = (Rinv)^T dLdS R^T
        // -----------------------------------------------------------------
        CommonOps.multTransA(Rinv, dLdS, W1);     // W1 = (Rinv)^T dLdS
        CommonOps.multTransB(W1, R, dLdD);        // dLdD = W1 R^T

        // -----------------------------------------------------------------
        // 2) Pack dLdD into source[] expected by chainGradient(...)
        //    source = [diag | upper(2x2 blocks) | lower(2x2 blocks)]
        // -----------------------------------------------------------------
        packDLdDToSource(dLdD, source);

        // -----------------------------------------------------------------
        // 3) chainGradient -> gradients w.r.t (a, r, theta)
        // -----------------------------------------------------------------
        parameter.chainGradient(source, localGrad);
        for (int i = 0; i < paramDim; i++) {
            target[offset + i] = localGrad[i];
        }

        // -----------------------------------------------------------------
        // 4) Compute gradient w.r.t. R (dense dim×dim), then flatten row-major
        //
        // Derivation (Frobenius inner product):
        //   S = R D R^{-1}
        //   dLdR = (dLdS) R^{-T} D^T  -  S^T (dLdS) R^{-T}
        //        = (dLdS) R^{-T} D^T  -  R^{-T} [ D^T ( R^T (dLdS) R^{-T} ) ]
        //
        // We evaluate both terms using only:
        //   - dense multiplies with R, Rinv
        //   - tri-diagonal multiplies with D^T in compressed form
        // -----------------------------------------------------------------
        computeDLdR(dLdS, blockDParams, R, Rinv, dLdR, W1, W2, W3);

        // Flatten dLdR row-major into the trailing block
        final double[] rData = dLdR.data;
        final int baseR = offset + paramDim;
        for (int i = 0; i < rDim; i++) {
            target[baseR + i] = rData[i];
        }
    }

    // ---------------------------------------------------------------------
    // Packing: dLdD (dense) -> source (diag + per-2x2 off-diags)
    // ---------------------------------------------------------------------

    private void packDLdDToSource(DenseMatrix64F dLdD, double[] sourceOut) {
        final double[] data = dLdD.data;

        // diag
        for (int i = 0; i < dim; i++) {
            sourceOut[i] = data[i * dim + i];
        }

        // upper/lower, only at 2x2 block starts, in the SAME order as parameter.chainGradient expects
        final int baseUpper = dim;
        final int baseLower = dim + num2x2Blocks;

        int idx2 = 0;
        for (int b = 0; b < numBlocks; b++) {
            int start = blockStarts[b];
            if (blockSizes[b] == 2) {
                sourceOut[baseUpper + idx2] = data[start * dim + (start + 1)];
                sourceOut[baseLower + idx2] = data[(start + 1) * dim + start];
                idx2++;
            }
        }
        if (idx2 != num2x2Blocks) {
            throw new IllegalStateException("Expected " + num2x2Blocks + " 2x2 blocks, saw " + idx2);
        }
    }

    // ---------------------------------------------------------------------
    // dL/dR for S = R D R^{-1}, using (dLdS, R, Rinv, D compressed)
    //
    // out_dLdR = (dLdS) R^{-T} D^T  -  R^{-T} [ D^T ( R^T (dLdS) R^{-T} ) ].
    // ---------------------------------------------------------------------

    private static void computeDLdR(DenseMatrix64F dLdS,
                                    double[] blockDParams,
                                    DenseMatrix64F R,
                                    DenseMatrix64F Rinv,
                                    DenseMatrix64F out_dLdR,
                                    DenseMatrix64F W1,
                                    DenseMatrix64F W2,
                                    DenseMatrix64F W3) {

        final int n = R.getNumRows();

        // Term 1: (dLdS) R^{-T} D^T
        // W1 = (dLdS) R^{-T} = dLdS (Rinv)^T
        CommonOps.multTransB(dLdS, Rinv, W1);
        // W2 = W1 D^T
        rightMultiplyByDtCompressed(W1, blockDParams, W2, n);

        // Term 2: R^{-T} [ D^T ( R^T (dLdS) R^{-T} ) ]
        // W1 = R^T dLdS
        CommonOps.multTransA(R, dLdS, W1);
        // W3 = R^T dLdS R^{-T}
        CommonOps.multTransB(W1, Rinv, W3);
        // W1 = D^T * W3
        leftMultiplyByDtCompressed(blockDParams, W3, W1, n);
        // W3 = R^{-T} * W1
        CommonOps.multTransA(Rinv, W1, W3);

        // out = term1 - term2
        CommonOps.subtract(W2, W3, out_dLdR);
    }

    // out = M * D^T, with D compressed as [diag | upper | lower]
    // D[i,i]   = d_i
    // D[i,i+1] = u_i
    // D[i+1,i] = l_i
    private static void rightMultiplyByDtCompressed(DenseMatrix64F M,
                                                    double[] blockDParams,
                                                    DenseMatrix64F out,
                                                    int n) {

        final int upperOffset = n;
        final int lowerOffset = n + (n - 1);

        final double[] m = M.data;
        final double[] o = out.data;

        for (int r = 0; r < n; r++) {
            final int rowBase = r * n;

            for (int c = 0; c < n; c++) {
                // D^T has:
                // (c,c)     = d_c
                // (c-1,c)   = D[c,c-1]   = l_{c-1}
                // (c+1,c)   = D[c,c+1]   = u_c
                double val = m[rowBase + c] * blockDParams[c]; // d_c

                if (c > 0) {
                    val += m[rowBase + (c - 1)] * blockDParams[lowerOffset + (c - 1)]; // l_{c-1}
                }
                if (c < n - 1) {
                    val += m[rowBase + (c + 1)] * blockDParams[upperOffset + c];       // u_c
                }

                o[rowBase + c] = val;
            }
        }
    }

    // out = D^T * M, with D compressed as [diag | upper | lower]
    private static void leftMultiplyByDtCompressed(double[] blockDParams,
                                                   DenseMatrix64F M,
                                                   DenseMatrix64F out,
                                                   int n) {

        final int upperOffset = n;
        final int lowerOffset = n + (n - 1);

        final double[] m = M.data;
        final double[] o = out.data;

        for (int r = 0; r < n; r++) {
            final double dr = blockDParams[r];
            final int rowBase = r * n;

            // diag contribution
            for (int c = 0; c < n; c++) {
                o[rowBase + c] = dr * m[rowBase + c];
            }

            // from (r-1): D^T[r,r-1] = D[r-1,r] = u_{r-1}
            if (r > 0) {
                final double urm1 = blockDParams[upperOffset + (r - 1)];
                final int rowPrev = (r - 1) * n;
                for (int c = 0; c < n; c++) {
                    o[rowBase + c] += urm1 * m[rowPrev + c];
                }
            }

            // from (r+1): D^T[r,r+1] = D[r+1,r] = l_r
            if (r < n - 1) {
                final double lr = blockDParams[lowerOffset + r];
                final int rowNext = (r + 1) * n;
                for (int c = 0; c < n; c++) {
                    o[rowBase + c] += lr * m[rowNext + c];
                }
            }
        }
    }
}
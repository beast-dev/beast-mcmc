package dr.evomodel.treedatalikelihood.continuous.backprop;

import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.util.Arrays;

/**
 * Exact block-space helper for Fréchet derivatives of exp(-D t).
 *
 * <p>The forward block integral is represented as a cached exact linear map on each
 * 1x1 / 2x2 block pair. This avoids numerical quadrature on the hot path and turns
 * repeated branch-local applications into tiny dense matrix-vector multiplies.</p>
 */
public final class BlockDiagonalFrechetHelper {

    private static final String QUADRATURE_DEBUG_PROPERTY =
            "beast.experimental.nativeUseQuadratureFrechet";

    private final int dim;
    private final BlockDiagonalExpSolver.BlockStructure structure;
    private final BlockDiagonalExpSolver expSolver;

    private final double[] blockDTransposeParams;
    private final double[] buf2a = new double[2];
    private final double[] buf4a = new double[4];
    private final double[] buf4b = new double[4];
    private final double[] buf4c = new double[4];
    private final DenseMatrix64F lastExp;
    private final BlockDiagonalFrechetExactPlan.Workspace workspace;

    private final double[] cachedPlanParams;
    private BlockDiagonalFrechetExactPlan.Plan cachedPlan;
    private double cachedPlanT;
    private boolean cachedPlanValid;

    public BlockDiagonalFrechetHelper(final BlockDiagonalExpSolver.BlockStructure structure) {
        this.structure = structure;
        this.dim = structure.getDim();
        this.expSolver = new BlockDiagonalExpSolver(structure);
        this.blockDTransposeParams = new double[3 * dim - 2];
        this.lastExp = new DenseMatrix64F(dim, dim);
        this.workspace = new BlockDiagonalFrechetExactPlan.Workspace();
        this.cachedPlanParams = new double[3 * dim - 2];
        this.cachedPlan = null;
        this.cachedPlanT = Double.NaN;
        this.cachedPlanValid = false;
    }

    public DenseMatrix64F computeExpInDBasis(final double[] blockDParams,
                                             final double t) {
        expSolver.compute(blockDParams, t, lastExp);
        return lastExp;
    }

    public void frechetAdjointExpInDBasis(final double[] blockDParams,
                                          final DenseMatrix64F xD,
                                          final double t,
                                          final DenseMatrix64F out) {
        validateInputs(blockDParams, xD, out);
        buildTransposeParams(blockDParams, blockDTransposeParams);
        computeForwardFrechetInDBasis(blockDTransposeParams, xD, t, out);
        CommonOps.transpose(out);
        CommonOps.scale(-t, out);
    }

    public void computeForwardFrechetInDBasis(final double[] blockDParams,
                                              final DenseMatrix64F eD,
                                              final double t,
                                              final DenseMatrix64F out) {
        validateInputs(blockDParams, eD, out);
        Arrays.fill(out.data, 0.0);
        if (Boolean.getBoolean(QUADRATURE_DEBUG_PROPERTY)) {
            computeForwardFrechetInDBasisQuadrature(blockDParams, eD, t, out);
            return;
        }
        getOrBuildPlan(blockDParams, t).apply(eD, out, workspace);
    }

    private void buildTransposeParams(final double[] src, final double[] dst) {
        final int diagLen = dim;
        final int offLen = dim - 1;
        System.arraycopy(src, 0, dst, 0, diagLen);
        System.arraycopy(src, diagLen + offLen, dst, diagLen, offLen);
        System.arraycopy(src, diagLen, dst, diagLen + offLen, offLen);
    }

    private void validateInputs(final double[] blockDParams,
                                final DenseMatrix64F in,
                                final DenseMatrix64F out) {
        BlockDiagonalExpSolver.validateCompressedParams(blockDParams, dim);
        validateMatrix(in, "input");
        validateMatrix(out, "out");
    }

    private void validateMatrix(final DenseMatrix64F x, final String name) {
        if (x.numRows != dim || x.numCols != dim) {
            throw new IllegalArgumentException(
                    name + " must be " + dim + "x" + dim + " but is " + x.numRows + "x" + x.numCols);
        }
    }

    private void computeForwardFrechetInDBasisQuadrature(final double[] blockDParams,
                                                         final DenseMatrix64F eD,
                                                         final double t,
                                                         final DenseMatrix64F out) {
        for (int bi = 0; bi < structure.getNumBlocks(); bi++) {
            final int startI = structure.getBlockStart(bi);
            final int sizeI = structure.getBlockSize(bi);
            for (int bj = 0; bj < structure.getNumBlocks(); bj++) {
                final int startJ = structure.getBlockStart(bj);
                final int sizeJ = structure.getBlockSize(bj);
                computeBlockFrechetIntegralQuadrature(blockDParams, eD, t, startI, sizeI, startJ, sizeJ, out);
            }
        }
    }

    private void computeBlockFrechetIntegralQuadrature(final double[] blockDParams,
                                                       final DenseMatrix64F eTilde,
                                                       final double t,
                                                       final int startI,
                                                       final int sizeI,
                                                       final int startJ,
                                                       final int sizeJ,
                                                       final DenseMatrix64F out) {
        final int upperOffset = dim;
        final int lowerOffset = dim + (dim - 1);

        if (sizeI == 1 && sizeJ == 1) {
            final double lambdaI = -t * blockDParams[startI];
            final double lambdaJ = -t * blockDParams[startJ];
            final double eIJ = eTilde.get(startI, startJ);
            final double value;
            if (Math.abs(lambdaI - lambdaJ) < 1.0e-10) {
                value = eIJ * Math.exp(lambdaI);
            } else {
                value = eIJ * (Math.exp(lambdaI) - Math.exp(lambdaJ)) / (lambdaI - lambdaJ);
            }
            out.set(startI, startJ, value);
            return;
        }

        if (sizeI == 1 && sizeJ == 2) {
            final double lambdaI = -t * blockDParams[startI];
            final double a = -t * blockDParams[startJ];
            final double d = -t * blockDParams[startJ + 1];
            final double b = -t * blockDParams[upperOffset + startJ];
            final double c = -t * blockDParams[lowerOffset + startJ];
            final double e1 = eTilde.get(startI, startJ);
            final double e2 = eTilde.get(startI, startJ + 1);

            BlockDiagonalFrechetIntegrator.integrate1x1_2x2(
                    lambdaI, a, b, c, d, e1, e2, buf2a, buf4a);
            out.set(startI, startJ, buf2a[0]);
            out.set(startI, startJ + 1, buf2a[1]);
            return;
        }

        if (sizeI == 2 && sizeJ == 1) {
            final double a = -t * blockDParams[startI];
            final double d = -t * blockDParams[startI + 1];
            final double b = -t * blockDParams[upperOffset + startI];
            final double c = -t * blockDParams[lowerOffset + startI];
            final double lambdaJ = -t * blockDParams[startJ];
            final double e1 = eTilde.get(startI, startJ);
            final double e2 = eTilde.get(startI + 1, startJ);

            BlockDiagonalFrechetIntegrator.integrate2x2_1x1(
                    a, b, c, d, lambdaJ, e1, e2, buf2a, buf4a);
            out.set(startI, startJ, buf2a[0]);
            out.set(startI + 1, startJ, buf2a[1]);
            return;
        }

        final double ai = -t * blockDParams[startI];
        final double di = -t * blockDParams[startI + 1];
        final double bi = -t * blockDParams[upperOffset + startI];
        final double ci = -t * blockDParams[lowerOffset + startI];

        final double aj = -t * blockDParams[startJ];
        final double dj = -t * blockDParams[startJ + 1];
        final double bj = -t * blockDParams[upperOffset + startJ];
        final double cj = -t * blockDParams[lowerOffset + startJ];

        final double e11 = eTilde.get(startI, startJ);
        final double e12 = eTilde.get(startI, startJ + 1);
        final double e21 = eTilde.get(startI + 1, startJ);
        final double e22 = eTilde.get(startI + 1, startJ + 1);

        BlockDiagonalFrechetIntegrator.integrate2x2_2x2(
                ai, bi, ci, di,
                aj, bj, cj, dj,
                e11, e12, e21, e22,
                buf4a, buf4b, buf4c
        );
        out.set(startI, startJ, buf4a[0]);
        out.set(startI, startJ + 1, buf4a[1]);
        out.set(startI + 1, startJ, buf4a[2]);
        out.set(startI + 1, startJ + 1, buf4a[3]);
    }

    private BlockDiagonalFrechetExactPlan.Plan getOrBuildPlan(final double[] blockDParams,
                                                              final double t) {
        if (!cachedPlanValid
                || Double.doubleToLongBits(t) != Double.doubleToLongBits(cachedPlanT)
                || !Arrays.equals(blockDParams, cachedPlanParams)) {
            cachedPlan = BlockDiagonalFrechetExactPlan.buildPlan(structure, blockDParams, t);
            System.arraycopy(blockDParams, 0, cachedPlanParams, 0, blockDParams.length);
            cachedPlanT = t;
            cachedPlanValid = true;
        }
        return cachedPlan;
    }
}

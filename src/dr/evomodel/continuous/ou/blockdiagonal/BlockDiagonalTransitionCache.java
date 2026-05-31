package dr.evomodel.continuous.ou.blockdiagonal;

import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalExpSolver;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.BlockDiagonalDecomposition;
import org.ejml.data.DenseMatrix64F;

import java.util.Arrays;

/**
 * Reusable branch-length dependent transition cache for block-diagonal OU paths.
 *
 * The decomposition owns R, R^{-1}, and D.  This cache owns only dt-dependent
 * derived values and scratch space.
 */
final class BlockDiagonalTransitionCache {

    private static final String DEBUG_NATIVE_R_CONSISTENCY_PROPERTY =
            "beast.debug.nativeRConsistency";

    final BlockDiagonalDecomposition decomposition;
    final double[] rData;
    final double[] rinvData;
    final double[] blockDParams;
    final int[] blockStarts;
    final int[] blockSizes;
    final double[] expD;
    final DenseMatrix64F rMatrix;
    final DenseMatrix64F rinvMatrix;
    final DenseMatrix64F transitionMatrix;

    private final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter;
    private final BlockDiagonalExpSolver expSolver;
    private final double[] cachedRData;
    private final double[] cachedRinvData;
    private final double[] cachedBlockDParams;
    private final double[] workMatrix;
    private final double[] debugBlockR;
    private final double[] debugBlockRinv;
    private boolean expCacheValid;
    private boolean basisCacheValid;
    private double cachedExpDt;
    private double cachedBasisDt;
    private int cachedRHash;
    private int cachedRinvHash;
    private int cachedBlockDHash;
    private long expParameterVersion;
    private long basisExpVersion;

    BlockDiagonalTransitionCache(final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
                                 final BlockDiagonalExpSolver expSolver) {
        this.blockParameter = blockParameter;
        this.expSolver = expSolver;
        final int d = blockParameter.getRowDimension();
        this.decomposition = blockParameter.createBlockDiagonalDecomposition();
        this.rData = decomposition.getR();
        this.rinvData = decomposition.getRInverse();
        this.blockDParams = decomposition.getBlockDiagonal();
        this.blockStarts = decomposition.getBlockStarts();
        this.blockSizes = decomposition.getBlockSizes();
        this.cachedRData = new double[d * d];
        this.cachedRinvData = new double[d * d];
        this.cachedBlockDParams = new double[blockParameter.getTridiagonalDDimension()];
        this.expD = new double[blockParameter.getCompressedDDimension()];
        this.rMatrix = DenseMatrix64F.wrap(d, d, rData);
        this.rinvMatrix = DenseMatrix64F.wrap(d, d, rinvData);
        this.transitionMatrix = new DenseMatrix64F(d, d);
        this.workMatrix = new double[d * d];
        this.debugBlockR = new double[d * d];
        this.debugBlockRinv = new double[d * d];
        this.expCacheValid = false;
        this.basisCacheValid = false;
        this.cachedExpDt = Double.NaN;
        this.cachedBasisDt = Double.NaN;
        this.cachedRHash = 0;
        this.cachedRinvHash = 0;
        this.cachedBlockDHash = 0;
        this.expParameterVersion = 0L;
        this.basisExpVersion = -1L;
    }

    void refresh(final double dt) {
        final int d = blockParameter.getRowDimension();
        blockParameter.fillBlockDiagonalDecomposition(decomposition);

        final boolean expNeedsRefresh = expNeedsRefresh(dt);
        if (Boolean.getBoolean(DEBUG_NATIVE_R_CONSISTENCY_PROPERTY)) {
            reportNativeRConsistency(dt, d);
        }
        final int rHash = Arrays.hashCode(rData);
        final int rinvHash = Arrays.hashCode(rinvData);
        final boolean rChanged = rHash != cachedRHash || !Arrays.equals(rData, cachedRData);
        final boolean rinvChanged = rinvHash != cachedRinvHash || !Arrays.equals(rinvData, cachedRinvData);
        final boolean basisNeedsRefresh = expNeedsRefresh
                || !basisCacheValid
                || Double.doubleToLongBits(dt) != Double.doubleToLongBits(cachedBasisDt)
                || basisExpVersion != expParameterVersion
                || rChanged
                || rinvChanged;
        if (!basisNeedsRefresh) {
            return;
        }

        if (expNeedsRefresh) {
            refreshExp(dt);
        }
        BlockDiagonalMatrixOps.fillTransitionMatrix(
                rData,
                expD,
                rinvData,
                d,
                blockStarts,
                blockSizes,
                workMatrix,
                transitionMatrix.data);
        System.arraycopy(rData, 0, cachedRData, 0, rData.length);
        System.arraycopy(rinvData, 0, cachedRinvData, 0, rinvData.length);
        cachedRHash = rHash;
        cachedRinvHash = rinvHash;
        cachedBasisDt = dt;
        basisExpVersion = expParameterVersion;
        basisCacheValid = true;
    }

    void refreshExpOnly(final double dt) {
        blockParameter.fillBlockDiagonalElements(blockDParams);
        if (expNeedsRefresh(dt)) {
            refreshExp(dt);
            basisCacheValid = false;
        }
    }

    private boolean expNeedsRefresh(final double dt) {
        final int blockDHash = Arrays.hashCode(blockDParams);
        return !expCacheValid
                || Double.doubleToLongBits(dt) != Double.doubleToLongBits(cachedExpDt)
                || blockDHash != cachedBlockDHash
                || !Arrays.equals(blockDParams, cachedBlockDParams);
    }

    private void refreshExp(final double dt) {
        expSolver.computeCompressed(blockDParams, dt, expD);
        System.arraycopy(blockDParams, 0, cachedBlockDParams, 0, blockDParams.length);
        cachedBlockDHash = Arrays.hashCode(blockDParams);
        expParameterVersion++;
        cachedExpDt = dt;
        expCacheValid = true;
    }

    private void reportNativeRConsistency(final double dt, final int d) {
        blockParameter.fillRAndRinv(debugBlockR, debugBlockRinv);
        double maxR = 0.0;
        double maxRinv = 0.0;
        double maxRinvVsTranspose = 0.0;
        for (int i = 0; i < d * d; ++i) {
            maxR = Math.max(maxR, Math.abs(rData[i] - debugBlockR[i]));
            maxRinv = Math.max(maxRinv, Math.abs(rinvData[i] - debugBlockRinv[i]));
        }
        for (int r = 0; r < d; ++r) {
            for (int c = 0; c < d; ++c) {
                final int rc = r * d + c;
                final int cr = c * d + r;
                maxRinvVsTranspose = Math.max(
                        maxRinvVsTranspose,
                        Math.abs(debugBlockRinv[rc] - debugBlockR[cr]));
            }
        }
        dr.evomodel.treedatalikelihood.continuous.canonical.CanonicalDiagnosticsLog.warning("nativeRConsistencyDebug dt=" + dt
                + " maxAbsRDiff=" + maxR
                + " maxAbsRinvDiff=" + maxRinv
                + " maxAbsRinvVsTranspose=" + maxRinvVsTranspose);
    }
}

package dr.evomodel.continuous.ou.orthogonalblockdiagonal;

import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalExpSolver;
import dr.evomodel.treedatalikelihood.continuous.canonical.math.MatrixOps;
import dr.inference.model.AbstractBlockDiagonalTwoByTwoMatrixParameter;
import dr.inference.model.OrthogonalMatrixProvider;
import org.ejml.data.DenseMatrix64F;

import java.util.Arrays;

/**
 * Reusable branch-length dependent basis cache for orthogonal block OU paths.
 */
final class OrthogonalBlockBasisCache {

    private static final String DEBUG_NATIVE_R_CONSISTENCY_PROPERTY =
            "beast.debug.nativeRConsistency";

    final double[] rData;
    final double[] rtData;
    final double[] blockDParams;
    final int[] blockStarts;
    final int[] blockSizes;
    final double[] expD;
    final DenseMatrix64F rMatrix;
    final DenseMatrix64F rtMatrix;
    final DenseMatrix64F transitionMatrix;

    private final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter;
    private final OrthogonalMatrixProvider orthogonalRotation;
    private final BlockDiagonalExpSolver expSolver;
    private final double[] cachedRData;
    private final double[] cachedRtData;
    private final double[] cachedBlockDParams;
    private final double[] workMatrix;
    private final double[] debugBlockR;
    private final double[] debugBlockRinv;
    private boolean expCacheValid;
    private boolean basisCacheValid;
    private double cachedExpDt;
    private double cachedBasisDt;
    private int cachedRHash;
    private int cachedRtHash;
    private int cachedBlockDHash;
    private long expParameterVersion;
    private long basisRVersion;
    private long basisRtVersion;
    private long basisExpVersion;

    OrthogonalBlockBasisCache(final AbstractBlockDiagonalTwoByTwoMatrixParameter blockParameter,
                              final OrthogonalMatrixProvider orthogonalRotation,
                              final BlockDiagonalExpSolver expSolver) {
        this.blockParameter = blockParameter;
        this.orthogonalRotation = orthogonalRotation;
        this.expSolver = expSolver;
        final int d = blockParameter.getRowDimension();
        this.rData = new double[d * d];
        this.rtData = new double[d * d];
        this.blockDParams = new double[blockParameter.getTridiagonalDDimension()];
        this.blockStarts = blockParameter.getBlockStarts();
        this.blockSizes = blockParameter.getBlockSizes();
        this.cachedRData = new double[d * d];
        this.cachedRtData = new double[d * d];
        this.cachedBlockDParams = new double[blockParameter.getTridiagonalDDimension()];
        this.expD = new double[blockParameter.getCompressedDDimension()];
        this.rMatrix = new DenseMatrix64F(d, d);
        this.rtMatrix = new DenseMatrix64F(d, d);
        this.transitionMatrix = new DenseMatrix64F(d, d);
        this.workMatrix = new double[d * d];
        this.debugBlockR = new double[d * d];
        this.debugBlockRinv = new double[d * d];
        this.expCacheValid = false;
        this.basisCacheValid = false;
        this.cachedExpDt = Double.NaN;
        this.cachedBasisDt = Double.NaN;
        this.cachedRHash = 0;
        this.cachedRtHash = 0;
        this.cachedBlockDHash = 0;
        this.expParameterVersion = 0L;
        this.basisRVersion = 0L;
        this.basisRtVersion = 0L;
        this.basisExpVersion = -1L;
    }

    void refresh(final double dt) {
        final int d = blockParameter.getRowDimension();
        orthogonalRotation.fillOrthogonalMatrix(rData);
        orthogonalRotation.fillOrthogonalTranspose(rtData);
        blockParameter.fillBlockDiagonalElements(blockDParams);

        final boolean expNeedsRefresh = expNeedsRefresh(dt);
        if (Boolean.getBoolean(DEBUG_NATIVE_R_CONSISTENCY_PROPERTY)) {
            reportNativeRConsistency(dt, d);
        }
        final int rHash = Arrays.hashCode(rData);
        final int rtHash = Arrays.hashCode(rtData);
        final boolean rChanged = rHash != cachedRHash || !Arrays.equals(rData, cachedRData);
        final boolean rtChanged = rtHash != cachedRtHash || !Arrays.equals(rtData, cachedRtData);
        final boolean basisNeedsRefresh = expNeedsRefresh
                || !basisCacheValid
                || Double.doubleToLongBits(dt) != Double.doubleToLongBits(cachedBasisDt)
                || basisExpVersion != expParameterVersion
                || rChanged
                || rtChanged;
        if (!basisNeedsRefresh) {
            return;
        }

        if (expNeedsRefresh) {
            refreshExp(dt);
        }
        if (rChanged) {
            basisRVersion++;
        }
        if (rtChanged) {
            basisRtVersion++;
        }

        MatrixOps.fromFlat(rData, rMatrix, d);
        MatrixOps.fromFlat(rtData, rtMatrix, d);
        OrthogonalBlockPreparedBasisBuilder.fillTransitionMatrix(
                rData,
                expD,
                rtData,
                d,
                blockStarts,
                blockSizes,
                workMatrix,
                transitionMatrix.data);
        System.arraycopy(rData, 0, cachedRData, 0, rData.length);
        System.arraycopy(rtData, 0, cachedRtData, 0, rtData.length);
        cachedRHash = rHash;
        cachedRtHash = rtHash;
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
        double maxRt = 0.0;
        double maxRinvVsTranspose = 0.0;
        for (int i = 0; i < d * d; ++i) {
            maxR = Math.max(maxR, Math.abs(rData[i] - debugBlockR[i]));
            maxRt = Math.max(maxRt, Math.abs(rtData[i] - debugBlockRinv[i]));
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
                + " maxAbsRtDiff=" + maxRt
                + " maxAbsRinvVsTranspose=" + maxRinvVsTranspose);
    }
}

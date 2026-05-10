package dr.evomodel.treedatalikelihood.discrete.beastBasedDiscreteTreeLikelihood.representations;

import dr.evomodel.substmodel.BaseSubstitutionModel;
import dr.evomodel.substmodel.EigenDecomposition;
import org.apache.commons.math.util.FastMath;

import java.util.Arrays;

/**
 * Spectral-basis partial representation with spectral branch propagation.
 *
 * Stores post-order partials as y = R^{-1} p and pre-order partials as z = R^T q.
 * Traversal stays in rotated coordinates; merge operations rotate only transient
 * working vectors.
 *
 * Branch coefficients are cached with CACHE_SLOTS slots per node so that all rate
 * categories coexist in the cache without evicting each other.  All cache storage
 * uses flat 1-D primitive arrays to avoid pointer-chasing through Object[][][] chains.
 *
 * Layout:
 *   coeffA / coeffB : double[], index = (node * CACHE_SLOTS + slot) * stateCount + state
 *   cachedBranchLength / cachedCoeffVersion : index = node * CACHE_SLOTS + slot
 */
public final class SpectralRotatedPartialsRepresentation
        implements BidirectionalRepresentation {

    private static final double IMAG_TOL = 1.0e-15;
    private static final double PAIR_TOL = 1.0e-12;

    private static final byte REAL_BLOCK    = 1;
    private static final byte COMPLEX_BLOCK = 2;

    // Number of branch-length slots cached per node.
    // 16 is more than enough for typical category counts (4–8).
    private static final int CACHE_SLOTS = 16;

    private final BaseSubstitutionModel substitutionModel;
    private final int stateCount;

    private final double[] matrixR;
    private final double[] matrixRInv;
    private final double[] eigenReal;
    private final double[] eigenImag;

    private final byte[] blockType;
    private final int[]  blockStart;
    private int blockCount;

    private final double[] tmpA;
    private final double[] tmpB;
    private final double[] tmpStandardA;
    private final double[] tmpStandardB;

    private boolean eigenDirty = true;
    private long    eigenVersion = 0L;

    // --- flat cache arrays ---
    // coeffA[node * CACHE_SLOTS * K + slot * K + state]  (K = stateCount)
    private double[] coeffA;
    private double[] coeffB;
    // cachedBranchLength[node * CACHE_SLOTS + slot]
    private double[] cachedBranchLength;
    // cachedCoeffVersion[node * CACHE_SLOTS + slot]
    private long[]   cachedCoeffVersion;
    // nextSlotByNode[node]  — round-robin eviction pointer
    private int[]    nextSlotByNode;
    // activeSlotByNode[node] — slot last found/written (for apply methods)
    private int[]    activeSlotByNode;
    private final int nodeCapacity;

    public SpectralRotatedPartialsRepresentation(BaseSubstitutionModel substitutionModel,
                                                  int nodeCount) {
        if (substitutionModel == null) {
            throw new IllegalArgumentException("substitutionModel must not be null");
        }

        this.substitutionModel = substitutionModel;
        this.stateCount = substitutionModel.getDataType().getStateCount();

        this.matrixR    = new double[stateCount * stateCount];
        this.matrixRInv = new double[stateCount * stateCount];
        this.eigenReal  = new double[stateCount];
        this.eigenImag  = new double[stateCount];

        this.blockType  = new byte[stateCount];
        this.blockStart = new int[stateCount];

        this.tmpA         = new double[stateCount];
        this.tmpB         = new double[stateCount];
        this.tmpStandardA = new double[stateCount];
        this.tmpStandardB = new double[stateCount];

        this.nodeCapacity       = nodeCount;
        this.coeffA             = new double[nodeCount * CACHE_SLOTS * stateCount];
        this.coeffB             = new double[nodeCount * CACHE_SLOTS * stateCount];
        this.cachedBranchLength = new double[nodeCount * CACHE_SLOTS];
        this.cachedCoeffVersion = new long[nodeCount * CACHE_SLOTS];
        this.nextSlotByNode     = new int[nodeCount];
        this.activeSlotByNode   = new int[nodeCount];

        Arrays.fill(cachedBranchLength, Double.NaN);
        Arrays.fill(cachedCoeffVersion, Long.MIN_VALUE);
    }

    @Override public String  getName()                { return "spectral-rotated-partials"; }
    @Override public int     getStateCount()          { return stateCount; }
    @Override public boolean supportsInternalScaling() { return true; }

    @Override
    public void markDirty() { eigenDirty = true; }

    @Override
    public void storeState() { /* caches are derived; nothing to snapshot */ }

    @Override
    public void restoreState() { eigenDirty = true; }

    @Override
    public void updateForLikelihood() { ensureEigenSystemCurrent(); }

    @Override
    public void initializeTipPartial(double[] standardTip, double[] outPartial) {
        ensureEigenSystemCurrent();
        multiplyMatrixVector(matrixRInv, standardTip, outPartial, stateCount);
    }

    @Override
    public void propagateToBranchTop(int nodeNumber,
                                     double branchLength,
                                     double[] childPartial,
                                     double[] outBranchTopPartial) {
        ensureEigenSystemCurrent();
        ensureBranchCoefficients(nodeNumber, branchLength);
        applyForwardCoefficients(nodeNumber, childPartial, outBranchTopPartial);
    }

    @Override
    public void combineBranchTopPartials(double[] leftBranchTopPartial,
                                         double[] rightBranchTopPartial,
                                         double[] outParentPartial) {
        ensureEigenSystemCurrent();
        multiplyMatrixVector(matrixR,    leftBranchTopPartial,  tmpStandardA, stateCount);
        multiplyMatrixVector(matrixR,    rightBranchTopPartial, tmpStandardB, stateCount);
        for (int i = 0; i < stateCount; i++) {
            tmpA[i] = tmpStandardA[i] * tmpStandardB[i];
        }
        multiplyMatrixVector(matrixRInv, tmpA, outParentPartial, stateCount);
    }

    @Override
    public double rootContribution(double[] rootFrequencies, double[] rootPartial) {
        ensureEigenSystemCurrent();
        multiplyTransposeMatrixVector(matrixR, rootFrequencies, tmpA, stateCount);
        double sum = 0.0;
        for (int i = 0; i < stateCount; i++) sum += tmpA[i] * rootPartial[i];
        return sum;
    }

    @Override
    public void initializeRootPartial(double[] rootFrequencies, double[] outPartial) {
        ensureEigenSystemCurrent();
        multiplyTransposeMatrixVector(matrixR, rootFrequencies, outPartial, stateCount);
    }

    @Override
    public void combineParentAndSibling(double[] parentNodePreOrder,
                                        double[] siblingBranchTopPostOrder,
                                        double[] outChildBranchTopPreOrder) {
        ensureEigenSystemCurrent();
        multiplyTransposeMatrixVector(matrixRInv, parentNodePreOrder,      tmpStandardA, stateCount);
        multiplyMatrixVector         (matrixR,    siblingBranchTopPostOrder, tmpStandardB, stateCount);
        for (int i = 0; i < stateCount; i++) {
            tmpStandardB[i] = tmpStandardA[i] * tmpStandardB[i];
        }
        multiplyTransposeMatrixVector(matrixR, tmpStandardB, outChildBranchTopPreOrder, stateCount);
    }

    @Override
    public void propagateToBranchBottom(int childNodeNumber,
                                        double branchLength,
                                        double[] childBranchTopPreOrder,
                                        double[] outChildNodePreOrder) {
        ensureEigenSystemCurrent();
        ensureBranchCoefficients(childNodeNumber, branchLength);
        applyTransposeCoefficients(childNodeNumber, childBranchTopPreOrder, outChildNodePreOrder);
    }

    @Override
    public void exportPostOrderPartial(double[] partial, double[] outPartial) {
        System.arraycopy(partial, 0, outPartial, 0, stateCount);
    }

    @Override
    public void exportPreOrderPartial(double[] preOrderPartial, double[] outPartial) {
        System.arraycopy(preOrderPartial, 0, outPartial, 0, stateCount);
    }

    // -------------------------------------------------------------------------

    private void ensureEigenSystemCurrent() {
        if (!eigenDirty) return;

        final EigenDecomposition eigen = substitutionModel.getEigenDecomposition();
        if (eigen == null) {
            Arrays.fill(matrixR,    0.0);
            Arrays.fill(matrixRInv, 0.0);
            Arrays.fill(eigenReal,  0.0);
            Arrays.fill(eigenImag,  0.0);
            blockCount = 0;
            eigenVersion++;
            eigenDirty = false;
            return;
        }

        loadEigenDecomposition(eigen, matrixR, matrixRInv, eigenReal, eigenImag, stateCount);
        buildBlockPlan(eigenReal, eigenImag, stateCount, blockType, blockStart);
        eigenVersion++;
        eigenDirty = false;
    }

    private void buildBlockPlan(double[] real, double[] imag, int dim,
                                byte[] outType, int[] outStart) {
        int b = 0, i = 0;
        while (i < dim) {
            if (Math.abs(imag[i]) <= IMAG_TOL) {
                outType[b]  = REAL_BLOCK;
                outStart[b] = i;
                b++; i++;
                continue;
            }
            if (i + 1 >= dim) throw new IllegalStateException(
                    "Malformed eigenvalue array: complex eigenvalue at final index " + i);
            final double a  = real[i];
            final double bi = imag[i];
            if (Math.abs(real[i + 1] - a) > PAIR_TOL) throw new IllegalStateException(
                    "Malformed eigenvalue array: real parts of conjugate pair do not match at " + i);
            if (Math.abs(imag[i + 1] + bi) > PAIR_TOL) throw new IllegalStateException(
                    "Malformed eigenvalue array: imaginary parts are not conjugates at " + i);
            outType[b]  = COMPLEX_BLOCK;
            outStart[b] = i;
            b++; i += 2;
        }
        blockCount = b;
    }

    private static void loadEigenDecomposition(EigenDecomposition eigen,
                                               double[] outR, double[] outRInv,
                                               double[] outReal, double[] outImag,
                                               int stateCount) {
        final double[] evec = eigen.getEigenVectors();
        final double[] ievc = eigen.getInverseEigenVectors();
        final double[] eval = eigen.getEigenValues();

        if (evec.length != stateCount * stateCount)
            throw new IllegalStateException("Unexpected eigenvector length " + evec.length);
        if (ievc.length != stateCount * stateCount)
            throw new IllegalStateException("Unexpected inverse eigenvector length " + ievc.length);

        System.arraycopy(evec, 0, outR,    0, outR.length);
        System.arraycopy(ievc, 0, outRInv, 0, outRInv.length);

        if (eval.length == stateCount) {
            System.arraycopy(eval, 0, outReal, 0, stateCount);
            Arrays.fill(outImag, 0.0);
        } else if (eval.length == 2 * stateCount) {
            System.arraycopy(eval, 0,          outReal, 0, stateCount);
            System.arraycopy(eval, stateCount, outImag, 0, stateCount);
        } else {
            throw new IllegalStateException("Unexpected eigenvalue array length " + eval.length);
        }
    }

    private void ensureBranchCoefficients(int nodeNumber, double branchLength) {
        if (nodeNumber >= nodeCapacity) {
            throw new IllegalArgumentException(
                    "nodeNumber " + nodeNumber + " exceeds nodeCapacity " + nodeCapacity);
        }

        final long   curVersion  = eigenVersion;
        final long   branchBits  = Double.doubleToLongBits(branchLength);
        final int    nodeSlotBase = nodeNumber * CACHE_SLOTS;

        // Fast path: scan all slots for this node.
        for (int s = 0; s < CACHE_SLOTS; s++) {
            final int idx = nodeSlotBase + s;
            if (cachedCoeffVersion[idx] == curVersion
                    && Double.doubleToLongBits(cachedBranchLength[idx]) == branchBits) {
                activeSlotByNode[nodeNumber] = s;
                return;
            }
        }

        // Cache miss: evict next round-robin slot.
        final int slot    = nextSlotByNode[nodeNumber];
        nextSlotByNode[nodeNumber]  = (slot + 1 == CACHE_SLOTS) ? 0 : slot + 1;
        activeSlotByNode[nodeNumber] = slot;

        final int coeffBase = (nodeSlotBase + slot) * stateCount;

        for (int b = 0; b < blockCount; b++) {
            final int i = blockStart[b];
            if (blockType[b] == REAL_BLOCK) {
                coeffA[coeffBase + i] = FastMath.exp(branchLength * eigenReal[i]);
                coeffB[coeffBase + i] = 0.0;
            } else {
                final double expat = FastMath.exp(branchLength * eigenReal[i]);
                final double bt    = branchLength * eigenImag[i];
                coeffA[coeffBase + i] = expat * Math.cos(bt);
                coeffB[coeffBase + i] = expat * Math.sin(bt);
            }
        }

        cachedBranchLength[nodeSlotBase + slot] = branchLength;
        cachedCoeffVersion[nodeSlotBase + slot]  = curVersion;
    }

    private void applyForwardCoefficients(int nodeNumber, double[] in, double[] out) {
        final int coeffBase = (nodeNumber * CACHE_SLOTS + activeSlotByNode[nodeNumber]) * stateCount;

        for (int b = 0; b < blockCount; b++) {
            final int i = blockStart[b];
            if (blockType[b] == REAL_BLOCK) {
                out[i] = coeffA[coeffBase + i] * in[i];
            } else {
                final double c = coeffA[coeffBase + i];
                final double s = coeffB[coeffBase + i];
                final double x = in[i];
                final double y = in[i + 1];
                out[i]     =  c * x + s * y;
                out[i + 1] = -s * x + c * y;
            }
        }
    }

    private void applyTransposeCoefficients(int nodeNumber, double[] in, double[] out) {
        final int coeffBase = (nodeNumber * CACHE_SLOTS + activeSlotByNode[nodeNumber]) * stateCount;

        for (int b = 0; b < blockCount; b++) {
            final int i = blockStart[b];
            if (blockType[b] == REAL_BLOCK) {
                out[i] = coeffA[coeffBase + i] * in[i];
            } else {
                final double c = coeffA[coeffBase + i];
                final double s = coeffB[coeffBase + i];
                final double x = in[i];
                final double y = in[i + 1];
                out[i]     = c * x - s * y;
                out[i + 1] = s * x + c * y;
            }
        }
    }

    private static void multiplyMatrixVector(double[] matrix, double[] vector,
                                             double[] out, int dim) {
        int base = 0;
        for (int i = 0; i < dim; i++) {
            double sum = 0.0;
            for (int j = 0; j < dim; j++) sum += matrix[base + j] * vector[j];
            out[i] = sum;
            base += dim;
        }
    }

    private static void multiplyTransposeMatrixVector(double[] matrix, double[] vector,
                                                      double[] out, int dim) {
        Arrays.fill(out, 0, dim, 0.0);
        int base = 0;
        for (int row = 0; row < dim; row++) {
            final double v = vector[row];
            for (int col = 0; col < dim; col++) out[col] += matrix[base + col] * v;
            base += dim;
        }
    }
}

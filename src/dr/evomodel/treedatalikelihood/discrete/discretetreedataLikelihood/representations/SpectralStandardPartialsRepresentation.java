package dr.evomodel.treedatalikelihood.discrete.discretetreedataLikelihood.representations;

import dr.evomodel.substmodel.BaseSubstitutionModel;
import dr.evomodel.substmodel.EigenDecomposition;
import org.apache.commons.math.util.FastMath;

import java.util.Arrays;

/**
 * Standard-basis partial representation with spectral branch propagation.
 *
 * Post-order : p_top    = P(t)   p_bottom   where P(t) = R exp(tD) R^{-1}
 * Pre-order  : q_bottom = P(t)^T q_top
 *
 * Messages are stored in the standard basis.
 * Propagation is performed through the eigensystem of the substitution model.
 *
 * Efficiency features (matching SpectralRotatedPartialsRepresentation):
 *   1. Flat 1-D primitive cache arrays (no Object[][] pointer-chasing).
 *   2. CACHE_SLOTS round-robin slots per node — all rate-categories coexist without eviction.
 *   3. Pre-computed R^T and R^{-T} (matrixRT / matrixRInvT) — avoids on-the-fly transposes.
 *   4. FastMath.exp for branch-coefficient computation.
 *   5. Real/complex block plan compiled once per eigen-invalidation.
 *
 * Complex eigenvalues follow the BEAST convention:
 *   eval[0..K-1]  = real parts,  eval[K..2K-1] = imaginary parts.
 *
 * @author Filippo Monti
 */
public final class SpectralStandardPartialsRepresentation
        implements BidirectionalRepresentation {

    private static final double IMAG_TOL = 1.0e-15;
    private static final double PAIR_TOL = 1.0e-12;

    private static final byte REAL_BLOCK    = 1;
    private static final byte COMPLEX_BLOCK = 2;

    private static final int CACHE_SLOTS = 16;

    private final BaseSubstitutionModel substitutionModel;
    private final int stateCount;

    private final double[] matrixR;
    private final double[] matrixRInv;
    private final double[] matrixRT;
    private final double[] matrixRInvT;
    private final double[] eigenReal;
    private final double[] eigenImag;

    private final byte[] blockType;
    private final int[]  blockStart;
    private int blockCount;

    private final double[] tmpA;
    private final double[] tmpB;

    private boolean eigenDirty = true;
    private long    eigenVersion = 0L;

    // Flat cache: index = (node * CACHE_SLOTS + slot) * stateCount + state
    private final double[] coeffA;
    private final double[] coeffB;
    // index = node * CACHE_SLOTS + slot
    private final double[] cachedBranchLength;
    private final long[]   cachedCoeffVersion;
    private final int[]    nextSlotByNode;
    private final int[]    activeSlotByNode;
    private final int      nodeCapacity;

    public SpectralStandardPartialsRepresentation(BaseSubstitutionModel substitutionModel,
                                                   int nodeCount) {
        if (substitutionModel == null) {
            throw new IllegalArgumentException("substitutionModel must not be null");
        }

        this.substitutionModel = substitutionModel;
        this.stateCount = substitutionModel.getDataType().getStateCount();

        this.matrixR    = new double[stateCount * stateCount];
        this.matrixRInv = new double[stateCount * stateCount];
        this.matrixRT   = new double[stateCount * stateCount];
        this.matrixRInvT = new double[stateCount * stateCount];
        this.eigenReal  = new double[stateCount];
        this.eigenImag  = new double[stateCount];

        this.blockType  = new byte[stateCount];
        this.blockStart = new int[stateCount];

        this.tmpA = new double[stateCount];
        this.tmpB = new double[stateCount];

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

    // ------------------------------------------------------------------
    // Shared metadata
    // ------------------------------------------------------------------

    @Override public String  getName()                  { return "spectral-standard-partials"; }
    @Override public int     getStateCount()            { return stateCount; }
    @Override public boolean supportsInternalScaling()  { return true; }
    @Override public boolean supportsSparseTipPartials() { return true; }

    // ------------------------------------------------------------------
    // Lifecycle
    // ------------------------------------------------------------------

    @Override public void markDirty()    { eigenDirty = true; }
    @Override public void storeState()   { /* caches are derived; nothing to snapshot */ }
    @Override public void restoreState() { eigenDirty = true; }
    @Override public void updateForLikelihood() { ensureEigenSystemCurrent(); }

    // ------------------------------------------------------------------
    // Post-order API
    // ------------------------------------------------------------------

    @Override
    public void initializeTipPartial(double[] standardTip, double[] outPartial) {
        System.arraycopy(standardTip, 0, outPartial, 0, stateCount);
    }

    @Override
    public void initializeSparseTipPartial(boolean[] stateSet, double[] outPartial) {
        initializeSparseTipPartial(stateSet, outPartial, 0);
    }

    @Override
    public void initializeSparseTipPartial(boolean[] stateSet, double[] outPartial, int outOffset) {
        for (int s = 0; s < stateCount; s++) {
            outPartial[outOffset + s] = stateSet[s] ? 1.0 : 0.0;
        }
    }

    @Override
    public void propagateSparseTipToBranchTop(int nodeNumber, double branchLength,
                                               boolean[] stateSet,
                                               double[] outBranchTopPartial, int outOffset) {
        ensureEigenSystemCurrent();
        ensureBranchCoefficients(nodeNumber, branchLength);

        // tmpA = R^{-1} * stateSet  (picks columns of R^{-1} for the observed states)
        multiplyMatrixStateSet(matrixRInv, stateSet, tmpA, 0, stateCount);
        // tmpB = exp(tD) tmpA
        applyForwardCoefficients(nodeNumber, tmpA, 0, tmpB, 0);
        // out = R tmpB
        multiplyMatrixVectorOffset(matrixR, tmpB, 0, outBranchTopPartial, outOffset, stateCount);
    }

    @Override
    public void propagateToBranchTop(int nodeNumber, double branchLength,
                                     double[] childPartial, double[] outBranchTopPartial) {
        propagateToBranchTop(nodeNumber, branchLength, childPartial, 0, outBranchTopPartial, 0);
    }

    @Override
    public void propagateToBranchTop(int nodeNumber, double branchLength,
                                     double[] childPartial, int childOffset,
                                     double[] outBranchTopPartial, int outOffset) {
        ensureEigenSystemCurrent();
        ensureBranchCoefficients(nodeNumber, branchLength);

        // tmpA = R^{-1} p
        multiplyMatrixVectorOffset(matrixRInv, childPartial, childOffset, tmpA, 0, stateCount);
        // tmpB = exp(tD) tmpA
        applyForwardCoefficients(nodeNumber, tmpA, 0, tmpB, 0);
        // out = R tmpB
        multiplyMatrixVectorOffset(matrixR, tmpB, 0, outBranchTopPartial, outOffset, stateCount);
    }

    @Override
    public void combineBranchTopPartials(double[] leftBranchTopPartial,
                                         double[] rightBranchTopPartial,
                                         double[] outParentPartial) {
        combineBranchTopPartials(leftBranchTopPartial, 0, rightBranchTopPartial, 0, outParentPartial, 0);
    }

    @Override
    public void combineBranchTopPartials(double[] leftBranchTopPartial, int leftOffset,
                                         double[] rightBranchTopPartial, int rightOffset,
                                         double[] outParentPartial, int outOffset) {
        for (int i = 0; i < stateCount; i++) {
            outParentPartial[outOffset + i] =
                    leftBranchTopPartial[leftOffset + i] * rightBranchTopPartial[rightOffset + i];
        }
    }

    @Override
    public double rootContribution(double[] rootFrequencies, double[] rootPartial) {
        return rootContribution(rootFrequencies, rootPartial, 0);
    }

    @Override
    public double rootContribution(double[] rootFrequencies, double[] rootPartial, int rootOffset) {
        double sum = 0.0;
        for (int i = 0; i < stateCount; i++) {
            sum += rootFrequencies[i] * rootPartial[rootOffset + i];
        }
        return sum;
    }

    // ------------------------------------------------------------------
    // Pre-order API
    // ------------------------------------------------------------------

    @Override
    public void initializeRootPartial(double[] rootFrequencies, double[] outPartial) {
        initializeRootPartial(rootFrequencies, outPartial, 0);
    }

    @Override
    public void initializeRootPartial(double[] rootFrequencies, double[] outPartial, int outOffset) {
        System.arraycopy(rootFrequencies, 0, outPartial, outOffset, stateCount);
    }

    @Override
    public void combineParentAndSibling(double[] parentPreOrderStandard, int parentOff,
                                        double[] siblingPostOrderStandard,
                                        double[] outChildBranchTopPreOrder, int outOff) {
        for (int i = 0; i < stateCount; i++) {
            outChildBranchTopPreOrder[outOff + i] =
                    parentPreOrderStandard[parentOff + i] * siblingPostOrderStandard[i];
        }
    }

    @Override
    public void propagateToBranchBottom(int nodeNumber, double branchLength,
                                        double[] branchTopPreOrder, int branchTopOffset,
                                        double[] out, int outOff) {
        ensureEigenSystemCurrent();
        ensureBranchCoefficients(nodeNumber, branchLength);

        // out = P(t)^T q  = R^{-T} exp(tD)^T R^T q
        // Step 1: tmpA = R^T q  (use pre-computed matrixRT for cache-friendly row access)
        multiplyMatrixVectorOffset(matrixRT, branchTopPreOrder, branchTopOffset, tmpA, 0, stateCount);
        // Step 2: tmpB = exp(tD)^T tmpA
        applyTransposeCoefficients(nodeNumber, tmpA, 0, tmpB, 0);
        // Step 3: out = R^{-T} tmpB  (use pre-computed matrixRInvT)
        multiplyMatrixVectorOffset(matrixRInvT, tmpB, 0, out, outOff, stateCount);
    }

    // ------------------------------------------------------------------
    // Export API
    // ------------------------------------------------------------------

    @Override
    public void exportPostOrderPartial(double[] partial, double[] outPartial) {
        exportPostOrderPartial(partial, 0, outPartial, 0);
    }

    @Override
    public void exportPostOrderPartial(double[] partial, int partialOffset,
                                       double[] outPartial, int outOffset) {
        System.arraycopy(partial, partialOffset, outPartial, outOffset, stateCount);
    }

    @Override
    public void exportPreOrderPartial(double[] partial, double[] outPartial) {
        exportPreOrderPartial(partial, 0, outPartial, 0);
    }

    @Override
    public void exportPreOrderPartial(double[] partial, int partialOffset,
                                      double[] outPartial, int outOffset) {
        System.arraycopy(partial, partialOffset, outPartial, outOffset, stateCount);
    }

    // ------------------------------------------------------------------
    // Eigensystem loading and block-plan compilation
    // ------------------------------------------------------------------

    private void ensureEigenSystemCurrent() {
        if (!eigenDirty) return;

        final EigenDecomposition eigen = substitutionModel.getEigenDecomposition();
        if (eigen == null) {
            Arrays.fill(matrixR,    0.0);
            Arrays.fill(matrixRInv, 0.0);
            Arrays.fill(matrixRT,   0.0);
            Arrays.fill(matrixRInvT, 0.0);
            Arrays.fill(eigenReal,  0.0);
            Arrays.fill(eigenImag,  0.0);
            blockCount = 0;
            eigenVersion++;
            eigenDirty = false;
            return;
        }

        loadEigenDecomposition(eigen, matrixR, matrixRInv, eigenReal, eigenImag, stateCount);
        transposeSquare(matrixR,    matrixRT,    stateCount);
        transposeSquare(matrixRInv, matrixRInvT, stateCount);
        buildBlockPlan(eigenReal, eigenImag, stateCount, blockType, blockStart);
        eigenVersion++;
        eigenDirty = false;
    }

    private void buildBlockPlan(double[] real, double[] imag, int dim,
                                byte[] outBlockType, int[] outBlockStart) {
        int b = 0, i = 0;
        while (i < dim) {
            if (Math.abs(imag[i]) <= IMAG_TOL) {
                outBlockType[b]  = REAL_BLOCK;
                outBlockStart[b] = i;
                b++; i++;
                continue;
            }
            if (i + 1 >= dim) throw new IllegalStateException(
                    "Malformed eigenvalue array: complex eigenvalue at final index " + i);
            final double a = real[i], bi = imag[i];
            if (Math.abs(real[i + 1] - a)  > PAIR_TOL) throw new IllegalStateException(
                    "Real parts of conjugate pair do not match at " + i);
            if (Math.abs(imag[i + 1] + bi) > PAIR_TOL) throw new IllegalStateException(
                    "Imaginary parts are not conjugates at " + i);
            outBlockType[b]  = COMPLEX_BLOCK;
            outBlockStart[b] = i;
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

        System.arraycopy(evec, 0, outR,    0, evec.length);
        System.arraycopy(ievc, 0, outRInv, 0, ievc.length);

        if (eval.length == stateCount) {
            System.arraycopy(eval, 0, outReal, 0, stateCount);
            Arrays.fill(outImag, 0.0);
        } else if (eval.length == 2 * stateCount) {
            System.arraycopy(eval, 0,          outReal, 0, stateCount);
            System.arraycopy(eval, stateCount,  outImag, 0, stateCount);
        } else {
            throw new IllegalStateException("Unexpected eigenvalue array length " + eval.length);
        }
    }

    private static void transposeSquare(double[] matrix, double[] transpose, int dim) {
        for (int i = 0; i < dim; i++) {
            final int rowBase = i * dim;
            for (int j = 0; j < dim; j++) {
                transpose[j * dim + i] = matrix[rowBase + j];
            }
        }
    }

    // ------------------------------------------------------------------
    // Flat branch-coefficient cache
    // ------------------------------------------------------------------

    private void ensureBranchCoefficients(int nodeNumber, double branchLength) {
        if (nodeNumber >= nodeCapacity) throw new IllegalArgumentException(
                "nodeNumber " + nodeNumber + " exceeds nodeCapacity " + nodeCapacity);

        final long curVersion = eigenVersion;
        final long branchBits = Double.doubleToLongBits(branchLength);
        final int  nodeSlotBase = nodeNumber * CACHE_SLOTS;

        for (int s = 0; s < CACHE_SLOTS; s++) {
            final int idx = nodeSlotBase + s;
            if (cachedCoeffVersion[idx] == curVersion
                    && Double.doubleToLongBits(cachedBranchLength[idx]) == branchBits) {
                activeSlotByNode[nodeNumber] = s;
                return;
            }
        }

        final int slot = nextSlotByNode[nodeNumber];
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

    // ------------------------------------------------------------------
    // Hot spectral application kernels (use flat cache via activeSlotByNode)
    // ------------------------------------------------------------------

    private void applyForwardCoefficients(int nodeNumber,
                                          double[] in,  int inOff,
                                          double[] out, int outOff) {
        final int coeffBase = (nodeNumber * CACHE_SLOTS + activeSlotByNode[nodeNumber]) * stateCount;

        for (int b = 0; b < blockCount; b++) {
            final int i = blockStart[b];
            if (blockType[b] == REAL_BLOCK) {
                out[outOff + i] = coeffA[coeffBase + i] * in[inOff + i];
            } else {
                final double c = coeffA[coeffBase + i];
                final double s = coeffB[coeffBase + i];
                final double x = in[inOff + i];
                final double y = in[inOff + i + 1];
                out[outOff + i]     =  c * x + s * y;
                out[outOff + i + 1] = -s * x + c * y;
            }
        }
    }

    private void applyTransposeCoefficients(int nodeNumber,
                                            double[] in,  int inOff,
                                            double[] out, int outOff) {
        final int coeffBase = (nodeNumber * CACHE_SLOTS + activeSlotByNode[nodeNumber]) * stateCount;

        for (int b = 0; b < blockCount; b++) {
            final int i = blockStart[b];
            if (blockType[b] == REAL_BLOCK) {
                out[outOff + i] = coeffA[coeffBase + i] * in[inOff + i];
            } else {
                final double c = coeffA[coeffBase + i];
                final double s = coeffB[coeffBase + i];
                final double x = in[inOff + i];
                final double y = in[inOff + i + 1];
                out[outOff + i]     = c * x - s * y;
                out[outOff + i + 1] = s * x + c * y;
            }
        }
    }

    // ------------------------------------------------------------------
    // Dense matrix-vector kernels
    // ------------------------------------------------------------------

    /**
     * out = matrix * stateSet, where stateSet is treated as a {0,1} indicator vector.
     * If exactly one state is set, copies the corresponding column of matrix (O(K)).
     * Otherwise sums the active columns (O(K * |stateSet|)).
     */
    private static void multiplyMatrixStateSet(double[] matrix, boolean[] stateSet,
                                               double[] out, int outOff, int dim) {
        // Fast path: single observed state → copy one column.
        int single = -1;
        for (int j = 0; j < dim; j++) {
            if (stateSet[j]) {
                if (single >= 0) { single = -1; break; }
                single = j;
            }
        }
        if (single >= 0) {
            for (int i = 0; i < dim; i++) out[outOff + i] = matrix[i * dim + single];
            return;
        }
        // General path: sum columns for all active states.
        Arrays.fill(out, outOff, outOff + dim, 0.0);
        int base = 0;
        for (int i = 0; i < dim; i++) {
            double sum = 0.0;
            for (int j = 0; j < dim; j++) { if (stateSet[j]) sum += matrix[base + j]; }
            out[outOff + i] = sum;
            base += dim;
        }
    }

    private static void multiplyMatrixVectorOffset(double[] matrix,
                                                   double[] vector, int vecOff,
                                                   double[] out,    int outOff,
                                                   int dim) {
        int base = 0;
        for (int i = 0; i < dim; i++) {
            double sum = 0.0;
            for (int j = 0; j < dim; j++) sum += matrix[base + j] * vector[vecOff + j];
            out[outOff + i] = sum;
            base += dim;
        }
    }
}

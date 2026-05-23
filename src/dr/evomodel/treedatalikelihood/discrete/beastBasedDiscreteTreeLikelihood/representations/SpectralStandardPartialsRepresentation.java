package dr.evomodel.treedatalikelihood.discrete.beastBasedDiscreteTreeLikelihood.representations;

import dr.evomodel.substmodel.BaseSubstitutionModel;
import dr.evomodel.substmodel.EigenDecomposition;

import java.util.Arrays;

/**
 * Standard-basis partial representation with spectral branch propagation.
 *
 * This class implements both post-order and pre-order traversals:
 *
 *   post-order : p_top    = P(t)     p_bottom
 *   pre-order  : q_bottom = P(t)^T   q_top
 *
 * where
 *
 *   P(t) = R exp(tD) R^{-1}.
 *
 * Messages are always stored in the standard basis.
 * Propagation is done through the eigensystem of the substitution model.
 *
 * Efficiency features:
 *
 * 1. One shared implementation for post- and pre-order.
 * 2. Eigensystem copied once per invalidation.
 * 3. Real/complex block layout parsed once per invalidation.
 * 4. Per-branch spectral coefficients cached and reused by both directions.
 * 5. Hot application kernels contain only multiply/add operations.
 *
 * Complex eigenvalues follow the BEAST convention:
 *
 *   eval[0 .. K-1]   = real parts
 *   eval[K .. 2K-1]  = imaginary parts
 */
/*
 * @author Filippo Monti
 */
public final class SpectralStandardPartialsRepresentation
        implements BidirectionalRepresentation {

    private static final double IMAG_TOL = 1.0e-15;
    private static final double PAIR_TOL = 1.0e-12;

    private static final byte REAL_BLOCK = 1;
    private static final byte COMPLEX_BLOCK = 2;

    private final BaseSubstitutionModel substitutionModel;
    private final int stateCount;

    // Current eigensystem cached locally
    private final double[] matrixR;
    private final double[] matrixRInv;
    private final double[] eigenReal;
    private final double[] eigenImag;

    // Block plan for exp(tD) action in transformed coordinates
    // blockType[b] in {REAL_BLOCK, COMPLEX_BLOCK}
    // blockStart[b] = starting index in transformed coordinates
    private final byte[] blockType;
    private final int[] blockStart;
    private int blockCount;

    // Scratch reused across calls
    private final double[] tmpA;
    private final double[] tmpB;

    // Eigensystem versioning
    private boolean eigenDirty = true;
    private long eigenVersion = 0L;

    // Per-node branch coefficient cache.
    // For a real block at index i:
    //   coeffA[node][i] = exp(t * real[i])
    //   coeffB[node][i] unused (= 0)
    //
    // For a complex pair block starting at i:
    //   coeffA[node][i] = exp(ta) cos(tb)
    //   coeffB[node][i] = exp(ta) sin(tb)
    //   slot i+1 is unused
    //
    // The cache is valid when:
    //   cachedCoeffVersion[node] == eigenVersion
    //   && cachedBranchLength[node] == branchLength
    private double[][] coeffAByNode;
    private double[][] coeffBByNode;
    private double[] cachedBranchLength;
    private long[] cachedCoeffVersion;

    public SpectralStandardPartialsRepresentation(BaseSubstitutionModel substitutionModel) {
        if (substitutionModel == null) {
            throw new IllegalArgumentException("substitutionModel must not be null");
        }

        this.substitutionModel = substitutionModel;
        this.stateCount = substitutionModel.getDataType().getStateCount();

        this.matrixR = new double[stateCount * stateCount];
        this.matrixRInv = new double[stateCount * stateCount];
        this.eigenReal = new double[stateCount];
        this.eigenImag = new double[stateCount];

        this.blockType = new byte[stateCount];
        this.blockStart = new int[stateCount];

        this.tmpA = new double[stateCount];
        this.tmpB = new double[stateCount];

        this.coeffAByNode = new double[0][];
        this.coeffBByNode = new double[0][];
        this.cachedBranchLength = new double[0];
        this.cachedCoeffVersion = new long[0];
    }

    // ------------------------------------------------------------------
    // Shared metadata
    // ------------------------------------------------------------------

    @Override
    public String getName() {
        return "spectral-standard-partials";
    }

    @Override
    public int getStateCount() {
        return stateCount;
    }

    @Override
    public boolean supportsInternalScaling() {
        return true;
    }

    // ------------------------------------------------------------------
    // Lifecycle
    // ------------------------------------------------------------------

    @Override
    public void markDirty() {
        eigenDirty = true;
    }

    @Override
    public void storeState() {
        // No explicit stored/restored branch state here.
        // All cached quantities are derived from the current substitution model state.
    }

    @Override
    public void restoreState() {
        // Conservative and safe.
        eigenDirty = true;
    }

    @Override
    public void updateForLikelihood() {
        ensureEigenSystemCurrent();
    }

    // ------------------------------------------------------------------
    // Post-order API
    // ------------------------------------------------------------------

    @Override
    public void initializeTipPartial(double[] standardTip, double[] outPartial) {
        initializeTipPartial(standardTip, outPartial, 0);
    }

    @Override
    public void initializeTipPartial(double[] standardTip, double[] outPartial, int outOffset) {
        for (int s = 0; s < stateCount; s++) {
            outPartial[outOffset + s] = standardTip[s];
        }
    }

    @Override
    public void propagateToBranchTop(int nodeNumber,
                                     double branchLength,
                                     double[] childPartial,
                                     double[] outBranchTopPartial) {
        propagateToBranchTop(nodeNumber, branchLength, childPartial, 0, outBranchTopPartial, 0);
    }

    @Override
    public void propagateToBranchTop(int nodeNumber,
                                     double branchLength,
                                     double[] childPartial,
                                     int childOffset,
                                     double[] outBranchTopPartial,
                                     int outOffset) {

        ensureEigenSystemCurrent();
        ensureBranchCoefficients(nodeNumber, branchLength);

        // tmpA = R^{-1} p
        multiplyMatrixVectorOffset(matrixRInv, childPartial, childOffset, tmpA, 0, stateCount);

        // tmpB = exp(tD) tmpA
        applyForwardCoefficients(nodeNumber, tmpA, tmpB);

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
    public void combineBranchTopPartials(double[] leftBranchTopPartial,
                                         int leftOffset,
                                         double[] rightBranchTopPartial,
                                         int rightOffset,
                                         double[] outParentPartial,
                                         int outOffset) {
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
        for (int s = 0; s < stateCount; s++) {
            outPartial[outOffset + s] = rootFrequencies[s];
        }
    }

    @Override
    public void combineParentAndSibling(double[] parentPreOrderStandard, int parentOff,
                                        double[] siblingPostOrderStandard,
                                        double[] outChildBranchTopPreOrder,
                                        int outOff) {
        for (int i = 0; i < stateCount; i++) {
            outChildBranchTopPreOrder[outOff + i] =
                    parentPreOrderStandard[parentOff + i] * siblingPostOrderStandard[i];
        }
    }

    @Override
    public void propagateToBranchBottom(int nodeNumber,
                                        double branchLength,
                                        double[] branchTopPreOrder,
                                        int branchTopOffset,
                                        double[] out, int outOff) {

        ensureEigenSystemCurrent();
        ensureBranchCoefficients(nodeNumber, branchLength);

        multiplyTransposeMatrixVectorOffset(matrixR, branchTopPreOrder, branchTopOffset, tmpA, 0, stateCount);
        applyTransposeCoefficients(nodeNumber, tmpA, tmpB);
        multiplyTransposeMatrixVectorOffset(matrixRInv, tmpB, 0, out, outOff, stateCount);
    }

    // ------------------------------------------------------------------
    // Common export
    // ------------------------------------------------------------------

    @Override
    public void exportPostOrderPartial(double[] partial, double[] outPartial) {
        exportPostOrderPartial(partial, 0, outPartial, 0);
    }

    @Override
    public void exportPostOrderPartial(double[] partial, int partialOffset,
                                       double[] outPartial, int outOffset) {
        for (int s = 0; s < stateCount; s++) {
            outPartial[outOffset + s] = partial[partialOffset + s];
        }
    }

    @Override
    public void exportPreOrderPartial(double[] partial, double[] outPartial) {
        exportPreOrderPartial(partial, 0, outPartial, 0);
    }

    @Override
    public void exportPreOrderPartial(double[] partial, int partialOffset,
                                      double[] outPartial, int outOffset) {
        for (int s = 0; s < stateCount; s++) {
            outPartial[outOffset + s] = partial[partialOffset + s];
        }
    }

    // ------------------------------------------------------------------
    // Eigensystem loading and block-plan compilation
    // ------------------------------------------------------------------
    private void ensureEigenSystemCurrent() {
        if (!eigenDirty) {
            return;
        }

        final EigenDecomposition eigen = substitutionModel.getEigenDecomposition();
        if (eigen == null) {
            Arrays.fill(matrixR, 0.0);
            Arrays.fill(matrixRInv, 0.0);
            Arrays.fill(eigenReal, 0.0);
            Arrays.fill(eigenImag, 0.0);
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

    private void buildBlockPlan(double[] real,
                                double[] imag,
                                int dim,
                                byte[] outBlockType,
                                int[] outBlockStart) {
        int b = 0;
        int i = 0;

        while (i < dim) {
            if (Math.abs(imag[i]) <= IMAG_TOL) {
                outBlockType[b] = REAL_BLOCK;
                outBlockStart[b] = i;
                b++;
                i++;
                continue;
            }

            if (i + 1 >= dim) {
                throw new IllegalStateException(
                        "Malformed eigenvalue array: complex eigenvalue at final index " + i);
            }

            final double a = real[i];
            final double bi = imag[i];

            if (Math.abs(real[i + 1] - a) > PAIR_TOL) {
                throw new IllegalStateException(
                        "Malformed eigenvalue array: real parts of conjugate pair do not match at indices "
                                + i + " and " + (i + 1));
            }

            if (Math.abs(imag[i + 1] + bi) > PAIR_TOL) {
                throw new IllegalStateException(
                        "Malformed eigenvalue array: imaginary parts are not conjugates at indices "
                                + i + " and " + (i + 1));
            }

            outBlockType[b] = COMPLEX_BLOCK;
            outBlockStart[b] = i;
            b++;
            i += 2;
        }

        blockCount = b;
    }

    private static void loadEigenDecomposition(EigenDecomposition eigen,
                                               double[] outR,
                                               double[] outRInv,
                                               double[] outEigenReal,
                                               double[] outEigenImag,
                                               int stateCount) {

        final double[] evec = eigen.getEigenVectors();
        final double[] ievc = eigen.getInverseEigenVectors();
        final double[] eval = eigen.getEigenValues();

        if (evec.length != stateCount * stateCount) {
            throw new IllegalStateException(
                    "Unexpected eigenvector length " + evec.length +
                            " for stateCount = " + stateCount);
        }

        if (ievc.length != stateCount * stateCount) {
            throw new IllegalStateException(
                    "Unexpected inverse eigenvector length " + ievc.length +
                            " for stateCount = " + stateCount);
        }

        for (int i = 0; i < outR.length; i++) {
            outR[i] = evec[i];
            outRInv[i] = ievc[i];
        }

        if (eval.length == stateCount) {
            for (int i = 0; i < stateCount; i++) {
                outEigenReal[i] = eval[i];
            }
            Arrays.fill(outEigenImag, 0.0);
            return;
        }

        if (eval.length == 2 * stateCount) {
            for (int i = 0; i < stateCount; i++) {
                outEigenReal[i] = eval[i];
                outEigenImag[i] = eval[stateCount + i];
            }
            return;
        }

        throw new IllegalStateException(
                "Unexpected eigenvalue array length " + eval.length +
                        " for stateCount = " + stateCount +
                        ". Expected either K or 2K.");
    }

    // ------------------------------------------------------------------
    // Branch coefficient cache
    // ------------------------------------------------------------------

    private void ensureNodeCapacity(int nodeNumber) {
        if (nodeNumber < coeffAByNode.length) {
            return;
        }

        final int newSize = Math.max(nodeNumber + 1, Math.max(4, coeffAByNode.length * 2));

        coeffAByNode = Arrays.copyOf(coeffAByNode, newSize);
        coeffBByNode = Arrays.copyOf(coeffBByNode, newSize);
        cachedBranchLength = Arrays.copyOf(cachedBranchLength, newSize);
        cachedCoeffVersion = Arrays.copyOf(cachedCoeffVersion, newSize);

        // Mark new entries invalid
        for (int i = nodeNumber; i < newSize; i++) {
            if (coeffAByNode[i] == null) {
                cachedCoeffVersion[i] = Long.MIN_VALUE;
                cachedBranchLength[i] = Double.NaN;
            }
        }
    }

    private void ensureBranchCoefficients(int nodeNumber, double branchLength) {
        if (nodeNumber < 0) {
            throw new IllegalArgumentException("nodeNumber must be nonnegative");
        }

        ensureNodeCapacity(nodeNumber);

        if (cachedCoeffVersion[nodeNumber] == eigenVersion
                && Double.doubleToLongBits(cachedBranchLength[nodeNumber]) == Double.doubleToLongBits(branchLength)) {
            return;
        }

        double[] coeffA = coeffAByNode[nodeNumber];
        double[] coeffB = coeffBByNode[nodeNumber];

        if (coeffA == null) {
            coeffA = new double[stateCount];
            coeffAByNode[nodeNumber] = coeffA;
        }
        if (coeffB == null) {
            coeffB = new double[stateCount];
            coeffBByNode[nodeNumber] = coeffB;
        }

        // Only block starts matter. Slots that are not block starts are ignored.
        for (int b = 0; b < blockCount; b++) {
            final int i = blockStart[b];

            if (blockType[b] == REAL_BLOCK) {
                coeffA[i] = Math.exp(branchLength * eigenReal[i]);
                coeffB[i] = 0.0;
            } else {
                final double a = eigenReal[i];
                final double beta = eigenImag[i];
                final double expat = Math.exp(branchLength * a);
                final double bt = branchLength * beta;
                coeffA[i] = expat * Math.cos(bt);
                coeffB[i] = expat * Math.sin(bt);
            }
        }

        cachedBranchLength[nodeNumber] = branchLength;
        cachedCoeffVersion[nodeNumber] = eigenVersion;
    }

    // ------------------------------------------------------------------
    // Hot spectral application kernels
    // ------------------------------------------------------------------

    private void applyForwardCoefficients(int nodeNumber, double[] in, double[] out) {
        final double[] coeffA = coeffAByNode[nodeNumber];
        final double[] coeffB = coeffBByNode[nodeNumber];

        for (int b = 0; b < blockCount; b++) {
            final int i = blockStart[b];

            if (blockType[b] == REAL_BLOCK) {
                out[i] = coeffA[i] * in[i];
            } else {
                final double c = coeffA[i];
                final double s = coeffB[i];
                final double x = in[i];
                final double y = in[i + 1];

                out[i]     = c * x + s * y;
                out[i + 1] = -s * x + c * y;
            }
        }
    }

    private void applyTransposeCoefficients(int nodeNumber, double[] in, double[] out) {
        final double[] coeffA = coeffAByNode[nodeNumber];
        final double[] coeffB = coeffBByNode[nodeNumber];

        for (int b = 0; b < blockCount; b++) {
            final int i = blockStart[b];

            if (blockType[b] == REAL_BLOCK) {
                out[i] = coeffA[i] * in[i];
            } else {
                final double c = coeffA[i];
                final double s = coeffB[i];
                final double x = in[i];
                final double y = in[i + 1];

                out[i]     = c * x - s * y;
                out[i + 1] = s * x + c * y;
            }
        }
    }

    // ------------------------------------------------------------------
    // Dense matrix-vector kernels
    // ------------------------------------------------------------------

    /**
     * out = matrix^T * vector
     *
     * matrix stored row-major, dimension dim x dim.
     */
    private static void multiplyTransposeMatrixVector(double[] matrix,
                                                      double[] vector,
                                                      double[] out,
                                                      int dim) {
        Arrays.fill(out, 0, dim, 0.0);

        int base = 0;
        for (int row = 0; row < dim; row++) {
            final double v = vector[row];
            for (int col = 0; col < dim; col++) {
                out[col] += matrix[base + col] * v;
            }
            base += dim;
        }
    }

    private static void multiplyMatrixVectorOffset(double[] matrix,
                                                   double[] vector, int vectorOff,
                                                   double[] out, int outOff,
                                                   int dim) {
        int base = 0;
        for (int i = 0; i < dim; i++) {
            double sum = 0.0;
            for (int j = 0; j < dim; j++) {
                sum += matrix[base + j] * vector[vectorOff + j];
            }
            out[outOff + i] = sum;
            base += dim;
        }
    }

    private static void multiplyTransposeMatrixVectorOffset(double[] matrix,
                                                            double[] vector, int vectorOff,
                                                            double[] out, int outOff,
                                                            int dim) {
        Arrays.fill(out, outOff, outOff + dim, 0.0);

        int base = 0;
        for (int row = 0; row < dim; row++) {
            final double v = vector[vectorOff + row];
            for (int col = 0; col < dim; col++) {
                out[outOff + col] += matrix[base + col] * v;
            }
            base += dim;
        }
    }
}

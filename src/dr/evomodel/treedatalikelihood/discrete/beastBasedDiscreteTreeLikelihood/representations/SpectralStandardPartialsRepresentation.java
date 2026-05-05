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
        eigenDirty = true;
        ensureEigenSystemCurrent();
    }

    // ------------------------------------------------------------------
    // Post-order API
    // ------------------------------------------------------------------

    @Override
    public void initializeTipPartial(double[] standardTip, double[] outPartial) {
        System.arraycopy(standardTip, 0, outPartial, 0, stateCount);
    }

    @Override
    public void propagateToBranchTop(int nodeNumber,
                                     double branchLength,
                                     double[] childPartial,
                                     double[] outBranchTopPartial) {

        ensureEigenSystemCurrent();
        ensureBranchCoefficients(nodeNumber, branchLength);

        // tmpA = R^{-1} p
        multiplyMatrixVector(matrixRInv, childPartial, tmpA, stateCount);

        // tmpB = exp(tD) tmpA
        applyForwardCoefficients(nodeNumber, tmpA, tmpB);

        // out = R tmpB
        multiplyMatrixVector(matrixR, tmpB, outBranchTopPartial, stateCount);
    }

    @Override
    public void combineBranchTopPartials(double[] leftBranchTopPartial,
                                         double[] rightBranchTopPartial,
                                         double[] outParentPartial) {
        for (int i = 0; i < stateCount; i++) {
            outParentPartial[i] = leftBranchTopPartial[i] * rightBranchTopPartial[i];
        }
    }

    @Override
    public double rootContribution(double[] rootFrequencies, double[] rootPartial) {
        double sum = 0.0;
        for (int i = 0; i < stateCount; i++) {
            sum += rootFrequencies[i] * rootPartial[i];
        }
        return sum;
    }

    // ------------------------------------------------------------------
    // Pre-order API
    // ------------------------------------------------------------------

    @Override
    public void initializeRootPartial(double[] rootFrequencies, double[] outPartial) {
        System.arraycopy(rootFrequencies, 0, outPartial, 0, stateCount);
    }

    @Override
    public void combineParentAndSibling(double[] parentPreOrderStandard,
                                        double[] siblingPostOrderStandard,
                                        double[] outChildBranchTopPreOrder) {
        for (int i = 0; i < stateCount; i++) {
            outChildBranchTopPreOrder[i] =
                    parentPreOrderStandard[i] * siblingPostOrderStandard[i];
        }
    }

    @Override
    public void propagateToBranchBottom(int nodeNumber,
                                        double branchLength,
                                        double[] branchTopPreOrder,
                                        double[] outBranchBottomPreOrder) {

        ensureEigenSystemCurrent();
        ensureBranchCoefficients(nodeNumber, branchLength);

        // tmpA = R^T q_top
        multiplyTransposeMatrixVector(matrixR, branchTopPreOrder, tmpA, stateCount);

        // tmpB = exp(tD)^T tmpA
        applyTransposeCoefficients(nodeNumber, tmpA, tmpB);

        // out = (R^{-1})^T tmpB
        multiplyTransposeMatrixVector(matrixRInv, tmpB, outBranchBottomPreOrder, stateCount);
    }

    // ------------------------------------------------------------------
    // Common standard-basis conversion
    // ------------------------------------------------------------------

    @Override
    public void toStandard(double[] partial, double[] outStandardPartial) {
        System.arraycopy(partial, 0, outStandardPartial, 0, stateCount);
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

        System.arraycopy(evec, 0, outR, 0, outR.length);
        System.arraycopy(ievc, 0, outRInv, 0, outRInv.length);

        if (eval.length == stateCount) {
            System.arraycopy(eval, 0, outEigenReal, 0, stateCount);
            Arrays.fill(outEigenImag, 0.0);
            return;
        }

        if (eval.length == 2 * stateCount) {
            System.arraycopy(eval, 0, outEigenReal, 0, stateCount);
            System.arraycopy(eval, stateCount, outEigenImag, 0, stateCount);
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
     * out = matrix * vector
     *
     * matrix stored row-major, dimension dim x dim.
     */
    private static void multiplyMatrixVector(double[] matrix,
                                             double[] vector,
                                             double[] out,
                                             int dim) {
        int base = 0;
        for (int i = 0; i < dim; i++) {
            double sum = 0.0;
            for (int j = 0; j < dim; j++) {
                sum += matrix[base + j] * vector[j];
            }
            out[i] = sum;
            base += dim;
        }
    }

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
}

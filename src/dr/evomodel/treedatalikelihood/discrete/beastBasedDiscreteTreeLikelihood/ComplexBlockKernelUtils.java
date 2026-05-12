package dr.evomodel.treedatalikelihood.discrete.beastBasedDiscreteTreeLikelihood;

import dr.evomodel.substmodel.EigenDecomposition;
import org.apache.commons.math.util.FastMath;

public final class ComplexBlockKernelUtils {

    private static final double EIGEN_TOLERANCE = 1.0e-12;
    private static final double[] PADE13 = new double[] {
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
    private static final double THETA13 = 5.371920351148152;

    private static final byte ENTRY_SCALAR = 1;
    private static final byte ENTRY_ONE_BY_TWO = 2;
    private static final byte ENTRY_TWO_BY_ONE = 3;
    private static final byte ENTRY_TWO_BY_TWO = 4;
    private static final byte ENTRY_GENERAL = 5;

    private ComplexBlockKernelUtils() { }

    // -------------------------------------------------------------------------
    // Public flat plan type — all storage pre-allocated at construction time.
    // Call fillPlan to populate; call applyPlan to use. Zero dynamic allocation.
    // -------------------------------------------------------------------------
    public static final class ComplexKernelPlan {
        // scratch used inside fillPlan (block structure and per-block exp/trig values)
        final int[]    blockStarts;  // [stateCount]
        final int[]    blockDims;    // [stateCount]
        final double[] expR;         // [stateCount]  exp(t * real_eigenvalue) per block
        final double[] cosB;         // [stateCount]  cos(t * imag_eigenvalue) per block
        final double[] sinB;         // [stateCount]  sin(t * imag_eigenvalue) per block

        // per-entry metadata; at most stateCount^2 entries
        int blockCount;
        int entryCount;
        final int[] eLeftStart;   // [stateCount^2]
        final int[] eLeftDim;     // [stateCount^2]
        final int[] eLeftBlock;   // [stateCount^2]
        final int[] eRightStart;  // [stateCount^2]
        final int[] eRightDim;    // [stateCount^2]
        final int[] eRightBlock;  // [stateCount^2]
        final byte[] eKind;       // [stateCount^2]
        final double[] eRealShift;
        final double[] eLeftImag;
        final double[] eRightImag;
        final double[] eInvDenom;
        final double[] ePlusImagShift;
        final double[] eMinusImagShift;
        final double[] eInvPlusDenom;
        final double[] eInvMinusDenom;

        // flat coefficient storage: entry e uses slots [e*16 .. e*16+15]
        // 16 is the max coefficient count for a 2×2 block pair
        final double[] coefficients;  // [stateCount^2 * 16]

        // Entries are stored grouped by kind so type-specific loops need no per-entry branch.
        // [scalarStart, scalarEnd)       → ENTRY_SCALAR
        // [oneByTwoStart, oneByTwoEnd)   → ENTRY_ONE_BY_TWO
        // [twoByOneStart, twoByOneEnd)   → ENTRY_TWO_BY_ONE
        // [twoByTwoStart, twoByTwoEnd)   → ENTRY_TWO_BY_TWO
        // [generalStart,  generalEnd)    → ENTRY_GENERAL
        int scalarStart, scalarEnd;
        int oneByTwoStart, oneByTwoEnd;
        int twoByOneStart, twoByOneEnd;
        int twoByTwoStart, twoByTwoEnd;
        int generalStart, generalEnd;

        public ComplexKernelPlan(int stateCount) {
            final int maxPairs = stateCount * stateCount;
            this.blockStarts  = new int[stateCount];
            this.blockDims    = new int[stateCount];
            this.expR         = new double[stateCount];
            this.cosB         = new double[stateCount];
            this.sinB         = new double[stateCount];
            this.eLeftStart   = new int[maxPairs];
            this.eLeftDim     = new int[maxPairs];
            this.eLeftBlock   = new int[maxPairs];
            this.eRightStart  = new int[maxPairs];
            this.eRightDim    = new int[maxPairs];
            this.eRightBlock  = new int[maxPairs];
            this.eKind        = new byte[maxPairs];
            this.eRealShift   = new double[maxPairs];
            this.eLeftImag    = new double[maxPairs];
            this.eRightImag   = new double[maxPairs];
            this.eInvDenom    = new double[maxPairs];
            this.ePlusImagShift  = new double[maxPairs];
            this.eMinusImagShift = new double[maxPairs];
            this.eInvPlusDenom   = new double[maxPairs];
            this.eInvMinusDenom  = new double[maxPairs];
            this.coefficients = new double[maxPairs * 16];
        }
    }

    /**
     * Fills a pre-allocated plan in-place — zero heap allocation.
     * The plan must have been constructed with the same stateCount.
     */
    public static void fillPlan(ComplexKernelPlan plan,
                                EigenDecomposition decomposition,
                                double time,
                                int stateCount) {
        fillStructure(plan, decomposition, stateCount);
        fillTimeDependentCoefficients(plan, decomposition, time, stateCount);
    }

    public static void fillStructure(ComplexKernelPlan plan,
                                     EigenDecomposition decomposition,
                                     int stateCount) {
        final double[] eigenValues = decomposition.getEigenValues();
        final int blockCount = buildEigenBlocks(eigenValues, plan.blockStarts, plan.blockDims, stateCount);
        plan.blockCount = blockCount;

        // Fill entries grouped by kind so downstream loops need no per-entry branch dispatch.
        int entryIndex = 0;
        plan.scalarStart   = entryIndex;
        entryIndex = fillEntriesByKind(plan, eigenValues, stateCount, blockCount, entryIndex, ENTRY_SCALAR);
        plan.scalarEnd     = entryIndex;
        plan.oneByTwoStart = entryIndex;
        entryIndex = fillEntriesByKind(plan, eigenValues, stateCount, blockCount, entryIndex, ENTRY_ONE_BY_TWO);
        plan.oneByTwoEnd   = entryIndex;
        plan.twoByOneStart = entryIndex;
        entryIndex = fillEntriesByKind(plan, eigenValues, stateCount, blockCount, entryIndex, ENTRY_TWO_BY_ONE);
        plan.twoByOneEnd   = entryIndex;
        plan.twoByTwoStart = entryIndex;
        entryIndex = fillEntriesByKind(plan, eigenValues, stateCount, blockCount, entryIndex, ENTRY_TWO_BY_TWO);
        plan.twoByTwoEnd   = entryIndex;
        plan.generalStart  = entryIndex;
        entryIndex = fillEntriesByKind(plan, eigenValues, stateCount, blockCount, entryIndex, ENTRY_GENERAL);
        plan.generalEnd    = entryIndex;
        plan.entryCount    = entryIndex;
    }

    private static int fillEntriesByKind(ComplexKernelPlan plan,
                                         double[] eigenValues,
                                         int stateCount,
                                         int blockCount,
                                         int entryIndex,
                                         byte kind) {
        for (int leftBlock = 0; leftBlock < blockCount; ++leftBlock) {
            final int leftStart = plan.blockStarts[leftBlock];
            final int leftDim   = plan.blockDims[leftBlock];
            for (int rightBlock = 0; rightBlock < blockCount; ++rightBlock) {
                final int rightStart = plan.blockStarts[rightBlock];
                final int rightDim   = plan.blockDims[rightBlock];

                final byte entryKind;
                if      (leftDim == 1 && rightDim == 1) entryKind = ENTRY_SCALAR;
                else if (leftDim == 1 && rightDim == 2) entryKind = ENTRY_ONE_BY_TWO;
                else if (leftDim == 2 && rightDim == 1) entryKind = ENTRY_TWO_BY_ONE;
                else if (leftDim == 2 && rightDim == 2) entryKind = ENTRY_TWO_BY_TWO;
                else                                    entryKind = ENTRY_GENERAL;

                if (entryKind != kind) continue;

                plan.eLeftStart[entryIndex]  = leftStart;
                plan.eLeftDim[entryIndex]    = leftDim;
                plan.eLeftBlock[entryIndex]  = leftBlock;
                plan.eRightStart[entryIndex] = rightStart;
                plan.eRightDim[entryIndex]   = rightDim;
                plan.eRightBlock[entryIndex] = rightBlock;
                plan.eKind[entryIndex]       = kind;
                fillEntryMetadata(plan, entryIndex, eigenValues, stateCount,
                        leftStart, leftDim, rightStart, rightDim);
                entryIndex++;
            }
        }
        return entryIndex;
    }

    public static void fillTimeDependentCoefficients(ComplexKernelPlan plan,
                                                     EigenDecomposition decomposition,
                                                     double time,
                                                     int stateCount) {
        final double[] eigenValues = decomposition.getEigenValues();

        // Pre-compute exp(t*real) and, for complex blocks, cos/sin(t*imag) once per block.
        for (int b = 0; b < plan.blockCount; ++b) {
            final int start = plan.blockStarts[b];
            plan.expR[b] = FastMath.exp(time * eigenValues[start]);
            if (plan.blockDims[b] == 2) {
                final double imag = eigenValues[start + stateCount];
                plan.cosB[b] = Math.cos(time * imag);
                plan.sinB[b] = Math.sin(time * imag);
            }
        }

        // Type-specific tight loops — no per-entry kind dispatch.
        final double[] c = plan.coefficients;

        for (int e = plan.scalarStart; e < plan.scalarEnd; e++) {
            fillScalarCoefficient(c, e * 16, plan.expR[plan.eLeftBlock[e]],
                    plan.expR[plan.eRightBlock[e]], plan.eInvDenom[e], time);
        }

        for (int e = plan.oneByTwoStart; e < plan.oneByTwoEnd; e++) {
            final int rb = plan.eRightBlock[e];
            fillOneByTwoCoefficients(c, e * 16,
                    plan.expR[plan.eLeftBlock[e]], plan.expR[rb],
                    plan.cosB[rb], plan.sinB[rb],
                    plan.eRealShift[e], plan.eRightImag[e], plan.eInvDenom[e], time);
        }

        for (int e = plan.twoByOneStart; e < plan.twoByOneEnd; e++) {
            final int lb = plan.eLeftBlock[e];
            fillTwoByOneCoefficients(c, e * 16,
                    plan.expR[lb], plan.expR[plan.eRightBlock[e]],
                    plan.cosB[lb], plan.sinB[lb],
                    plan.eRealShift[e], plan.eLeftImag[e], plan.eInvDenom[e], time);
        }

        for (int e = plan.twoByTwoStart; e < plan.twoByTwoEnd; e++) {
            final int lb = plan.eLeftBlock[e];
            final int rb = plan.eRightBlock[e];
            fillTwoByTwoCoefficients(c, e * 16,
                    plan.expR[lb], plan.expR[rb],
                    plan.cosB[lb], plan.sinB[lb],
                    plan.cosB[rb], plan.sinB[rb],
                    plan.eRealShift[e],
                    plan.ePlusImagShift[e], plan.eMinusImagShift[e],
                    plan.eInvPlusDenom[e], plan.eInvMinusDenom[e], time);
        }

        for (int e = plan.generalStart; e < plan.generalEnd; e++) {
            fillGeneralCoefficients(c, e * 16, eigenValues,
                    plan.eLeftStart[e], plan.eLeftDim[e],
                    plan.eRightStart[e], plan.eRightDim[e],
                    stateCount, time);
        }
    }

    public static void applyPlan(ComplexKernelPlan plan,
                                 double[] transformed,
                                 double[] eigenBasisGradient,
                                 Workspace workspace,
                                 int stateCount) {
        final double[] c = plan.coefficients;

        // SCALAR
        for (int e = plan.scalarStart; e < plan.scalarEnd; e++) {
            final int leftStart  = plan.eLeftStart[e];
            final int rightStart = plan.eRightStart[e];
            eigenBasisGradient[leftStart * stateCount + rightStart] +=
                    c[e * 16] * transformed[leftStart * stateCount + rightStart];
        }

        // ONE_BY_TWO: scalar left, 2-wide right
        for (int e = plan.oneByTwoStart; e < plan.oneByTwoEnd; e++) {
            final int leftStart  = plan.eLeftStart[e];
            final int rightStart = plan.eRightStart[e];
            final int coeffBase  = e * 16;
            final double in0 = transformed[leftStart * stateCount + rightStart];
            final double in1 = transformed[leftStart * stateCount + rightStart + 1];
            final int base = leftStart * stateCount + rightStart;
            eigenBasisGradient[base]     += c[coeffBase]     * in0 + c[coeffBase + 1] * in1;
            eigenBasisGradient[base + 1] += c[coeffBase + 2] * in0 + c[coeffBase + 3] * in1;
        }

        // TWO_BY_ONE: 2-tall left, scalar right
        for (int e = plan.twoByOneStart; e < plan.twoByOneEnd; e++) {
            final int leftStart  = plan.eLeftStart[e];
            final int rightStart = plan.eRightStart[e];
            final int coeffBase  = e * 16;
            final double in0 = transformed[leftStart * stateCount + rightStart];
            final double in1 = transformed[(leftStart + 1) * stateCount + rightStart];
            final int base = leftStart * stateCount + rightStart;
            eigenBasisGradient[base]              += c[coeffBase]     * in0 + c[coeffBase + 1] * in1;
            eigenBasisGradient[base + stateCount] += c[coeffBase + 2] * in0 + c[coeffBase + 3] * in1;
        }

        // TWO_BY_TWO: 2×2 block pair — 4×4 coefficient matrix applied to 4-element input
        for (int e = plan.twoByTwoStart; e < plan.twoByTwoEnd; e++) {
            final int leftStart  = plan.eLeftStart[e];
            final int rightStart = plan.eRightStart[e];
            final int coeffBase  = e * 16;
            final double in0 = transformed[leftStart * stateCount + rightStart];
            final double in1 = transformed[leftStart * stateCount + rightStart + 1];
            final double in2 = transformed[(leftStart + 1) * stateCount + rightStart];
            final double in3 = transformed[(leftStart + 1) * stateCount + rightStart + 1];
            final int base = leftStart * stateCount + rightStart;
            eigenBasisGradient[base]                 += c[coeffBase]      * in0 + c[coeffBase + 1]  * in1 + c[coeffBase + 2]  * in2 + c[coeffBase + 3]  * in3;
            eigenBasisGradient[base + 1]             += c[coeffBase + 4]  * in0 + c[coeffBase + 5]  * in1 + c[coeffBase + 6]  * in2 + c[coeffBase + 7]  * in3;
            eigenBasisGradient[base + stateCount]     += c[coeffBase + 8]  * in0 + c[coeffBase + 9]  * in1 + c[coeffBase + 10] * in2 + c[coeffBase + 11] * in3;
            eigenBasisGradient[base + stateCount + 1] += c[coeffBase + 12] * in0 + c[coeffBase + 13] * in1 + c[coeffBase + 14] * in2 + c[coeffBase + 15] * in3;
        }

        // GENERAL fallback (very rare)
        for (int e = plan.generalStart; e < plan.generalEnd; e++) {
            final int leftStart  = plan.eLeftStart[e];
            final int leftDim    = plan.eLeftDim[e];
            final int rightStart = plan.eRightStart[e];
            final int rightDim   = plan.eRightDim[e];
            final int size       = leftDim * rightDim;
            final int coeffBase  = e * 16;
            final double[] in    = workspace.kernelInput[size];
            for (int i = 0; i < leftDim; ++i) {
                final int rowOff = (leftStart + i) * stateCount;
                for (int j = 0; j < rightDim; ++j) {
                    in[i * rightDim + j] = transformed[rowOff + rightStart + j];
                }
            }
            for (int row = 0; row < size; row++) {
                double sum = 0.0;
                final int rowBase = coeffBase + row * size;
                for (int col = 0; col < size; col++) {
                    sum += c[rowBase + col] * in[col];
                }
                eigenBasisGradient[(leftStart + row / rightDim) * stateCount + rightStart + row % rightDim] += sum;
            }
        }
    }

    public static void applyPlanToOuterProduct(ComplexKernelPlan plan,
                                               double[] leftVector,
                                               double[] rightVector,
                                               double scale,
                                               double[] eigenBasisGradient,
                                               int stateCount) {
        final double[] c = plan.coefficients;

        // SCALAR: each entry is a single outer-product dot with one coefficient
        for (int e = plan.scalarStart; e < plan.scalarEnd; e++) {
            final int leftStart  = plan.eLeftStart[e];
            final int rightStart = plan.eRightStart[e];
            eigenBasisGradient[leftStart * stateCount + rightStart] +=
                    c[e * 16] * scale * leftVector[leftStart] * rightVector[rightStart];
        }

        // ONE_BY_TWO: scalar left, 2-wide right — 2-element matvec
        for (int e = plan.oneByTwoStart; e < plan.oneByTwoEnd; e++) {
            final int leftStart  = plan.eLeftStart[e];
            final int rightStart = plan.eRightStart[e];
            final int coeffBase  = e * 16;
            final double x   = scale * leftVector[leftStart];
            final double in0 = x * rightVector[rightStart];
            final double in1 = x * rightVector[rightStart + 1];
            final int base   = leftStart * stateCount + rightStart;
            eigenBasisGradient[base]     += c[coeffBase]     * in0 + c[coeffBase + 1] * in1;
            eigenBasisGradient[base + 1] += c[coeffBase + 2] * in0 + c[coeffBase + 3] * in1;
        }

        // TWO_BY_ONE: 2-tall left, scalar right — 2-element matvec
        for (int e = plan.twoByOneStart; e < plan.twoByOneEnd; e++) {
            final int leftStart  = plan.eLeftStart[e];
            final int rightStart = plan.eRightStart[e];
            final int coeffBase  = e * 16;
            final double y   = rightVector[rightStart];
            final double in0 = scale * leftVector[leftStart]     * y;
            final double in1 = scale * leftVector[leftStart + 1] * y;
            final int base   = leftStart * stateCount + rightStart;
            eigenBasisGradient[base]              += c[coeffBase]     * in0 + c[coeffBase + 1] * in1;
            eigenBasisGradient[base + stateCount] += c[coeffBase + 2] * in0 + c[coeffBase + 3] * in1;
        }

        // TWO_BY_TWO: 2×2 outer product — 4×4 coefficient matrix applied to rank-1 input
        for (int e = plan.twoByTwoStart; e < plan.twoByTwoEnd; e++) {
            final int leftStart  = plan.eLeftStart[e];
            final int rightStart = plan.eRightStart[e];
            final int coeffBase  = e * 16;
            final double x0 = scale * leftVector[leftStart];
            final double x1 = scale * leftVector[leftStart + 1];
            final double y0 = rightVector[rightStart];
            final double y1 = rightVector[rightStart + 1];
            final double in0 = x0 * y0;
            final double in1 = x0 * y1;
            final double in2 = x1 * y0;
            final double in3 = x1 * y1;
            final int base = leftStart * stateCount + rightStart;
            eigenBasisGradient[base]                  += c[coeffBase]      * in0 + c[coeffBase + 1]  * in1 + c[coeffBase + 2]  * in2 + c[coeffBase + 3]  * in3;
            eigenBasisGradient[base + 1]              += c[coeffBase + 4]  * in0 + c[coeffBase + 5]  * in1 + c[coeffBase + 6]  * in2 + c[coeffBase + 7]  * in3;
            eigenBasisGradient[base + stateCount]     += c[coeffBase + 8]  * in0 + c[coeffBase + 9]  * in1 + c[coeffBase + 10] * in2 + c[coeffBase + 11] * in3;
            eigenBasisGradient[base + stateCount + 1] += c[coeffBase + 12] * in0 + c[coeffBase + 13] * in1 + c[coeffBase + 14] * in2 + c[coeffBase + 15] * in3;
        }

        // GENERAL fallback (very rare — block dimension > 2)
        for (int e = plan.generalStart; e < plan.generalEnd; e++) {
            applyGeneralOuterProductBlock(plan, e, leftVector, rightVector, scale, eigenBasisGradient, stateCount);
        }
    }

    /**
     * Single-pass fused fill+apply for one (leftVector ⊗ rightVector) pattern.
     * Computes each entry's coefficients on-the-fly and immediately accumulates
     * into eigenBasisGradient, avoiding the K²×16 intermediate coefficient store/load.
     * For patternCount=1 (phylogeography) this eliminates ~100 MB of memory traffic
     * per gradient evaluation at K=26/598 branches.
     */
    public static void fillAndApplyToOuterProduct(ComplexKernelPlan plan,
                                                   EigenDecomposition decomposition,
                                                   double time,
                                                   double[] leftVector,
                                                   double[] rightVector,
                                                   double scale,
                                                   double[] eigenBasisGradient,
                                                   int stateCount) {
        final double[] eigenValues = decomposition.getEigenValues();

        // Pre-compute exp/cos/sin per block — identical to fillTimeDependentCoefficients header.
        for (int b = 0; b < plan.blockCount; ++b) {
            final int start = plan.blockStarts[b];
            plan.expR[b] = FastMath.exp(time * eigenValues[start]);
            if (plan.blockDims[b] == 2) {
                final double imag = eigenValues[start + stateCount];
                plan.cosB[b] = Math.cos(time * imag);
                plan.sinB[b] = Math.sin(time * imag);
            }
        }

        // SCALAR: coefficient is a single scalar, used immediately.
        for (int e = plan.scalarStart; e < plan.scalarEnd; e++) {
            final int leftStart  = plan.eLeftStart[e];
            final int rightStart = plan.eRightStart[e];
            final double expLeft  = plan.expR[plan.eLeftBlock[e]];
            final double expRight = plan.expR[plan.eRightBlock[e]];
            final double invDenom  = plan.eInvDenom[e];
            final double coeff = invDenom == 0.0
                    ? time * expLeft
                    : (expLeft - expRight) * invDenom;
            eigenBasisGradient[leftStart * stateCount + rightStart] +=
                    coeff * scale * leftVector[leftStart] * rightVector[rightStart];
        }

        // ONE_BY_TWO: 2×2 rotation kernel [eic, eis; -eis, eic] applied to a rank-1 input.
        for (int e = plan.oneByTwoStart; e < plan.oneByTwoEnd; e++) {
            final int leftStart  = plan.eLeftStart[e];
            final int rightStart = plan.eRightStart[e];
            final int rb         = plan.eRightBlock[e];
            final double expLeft  = plan.expR[plan.eLeftBlock[e]];
            final double expShift = plan.expR[rb] / expLeft;
            final double cosRight = plan.cosB[rb];
            final double sinRight = plan.sinB[rb];
            final double realShift = plan.eRealShift[e];
            final double rightImag = plan.eRightImag[e];
            final double invDenom  = plan.eInvDenom[e];
            final double ic, is;
            if (invDenom == 0.0) {
                ic = time; is = 0.0;
            } else {
                ic = (expShift * (realShift * cosRight + rightImag * sinRight) - realShift) * invDenom;
                is = (expShift * (realShift * sinRight - rightImag * cosRight) + rightImag) * invDenom;
            }
            final double eic = expLeft * ic;
            final double eis = expLeft * is;
            final double x   = scale * leftVector[leftStart];
            final double in0 = x * rightVector[rightStart];
            final double in1 = x * rightVector[rightStart + 1];
            final int base = leftStart * stateCount + rightStart;
            eigenBasisGradient[base]     += eic * in0 + eis * in1;
            eigenBasisGradient[base + 1] += -eis * in0 + eic * in1;
        }

        // TWO_BY_ONE: 2×2 rotation kernel [A, B; -B, A] where A=cosLeft*eic+sinLeft*eis.
        for (int e = plan.twoByOneStart; e < plan.twoByOneEnd; e++) {
            final int leftStart  = plan.eLeftStart[e];
            final int rightStart = plan.eRightStart[e];
            final int lb         = plan.eLeftBlock[e];
            final double expLeft  = plan.expR[lb];
            final double expShift = plan.expR[plan.eRightBlock[e]] / expLeft;
            final double cosLeft  = plan.cosB[lb];
            final double sinLeft  = plan.sinB[lb];
            final double realShift = plan.eRealShift[e];
            final double leftImag  = plan.eLeftImag[e];
            final double invDenom  = plan.eInvDenom[e];
            final double ic, is;
            if (invDenom == 0.0) {
                ic = time; is = 0.0;
            } else {
                ic = (expShift * (realShift * cosLeft + leftImag * sinLeft) - realShift) * invDenom;
                is = (expShift * (realShift * sinLeft - leftImag * cosLeft) + leftImag) * invDenom;
            }
            final double eic = expLeft * ic;
            final double eis = expLeft * is;
            final double A   = cosLeft * eic + sinLeft * eis;
            final double B   = cosLeft * eis - sinLeft * eic;
            final double y   = rightVector[rightStart];
            final double in0 = scale * leftVector[leftStart]     * y;
            final double in1 = scale * leftVector[leftStart + 1] * y;
            final int base = leftStart * stateCount + rightStart;
            eigenBasisGradient[base]              +=  A * in0 + B * in1;
            eigenBasisGradient[base + stateCount] += -B * in0 + A * in1;
        }

        // TWO_BY_TWO: factored 4×4 matvec — only 8 multiplications for the apply step
        // (vs 16 in the naive matvec), using sum/difference identities on the outer-product inputs.
        for (int e = plan.twoByTwoStart; e < plan.twoByTwoEnd; e++) {
            final int leftStart  = plan.eLeftStart[e];
            final int rightStart = plan.eRightStart[e];
            final int lb = plan.eLeftBlock[e];
            final int rb = plan.eRightBlock[e];
            final double expLeft  = plan.expR[lb];
            final double expShift = plan.expR[rb] / expLeft;
            final double cosLeft  = plan.cosB[lb];
            final double sinLeft  = plan.sinB[lb];
            final double cosRight = plan.cosB[rb];
            final double sinRight = plan.sinB[rb];
            final double realShift      = plan.eRealShift[e];
            final double plusImagShift  = plan.ePlusImagShift[e];
            final double minusImagShift = plan.eMinusImagShift[e];
            final double invPlusDenom   = plan.eInvPlusDenom[e];
            final double invMinusDenom  = plan.eInvMinusDenom[e];

            // Cross-block angle-addition identities
            final double cosPlusImag  = cosLeft * cosRight - sinLeft * sinRight;
            final double sinPlusImag  = sinLeft * cosRight + cosLeft * sinRight;
            final double cosMinusImag = cosRight * cosLeft + sinRight * sinLeft;
            final double sinMinusImag = sinRight * cosLeft - cosRight * sinLeft;

            final double pi0, pi1;
            if (invPlusDenom == 0.0) {
                pi0 = time; pi1 = 0.0;
            } else {
                final double expCosPl = expShift * cosPlusImag;
                final double expSinPl = expShift * sinPlusImag;
                pi0 = (realShift * (expCosPl - 1.0) + plusImagShift  * expSinPl) * invPlusDenom;
                pi1 = (realShift * expSinPl - plusImagShift  * (expCosPl - 1.0)) * invPlusDenom;
            }

            final double mi0, mi1;
            if (invMinusDenom == 0.0) {
                mi0 = time; mi1 = 0.0;
            } else {
                final double expCosMi = expShift * cosMinusImag;
                final double expSinMi = expShift * sinMinusImag;
                mi0 = (realShift * (expCosMi - 1.0) + minusImagShift * expSinMi) * invMinusDenom;
                mi1 = (realShift * expSinMi - minusImagShift * (expCosMi - 1.0)) * invMinusDenom;
            }

            // expNeg / expPos complex exponentials for the left block
            final double enR =  expLeft * cosLeft;
            final double enI = -expLeft * sinLeft;
            final double epI =  expLeft * sinLeft;

            // plus = expNeg * plusIntegral,  minus = expPos * minusIntegral  (complex multiply)
            final double plR = enR * pi0 - enI * pi1;
            final double plI = enR * pi1 + enI * pi0;
            final double miR = enR * mi0 - epI * mi1;
            final double miI = enR * mi1 + epI * mi0;

            // Outer-product inputs (rank-1 2×2 matrix, stored as 4 scalars)
            final double x0 = scale * leftVector[leftStart];
            final double x1 = scale * leftVector[leftStart + 1];
            final double y0 = rightVector[rightStart];
            final double y1 = rightVector[rightStart + 1];

            // Sum/difference factoring: reduces 4×4 matvec from 16 to 8 multiplications.
            // in0=x0y0, in1=x0y1, in2=x1y0, in3=x1y1
            final double sPlusV  = x0 * y0 + x1 * y1;   // in0 + in3
            final double sMinusV = x0 * y0 - x1 * y1;   // in0 - in3
            final double tPlusV  = x0 * y1 + x1 * y0;   // in1 + in2
            final double tMinusV = x0 * y1 - x1 * y0;   // in1 - in2

            final double A = miR * sPlusV  + miI * tMinusV;
            final double B = plR * sMinusV + plI * tPlusV;
            final double C = miR * tMinusV - miI * sPlusV;
            final double D = plR * tPlusV  - plI * sMinusV;

            final int base = leftStart * stateCount + rightStart;
            eigenBasisGradient[base]                  += (A + B) * 0.5;
            eigenBasisGradient[base + 1]              += (C + D) * 0.5;
            eigenBasisGradient[base + stateCount]     += (D - C) * 0.5;
            eigenBasisGradient[base + stateCount + 1] += (A - B) * 0.5;
        }

        // GENERAL fallback (very rare): fill coefficient slots then apply.
        if (plan.generalStart < plan.generalEnd) {
            for (int e = plan.generalStart; e < plan.generalEnd; e++) {
                fillGeneralCoefficients(plan.coefficients, e * 16, eigenValues,
                        plan.eLeftStart[e], plan.eLeftDim[e],
                        plan.eRightStart[e], plan.eRightDim[e],
                        stateCount, time);
            }
            for (int e = plan.generalStart; e < plan.generalEnd; e++) {
                applyGeneralOuterProductBlock(plan, e, leftVector, rightVector,
                        scale, eigenBasisGradient, stateCount);
            }
        }
    }

    private static void applyGeneralOuterProductBlock(ComplexKernelPlan plan,
                                                      int entryIndex,
                                                      double[] leftVector,
                                                      double[] rightVector,
                                                      double scale,
                                                      double[] eigenBasisGradient,
                                                      int stateCount) {
        final int leftStart = plan.eLeftStart[entryIndex];
        final int rightStart = plan.eRightStart[entryIndex];
        final int leftDim = plan.eLeftDim[entryIndex];
        final int rightDim = plan.eRightDim[entryIndex];
        final int size = leftDim * rightDim;
        final int coeffBase = entryIndex * 16;

        for (int out = 0; out < size; out++) {
            double sum = 0.0;
            for (int inLeft = 0; inLeft < leftDim; inLeft++) {
                final double x = scale * leftVector[leftStart + inLeft];
                for (int inRight = 0; inRight < rightDim; inRight++) {
                    final int in = inLeft * rightDim + inRight;
                    sum += plan.coefficients[coeffBase + out * size + in] * x * rightVector[rightStart + inRight];
                }
            }
            final int outLeft = out / rightDim;
            final int outRight = out - outLeft * rightDim;
            eigenBasisGradient[(leftStart + outLeft) * stateCount + rightStart + outRight] += sum;
        }
    }

    public static final class Workspace {
        final double[][] kernelInput;

        public Workspace() {
            this.kernelInput = new double[5][];
            this.kernelInput[1] = new double[1];
            this.kernelInput[2] = new double[2];
            this.kernelInput[4] = new double[4];
        }
    }

    // -------------------------------------------------------------------------
    // Internal helpers
    // -------------------------------------------------------------------------

    private static int buildEigenBlocks(double[] eigenValues, int[] blockStarts, int[] blockDims, int stateCount) {
        int blockCount = 0;
        for (int i = 0; i < stateCount; ++i) {
            blockStarts[blockCount] = i;
            if (eigenValues.length > stateCount && Math.abs(eigenValues[stateCount + i]) > EIGEN_TOLERANCE) {
                blockDims[blockCount] = 2;
                i++;
            } else {
                blockDims[blockCount] = 1;
            }
            blockCount++;
        }
        return blockCount;
    }

    private static void fillEntryMetadata(ComplexKernelPlan plan,
                                          int entryIndex,
                                          double[] eigenValues,
                                          int stateCount,
                                          int leftStart,
                                          int leftDim,
                                          int rightStart,
                                          int rightDim) {
        final double leftReal = eigenValues[leftStart];
        final double rightReal = eigenValues[rightStart];
        final double realShift = rightReal - leftReal;
        plan.eRealShift[entryIndex] = realShift;

        if (leftDim == 1 && rightDim == 1) {
            plan.eKind[entryIndex] = ENTRY_SCALAR;
            final double delta = leftReal - rightReal;
            plan.eInvDenom[entryIndex] = Math.abs(delta) < EIGEN_TOLERANCE ? 0.0 : 1.0 / delta;
            return;
        }

        if (leftDim == 1 && rightDim == 2) {
            plan.eKind[entryIndex] = ENTRY_ONE_BY_TWO;
            final double rightImag = eigenValues[rightStart + stateCount];
            plan.eRightImag[entryIndex] = rightImag;
            final double denom = realShift * realShift + rightImag * rightImag;
            plan.eInvDenom[entryIndex] = denom < EIGEN_TOLERANCE ? 0.0 : 1.0 / denom;
            return;
        }

        if (leftDim == 2 && rightDim == 1) {
            plan.eKind[entryIndex] = ENTRY_TWO_BY_ONE;
            final double leftImag = eigenValues[leftStart + stateCount];
            plan.eLeftImag[entryIndex] = leftImag;
            final double denom = realShift * realShift + leftImag * leftImag;
            plan.eInvDenom[entryIndex] = denom < EIGEN_TOLERANCE ? 0.0 : 1.0 / denom;
            return;
        }

        if (leftDim == 2 && rightDim == 2) {
            plan.eKind[entryIndex] = ENTRY_TWO_BY_TWO;
            final double leftImag = eigenValues[leftStart + stateCount];
            final double rightImag = eigenValues[rightStart + stateCount];
            final double plusImagShift = leftImag + rightImag;
            final double minusImagShift = rightImag - leftImag;
            final double plusDenom = realShift * realShift + plusImagShift * plusImagShift;
            final double minusDenom = realShift * realShift + minusImagShift * minusImagShift;
            plan.eLeftImag[entryIndex] = leftImag;
            plan.eRightImag[entryIndex] = rightImag;
            plan.ePlusImagShift[entryIndex] = plusImagShift;
            plan.eMinusImagShift[entryIndex] = minusImagShift;
            plan.eInvPlusDenom[entryIndex] = plusDenom < EIGEN_TOLERANCE ? 0.0 : 1.0 / plusDenom;
            plan.eInvMinusDenom[entryIndex] = minusDenom < EIGEN_TOLERANCE ? 0.0 : 1.0 / minusDenom;
            return;
        }

        plan.eKind[entryIndex] = ENTRY_GENERAL;
    }

    // 1×1 scalar block pair. Uses pre-computed expLeft, expRight.
    private static void fillScalarCoefficient(double[] c, int off,
                                              double expLeft, double expRight,
                                              double inverseDelta,
                                              double time) {
        c[off] = inverseDelta == 0.0
                ? time * expLeft
                : (expLeft - expRight) * inverseDelta;
    }

    // 1×2 block pair (scalar left, complex right).
    // expLeft  = exp(t*leftEigenvalue)
    // expRight = exp(t*rightReal),  cosRight = cos(t*rightImag),  sinRight = sin(t*rightImag)
    private static void fillOneByTwoCoefficients(double[] c, int off,
                                                  double expLeft, double expRight,
                                                  double cosRight, double sinRight,
                                                  double realShift, double rightImag,
                                                  double inverseDenom,
                                                  double time) {
        final double expShift  = expRight / expLeft;  // exp(t * realShift)
        // Inline fillScaledRotationIntegralPrecomputed(realShift, rightImag, time, expShift, cosRight, sinRight)
        final double ic, is;
        if (inverseDenom == 0.0) {
            ic = time;
            is = 0.0;
        } else {
            ic = (expShift * (realShift * cosRight + rightImag * sinRight) - realShift) * inverseDenom;
            is = (expShift * (realShift * sinRight - rightImag * cosRight) + rightImag) * inverseDenom;
        }
        // Multiply by expLeft and write rotation matrix [ic, is; -is, ic] into c[off..off+3]
        c[off]     =  expLeft * ic;
        c[off + 1] =  expLeft * is;
        c[off + 2] = -expLeft * is;
        c[off + 3] =  expLeft * ic;
    }

    // 2×1 block pair (complex left, scalar right).
    // expLeft  = exp(t*leftReal),   cosLeft = cos(t*leftImag),   sinLeft = sin(t*leftImag)
    // expRight = exp(t*rightEigenvalue)
    private static void fillTwoByOneCoefficients(double[] c, int off,
                                                  double expLeft, double expRight,
                                                  double cosLeft, double sinLeft,
                                                  double realShift, double leftImag,
                                                  double inverseDenom,
                                                  double time) {
        final double expShift  = expRight / expLeft;  // exp(t * realShift)
        // Inline fillScaledRotationIntegralPrecomputed(realShift, leftImag, time, expShift, cosLeft, sinLeft)
        final double ic, is;
        if (inverseDenom == 0.0) {
            ic = time;
            is = 0.0;
        } else {
            ic = (expShift * (realShift * cosLeft + leftImag * sinLeft) - realShift) * inverseDenom;
            is = (expShift * (realShift * sinLeft - leftImag * cosLeft) + leftImag) * inverseDenom;
        }
        // leftExp matrix: L = exp(t*leftReal) * [cos  -sin ; sin  cos]
        final double l00 =  expLeft * cosLeft;
        final double l01 = -expLeft * sinLeft;
        final double l10 =  expLeft * sinLeft;
        final double l11 =  expLeft * cosLeft;
        // out = L * [ic, is; -is, ic]  (2×2 matrix stored row-major)
        c[off]     = l00 * ic + l01 * (-is);
        c[off + 1] = l00 * is + l01 * ic;
        c[off + 2] = l10 * ic + l11 * (-is);
        c[off + 3] = l10 * is + l11 * ic;
    }

    // 2×2 block pair (complex left, complex right).
    // All exp/trig values are pre-computed; angle-addition identities used for cross-block trig.
    // Fully inlined with scalar temporaries — zero local array allocations.
    private static void fillTwoByTwoCoefficients(double[] c, int off,
                                                  double expLeft, double expRight,
                                                  double cosLeft, double sinLeft,
                                                  double cosRight, double sinRight,
                                                  double realShift,
                                                  double plusImagShift, double minusImagShift,
                                                  double inversePlusDenom, double inverseMinusDenom,
                                                  double time) {
        final double expShift  = expRight / expLeft;   // exp(t*(rightReal - leftReal))

        // cos/sin of t*(leftImag + rightImag) via angle addition
        final double cosPlusImag  = cosLeft * cosRight - sinLeft * sinRight;
        final double sinPlusImag  = sinLeft * cosRight + cosLeft * sinRight;
        // cos/sin of t*(rightImag - leftImag) via angle addition
        final double cosMinusImag = cosRight * cosLeft + sinRight * sinLeft;
        final double sinMinusImag = sinRight * cosLeft - cosRight * sinLeft;

        // plusIntegral = integralComplexPrecomputed(realShift, leftImag + rightImag, ...)
        final double pi0, pi1;
        if (inversePlusDenom == 0.0) {
            pi0 = time; pi1 = 0.0;
        } else {
            final double expCos = expShift * cosPlusImag;
            final double expSin = expShift * sinPlusImag;
            pi0 = (realShift * (expCos - 1.0) + plusImagShift * expSin) * inversePlusDenom;
            pi1 = (realShift * expSin - plusImagShift * (expCos - 1.0)) * inversePlusDenom;
        }

        // minusIntegral = integralComplexPrecomputed(realShift, rightImag - leftImag, ...)
        final double mi0, mi1;
        if (inverseMinusDenom == 0.0) {
            mi0 = time; mi1 = 0.0;
        } else {
            final double expCos = expShift * cosMinusImag;
            final double expSin = expShift * sinMinusImag;
            mi0 = (realShift * (expCos - 1.0) + minusImagShift * expSin) * inverseMinusDenom;
            mi1 = (realShift * expSin - minusImagShift * (expCos - 1.0)) * inverseMinusDenom;
        }

        // expNeg = expComplex(leftReal, -leftImag, t) = {eL*cL, -eL*sL}
        // expPos = expComplex(leftReal, +leftImag, t) = {eL*cL, +eL*sL}
        final double enR = expLeft * cosLeft;
        final double enI = -expLeft * sinLeft;
        // epR = enR (same real part), epI = -enI
        final double epI = expLeft * sinLeft;

        // plus = complexMatrix(multiplyComplex(expNeg, plusIntegral))
        // multiplyComplex({enR,enI}, {pi0,pi1}) = {enR*pi0 - enI*pi1, enR*pi1 + enI*pi0}
        final double plR = enR * pi0 - enI * pi1;
        final double plI = enR * pi1 + enI * pi0;
        // complexMatrix({plR,plI}) stored as [plR, plI, -plI, plR]

        // minus = complexMatrix(multiplyComplex(expPos, minusIntegral))
        // multiplyComplex({enR,epI}, {mi0,mi1}) = {enR*mi0 - epI*mi1, enR*mi1 + epI*mi0}
        final double miR = enR * mi0 - epI * mi1;
        final double miI = enR * mi1 + epI * mi0;
        // complexMatrix({miR,miI}) stored as [miR, miI, -miI, miR]

        // Apply the 4 basis columns (unrolled):
        //   col k: u=bC[k][0], v=bC[k][1], p=bC[k][2], q=bC[k][3]
        //   outU = miR*u + miI*v,  outV = -miI*u + miR*v
        //   outP = plR*p + plI*q,  outQ = -plI*p + plR*q
        //   c[k] = outU+outP,  c[4+k] = outV+outQ,  c[8+k] = -outV+outQ,  c[12+k] = outU-outP

        // col 0: u=0.5, v=0, p=0.5, q=0
        final double u0 =  miR * 0.5,  v0 = -miI * 0.5,  p0 =  plR * 0.5,  q0 = -plI * 0.5;
        c[off]      = u0 + p0;  c[off + 4]  = v0 + q0;  c[off + 8]  = -v0 + q0;  c[off + 12] = u0 - p0;

        // col 1: u=0, v=0.5, p=0, q=0.5
        final double u1 =  miI * 0.5,  v1 =  miR * 0.5,  p1 =  plI * 0.5,  q1 =  plR * 0.5;
        c[off + 1]  = u1 + p1;  c[off + 5]  = v1 + q1;  c[off + 9]  = -v1 + q1;  c[off + 13] = u1 - p1;

        // col 2: u=0, v=-0.5, p=0, q=0.5
        final double u2 = -miI * 0.5,  v2 = -miR * 0.5,  p2 =  plI * 0.5,  q2 =  plR * 0.5;
        c[off + 2]  = u2 + p2;  c[off + 6]  = v2 + q2;  c[off + 10] = -v2 + q2;  c[off + 14] = u2 - p2;

        // col 3: u=0.5, v=0, p=-0.5, q=0
        final double u3 =  miR * 0.5,  v3 = -miI * 0.5,  p3 = -plR * 0.5,  q3 =  plI * 0.5;
        c[off + 3]  = u3 + p3;  c[off + 7]  = v3 + q3;  c[off + 11] = -v3 + q3;  c[off + 15] = u3 - p3;
    }

    // General fallback via Padé-13 for block dimensions > 2 (extremely rare).
    private static void fillGeneralCoefficients(double[] coefficients, int offset,
                                                double[] eigenValues,
                                                int leftStart, int leftDim,
                                                int rightStart, int rightDim,
                                                int stateCount, double time) {
        final int size = leftDim * rightDim;
        double[][] aTranspose = new double[leftDim][leftDim];
        double[][] leftExp    = new double[leftDim][leftDim];
        double[][] b          = new double[rightDim][rightDim];
        fillTransposeEigenBlock(aTranspose, leftExp, leftStart, leftDim, eigenValues, stateCount, time);
        fillOriginalEigenBlock(b, rightStart, rightDim, eigenValues, stateCount);

        double[][] generator = new double[size][size];
        fillRowMajorKroneckerGenerator(generator, aTranspose, b, leftDim, rightDim);

        double[][] phi   = phiMatrix(generator, time);
        double[] scratch = new double[size * size];
        leftMultiplyByExp(scratch, leftExp, phi, leftDim, rightDim);
        System.arraycopy(scratch, 0, coefficients, offset, scratch.length);
    }

    private static void zeroSquare(double[][] matrix, int dim) {
        for (int i = 0; i < dim; ++i) {
            java.util.Arrays.fill(matrix[i], 0, dim, 0.0);
        }
    }

    private static void fillTransposeEigenBlock(double[][] transposeBlock,
                                                double[][] expTransposeBlock,
                                                int eigenStart,
                                                int blockDim,
                                                double[] eigenValues,
                                                int stateCount,
                                                double time) {
        final double real = eigenValues[eigenStart];
        transposeBlock[0][0] = real;
        if (blockDim == 1) {
            expTransposeBlock[0][0] = Math.exp(time * real);
            return;
        }

        final double imag = eigenValues[eigenStart + stateCount];
        transposeBlock[0][1] = -imag;
        transposeBlock[1][0] =  imag;
        transposeBlock[1][1] =  real;

        final double exp = Math.exp(time * real);
        final double cos = Math.cos(time * imag);
        final double sin = Math.sin(time * imag);
        expTransposeBlock[0][0] =  exp * cos;
        expTransposeBlock[0][1] = -exp * sin;
        expTransposeBlock[1][0] =  exp * sin;
        expTransposeBlock[1][1] =  exp * cos;
    }

    private static void fillOriginalEigenBlock(double[][] block,
                                               int eigenStart,
                                               int blockDim,
                                               double[] eigenValues,
                                               int stateCount) {
        final double real = eigenValues[eigenStart];
        block[0][0] = real;
        if (blockDim == 1) {
            return;
        }

        final double imag = eigenValues[eigenStart + stateCount];
        block[0][1] =  imag;
        block[1][0] = -imag;
        block[1][1] =  real;
    }

    private static void fillRowMajorKroneckerGenerator(double[][] generator,
                                                       double[][] aTranspose,
                                                       double[][] b,
                                                       int leftDim,
                                                       int rightDim) {
        zeroSquare(generator, generator.length);

        for (int outLeft = 0; outLeft < leftDim; ++outLeft) {
            for (int inLeft = 0; inLeft < leftDim; ++inLeft) {
                final double value = -aTranspose[outLeft][inLeft];
                if (value == 0.0) continue;
                final int outBase = outLeft * rightDim;
                final int inBase  = inLeft  * rightDim;
                for (int right = 0; right < rightDim; ++right) {
                    generator[outBase + right][inBase + right] += value;
                }
            }
        }

        for (int left = 0; left < leftDim; ++left) {
            final int base = left * rightDim;
            for (int outRight = 0; outRight < rightDim; ++outRight) {
                final int outIndex = base + outRight;
                for (int inRight = 0; inRight < rightDim; ++inRight) {
                    generator[outIndex][base + inRight] += b[outRight][inRight];
                }
            }
        }
    }

    private static double[][] phiMatrix(double[][] generator, double time) {
        final int size      = generator.length;
        final int blockSize = 2 * size;
        double[][] block = new double[blockSize][blockSize];
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                block[i][j] = time * generator[i][j];
            }
            block[i][size + i] = 1.0;
        }

        double[][] expBlock = expmPade13(block);
        double[][] phi = new double[size][size];
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                phi[i][j] = time * expBlock[i][size + j];
            }
        }
        return phi;
    }

    private static void leftMultiplyByExp(double[] coefficients,
                                          double[][] leftExp,
                                          double[][] phi,
                                          int leftDim,
                                          int rightDim) {
        final int size = leftDim * rightDim;
        if (leftDim == 1) {
            final double scale = leftExp[0][0];
            for (int i = 0; i < size; ++i) {
                final int off = i * size;
                for (int j = 0; j < size; ++j) {
                    coefficients[off + j] = scale * phi[i][j];
                }
            }
            return;
        }

        for (int outLeft = 0; outLeft < leftDim; ++outLeft) {
            final int outBase = outLeft * rightDim;
            for (int inLeft = 0; inLeft < leftDim; ++inLeft) {
                final double scale = leftExp[outLeft][inLeft];
                if (scale == 0.0) continue;
                final int midBase = inLeft * rightDim;
                for (int right = 0; right < rightDim; ++right) {
                    final int outRow     = outBase + right;
                    final int phiRow     = midBase + right;
                    final int coeffOffset = outRow * size;
                    for (int col = 0; col < size; ++col) {
                        coefficients[coeffOffset + col] += scale * phi[phiRow][col];
                    }
                }
            }
        }
    }

    private static double[][] expmPade13(double[][] matrix) {
        final int n = matrix.length;
        double[][] ident = identity(n);
        final double norm  = matrixOneNorm(matrix);
        final int scale = norm > 0.0 ? Math.max(0, (int) Math.ceil(log2(norm / THETA13))) : 0;
        double[][] a = scale(matrix, 1.0 / Math.pow(2.0, scale));

        double[][] a2 = multiply(a, a);
        double[][] a4 = multiply(a2, a2);
        double[][] a6 = multiply(a4, a2);

        double[][] uInner = add(
                add(multiply(a6, add(add(scale(a6, PADE13[13]), scale(a4, PADE13[11])), scale(a2, PADE13[9]))), scale(a6, PADE13[7])),
                add(scale(a4, PADE13[5]), add(scale(a2, PADE13[3]), scale(ident, PADE13[1])))
        );
        double[][] u = multiply(a, uInner);

        double[][] v = add(
                add(multiply(a6, add(add(scale(a6, PADE13[12]), scale(a4, PADE13[10])), scale(a2, PADE13[8]))), scale(a6, PADE13[6])),
                add(scale(a4, PADE13[4]), add(scale(a2, PADE13[2]), scale(ident, PADE13[0])))
        );

        double[][] result = solve(subtract(v, u), add(v, u));
        for (int i = 0; i < scale; ++i) {
            result = multiply(result, result);
        }
        return result;
    }

    private static double matrixOneNorm(double[][] matrix) {
        double max = 0.0;
        for (int j = 0; j < matrix[0].length; j++) {
            double sum = 0.0;
            for (double[] row : matrix) {
                sum += Math.abs(row[j]);
            }
            max = Math.max(max, sum);
        }
        return max;
    }

    private static double[][] identity(int n) {
        double[][] out = new double[n][n];
        for (int i = 0; i < n; ++i) {
            out[i][i] = 1.0;
        }
        return out;
    }

    private static double[][] scale(double[][] matrix, double factor) {
        final int n = matrix.length;
        final int m = matrix[0].length;
        double[][] out = new double[n][m];
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                out[i][j] = factor * matrix[i][j];
            }
        }
        return out;
    }

    private static double[][] multiply(double[][] a, double[][] b) {
        final int n     = a.length;
        final int m     = b[0].length;
        final int inner = b.length;
        double[][] out = new double[n][m];
        for (int i = 0; i < n; ++i) {
            for (int k = 0; k < inner; ++k) {
                final double aik = a[i][k];
                if (aik == 0.0) continue;
                for (int j = 0; j < m; ++j) {
                    out[i][j] += aik * b[k][j];
                }
            }
        }
        return out;
    }

    private static double[][] add(double[][] a, double[][] b) {
        double[][] out = new double[a.length][a[0].length];
        for (int i = 0; i < a.length; ++i) {
            for (int j = 0; j < a[i].length; ++j) {
                out[i][j] = a[i][j] + b[i][j];
            }
        }
        return out;
    }

    private static double[][] subtract(double[][] a, double[][] b) {
        double[][] out = new double[a.length][a[0].length];
        for (int i = 0; i < a.length; ++i) {
            for (int j = 0; j < a[i].length; ++j) {
                out[i][j] = a[i][j] - b[i][j];
            }
        }
        return out;
    }

    private static double[][] solve(double[][] a, double[][] b) {
        final int n = a.length;
        final int m = b[0].length;
        double[][] aug = new double[n][n + m];
        for (int i = 0; i < n; ++i) {
            System.arraycopy(a[i], 0, aug[i], 0, n);
            System.arraycopy(b[i], 0, aug[i], n, m);
        }

        for (int col = 0; col < n; ++col) {
            int pivot = col;
            double max = Math.abs(aug[col][col]);
            for (int row = col + 1; row < n; ++row) {
                final double value = Math.abs(aug[row][col]);
                if (value > max) { max = value; pivot = row; }
            }
            if (pivot != col) {
                double[] tmp = aug[col];
                aug[col]   = aug[pivot];
                aug[pivot] = tmp;
            }

            final double pivotValue = aug[col][col];
            for (int j = col; j < n + m; ++j) {
                aug[col][j] /= pivotValue;
            }

            for (int row = 0; row < n; ++row) {
                if (row == col) continue;
                final double factor = aug[row][col];
                if (factor == 0.0) continue;
                for (int j = col; j < n + m; ++j) {
                    aug[row][j] -= factor * aug[col][j];
                }
            }
        }

        double[][] out = new double[n][m];
        for (int i = 0; i < n; ++i) {
            System.arraycopy(aug[i], n, out[i], 0, m);
        }
        return out;
    }

    private static double log2(double x) {
        return Math.log(x) / Math.log(2.0);
    }
}

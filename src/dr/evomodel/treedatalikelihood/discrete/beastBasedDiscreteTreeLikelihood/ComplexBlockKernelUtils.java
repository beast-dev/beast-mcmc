package dr.evomodel.treedatalikelihood.discrete.beastBasedDiscreteTreeLikelihood;

import dr.evomodel.substmodel.EigenDecomposition;
import org.apache.commons.math.util.FastMath;

/**
 * Fréchet-integral kernel for substitution models whose Jordan blocks are 1×1 (real
 * eigenvalue) or 2×2 (complex-conjugate pair).  Block dimension > 2 is not supported
 * and will throw {@link IllegalStateException} at plan-fill time.
 */
public final class ComplexBlockKernelUtils {

    private static final double EIGEN_TOLERANCE = 1.0e-12;

    private static final byte ENTRY_SCALAR     = 1;
    private static final byte ENTRY_ONE_BY_TWO = 2;
    private static final byte ENTRY_TWO_BY_ONE = 3;
    private static final byte ENTRY_TWO_BY_TWO = 4;

    private ComplexBlockKernelUtils() { }

    // -------------------------------------------------------------------------
    // Public flat plan type — metadata storage pre-allocated at construction time;
    // coefficient storage is sized exactly after the block structure is known.
    // -------------------------------------------------------------------------
    public static final class ComplexKernelPlan {
        private final int stateCount;
        final int[]    blockStarts;
        final int[]    blockDims;
        final double[] expR;
        final double[] cosB;
        final double[] sinB;

        int blockCount;
        int entryCount;
        final int[]    eLeftStart;
        final int[]    eLeftBlock;
        final int[]    eRightStart;
        final int[]    eRightBlock;
        final double[] eRealShift;
        final double[] eLeftImag;
        final double[] eRightImag;
        final double[] eInvDenom;
        final double[] ePlusImagShift;
        final double[] eMinusImagShift;
        final double[] eInvPlusDenom;
        final double[] eInvMinusDenom;

        // Kind-specialized coefficient storage.  Offsets are dense within each
        // grouped entry range, avoiding the 16-double stride for smaller blocks.
        double[] scalarCoefficients;
        double[] oneByTwoCoefficients;
        double[] twoByOneCoefficients;
        double[] twoByTwoCoefficients;

        // Entry ranges grouped by kind — type-specific loops, no per-entry branch dispatch.
        int scalarStart,   scalarEnd;
        int oneByTwoStart, oneByTwoEnd;
        int twoByOneStart, twoByOneEnd;
        int twoByTwoStart, twoByTwoEnd;

        public ComplexKernelPlan(int stateCount) {
            if (stateCount <= 0) {
                throw new IllegalArgumentException("stateCount must be positive");
            }
            this.stateCount = stateCount;
            final int maxPairs = stateCount * stateCount;
            this.blockStarts  = new int[stateCount];
            this.blockDims    = new int[stateCount];
            this.expR         = new double[stateCount];
            this.cosB         = new double[stateCount];
            this.sinB         = new double[stateCount];
            this.eLeftStart   = new int[maxPairs];
            this.eLeftBlock   = new int[maxPairs];
            this.eRightStart  = new int[maxPairs];
            this.eRightBlock  = new int[maxPairs];
            this.eRealShift   = new double[maxPairs];
            this.eLeftImag    = new double[maxPairs];
            this.eRightImag   = new double[maxPairs];
            this.eInvDenom    = new double[maxPairs];
            this.ePlusImagShift  = new double[maxPairs];
            this.eMinusImagShift = new double[maxPairs];
            this.eInvPlusDenom   = new double[maxPairs];
            this.eInvMinusDenom  = new double[maxPairs];
        }
    }

    public static void fillStructure(ComplexKernelPlan plan,
                                     EigenDecomposition decomposition,
                                     int stateCount) {
        checkStateCount(plan, stateCount);
        final double[] eigenValues = decomposition.getEigenValues();
        final int blockCount = buildEigenBlocks(eigenValues, plan.blockStarts, plan.blockDims, stateCount);
        plan.blockCount = blockCount;

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
        plan.entryCount    = entryIndex;

        ensureCoefficientCapacity(plan);
    }

    private static void ensureCoefficientCapacity(ComplexKernelPlan plan) {
        final int scalarCount   = plan.scalarEnd   - plan.scalarStart;
        final int oneByTwoCount = plan.oneByTwoEnd - plan.oneByTwoStart;
        final int twoByOneCount = plan.twoByOneEnd - plan.twoByOneStart;
        final int twoByTwoCount = plan.twoByTwoEnd - plan.twoByTwoStart;

        if (plan.scalarCoefficients == null || plan.scalarCoefficients.length != scalarCount) {
            plan.scalarCoefficients = new double[scalarCount];
        }
        if (plan.oneByTwoCoefficients == null || plan.oneByTwoCoefficients.length != oneByTwoCount * 4) {
            plan.oneByTwoCoefficients = new double[oneByTwoCount * 4];
        }
        if (plan.twoByOneCoefficients == null || plan.twoByOneCoefficients.length != twoByOneCount * 4) {
            plan.twoByOneCoefficients = new double[twoByOneCount * 4];
        }
        if (plan.twoByTwoCoefficients == null || plan.twoByTwoCoefficients.length != twoByTwoCount * 16) {
            plan.twoByTwoCoefficients = new double[twoByTwoCount * 16];
        }
    }

    private static void checkStateCount(ComplexKernelPlan plan, int stateCount) {
        if (stateCount != plan.stateCount) {
            throw new IllegalArgumentException("stateCount does not match plan");
        }
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
                else throw new IllegalStateException("Block dimension > 2 is not supported");

                if (entryKind != kind) continue;

                plan.eLeftStart[entryIndex]  = leftStart;
                plan.eLeftBlock[entryIndex]  = leftBlock;
                plan.eRightStart[entryIndex] = rightStart;
                plan.eRightBlock[entryIndex] = rightBlock;
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
        checkStateCount(plan, stateCount);
        final double[] eigenValues = decomposition.getEigenValues();

        for (int b = 0; b < plan.blockCount; ++b) {
            final int start = plan.blockStarts[b];
            plan.expR[b] = FastMath.exp(time * eigenValues[start]);
            if (plan.blockDims[b] == 2) {
                final double imag = eigenValues[start + stateCount];
                plan.cosB[b] = Math.cos(time * imag);
                plan.sinB[b] = Math.sin(time * imag);
            }
        }

        for (int e = plan.scalarStart, k = 0; e < plan.scalarEnd; e++, k++) {
            fillScalarCoefficient(plan.scalarCoefficients, k,
                    plan.expR[plan.eLeftBlock[e]],
                    plan.expR[plan.eRightBlock[e]], plan.eInvDenom[e], time);
        }

        for (int e = plan.oneByTwoStart, k = 0; e < plan.oneByTwoEnd; e++, k += 4) {
            final int rb = plan.eRightBlock[e];
            fillOneByTwoCoefficients(plan.oneByTwoCoefficients, k,
                    plan.expR[plan.eLeftBlock[e]], plan.expR[rb],
                    plan.cosB[rb], plan.sinB[rb],
                    plan.eRealShift[e], plan.eRightImag[e], plan.eInvDenom[e], time);
        }

        for (int e = plan.twoByOneStart, k = 0; e < plan.twoByOneEnd; e++, k += 4) {
            final int lb = plan.eLeftBlock[e];
            fillTwoByOneCoefficients(plan.twoByOneCoefficients, k,
                    plan.expR[lb], plan.expR[plan.eRightBlock[e]],
                    plan.cosB[lb], plan.sinB[lb],
                    plan.eRealShift[e], plan.eLeftImag[e], plan.eInvDenom[e], time);
        }

        for (int e = plan.twoByTwoStart, k = 0; e < plan.twoByTwoEnd; e++, k += 16) {
            final int lb = plan.eLeftBlock[e];
            final int rb = plan.eRightBlock[e];
            fillTwoByTwoCoefficients(plan.twoByTwoCoefficients, k,
                    plan.expR[lb], plan.expR[rb],
                    plan.cosB[lb], plan.sinB[lb],
                    plan.cosB[rb], plan.sinB[rb],
                    plan.eRealShift[e],
                    plan.ePlusImagShift[e], plan.eMinusImagShift[e],
                    plan.eInvPlusDenom[e], plan.eInvMinusDenom[e], time);
        }
    }

    public static void applyPlanToOuterProduct(ComplexKernelPlan plan,
                                               double[] leftVector,
                                               double[] rightVector,
                                               double scale,
                                               double[] eigenBasisGradient,
                                               int stateCount) {
        checkStateCount(plan, stateCount);
        for (int e = plan.scalarStart, k = 0; e < plan.scalarEnd; e++, k++) {
            final int leftStart  = plan.eLeftStart[e];
            final int rightStart = plan.eRightStart[e];
            eigenBasisGradient[leftStart * stateCount + rightStart] +=
                    plan.scalarCoefficients[k] *
                            scale * leftVector[leftStart] * rightVector[rightStart];
        }

        final double[] c12 = plan.oneByTwoCoefficients;
        for (int e = plan.oneByTwoStart, k = 0; e < plan.oneByTwoEnd; e++, k += 4) {
            final int leftStart  = plan.eLeftStart[e];
            final int rightStart = plan.eRightStart[e];
            final double x   = scale * leftVector[leftStart];
            final double in0 = x * rightVector[rightStart];
            final double in1 = x * rightVector[rightStart + 1];
            final int base   = leftStart * stateCount + rightStart;
            eigenBasisGradient[base]     += c12[k]     * in0 + c12[k + 1] * in1;
            eigenBasisGradient[base + 1] += c12[k + 2] * in0 + c12[k + 3] * in1;
        }

        final double[] c21 = plan.twoByOneCoefficients;
        for (int e = plan.twoByOneStart, k = 0; e < plan.twoByOneEnd; e++, k += 4) {
            final int leftStart  = plan.eLeftStart[e];
            final int rightStart = plan.eRightStart[e];
            final double y   = rightVector[rightStart];
            final double in0 = scale * leftVector[leftStart]     * y;
            final double in1 = scale * leftVector[leftStart + 1] * y;
            final int base   = leftStart * stateCount + rightStart;
            eigenBasisGradient[base]              += c21[k]     * in0 + c21[k + 1] * in1;
            eigenBasisGradient[base + stateCount] += c21[k + 2] * in0 + c21[k + 3] * in1;
        }

        final double[] c22 = plan.twoByTwoCoefficients;
        for (int e = plan.twoByTwoStart, k = 0; e < plan.twoByTwoEnd; e++, k += 16) {
            final int leftStart  = plan.eLeftStart[e];
            final int rightStart = plan.eRightStart[e];
            final double x0 = scale * leftVector[leftStart];
            final double x1 = scale * leftVector[leftStart + 1];
            final double y0 = rightVector[rightStart];
            final double y1 = rightVector[rightStart + 1];
            final double in0 = x0 * y0;
            final double in1 = x0 * y1;
            final double in2 = x1 * y0;
            final double in3 = x1 * y1;
            final int base = leftStart * stateCount + rightStart;
            eigenBasisGradient[base]                  += c22[k]      * in0 + c22[k + 1]  * in1 + c22[k + 2]  * in2 + c22[k + 3]  * in3;
            eigenBasisGradient[base + 1]              += c22[k + 4]  * in0 + c22[k + 5]  * in1 + c22[k + 6]  * in2 + c22[k + 7]  * in3;
            eigenBasisGradient[base + stateCount]     += c22[k + 8]  * in0 + c22[k + 9]  * in1 + c22[k + 10] * in2 + c22[k + 11] * in3;
            eigenBasisGradient[base + stateCount + 1] += c22[k + 12] * in0 + c22[k + 13] * in1 + c22[k + 14] * in2 + c22[k + 15] * in3;
        }
    }

    /**
     * Single-pass fused fill+apply for one (leftVector ⊗ rightVector) pattern.
     * Computes each entry's coefficients on-the-fly and immediately accumulates
     * into eigenBasisGradient, avoiding the K²×16 intermediate coefficient store/load.
     */
    public static void fillAndApplyToOuterProduct(ComplexKernelPlan plan,
                                                   EigenDecomposition decomposition,
                                                   double time,
                                                   double[] leftVector,
                                                   double[] rightVector,
                                                   double scale,
                                                   double[] eigenBasisGradient,
                                                   int stateCount) {
        checkStateCount(plan, stateCount);
        final double[] eigenValues = decomposition.getEigenValues();

        for (int b = 0; b < plan.blockCount; ++b) {
            final int start = plan.blockStarts[b];
            plan.expR[b] = FastMath.exp(time * eigenValues[start]);
            if (plan.blockDims[b] == 2) {
                final double imag = eigenValues[start + stateCount];
                plan.cosB[b] = Math.cos(time * imag);
                plan.sinB[b] = Math.sin(time * imag);
            }
        }

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

            final double enR =  expLeft * cosLeft;
            final double enI = -expLeft * sinLeft;
            final double epI =  expLeft * sinLeft;

            final double plR = enR * pi0 - enI * pi1;
            final double plI = enR * pi1 + enI * pi0;
            final double miR = enR * mi0 - epI * mi1;
            final double miI = enR * mi1 + epI * mi0;

            final double x0 = scale * leftVector[leftStart];
            final double x1 = scale * leftVector[leftStart + 1];
            final double y0 = rightVector[rightStart];
            final double y1 = rightVector[rightStart + 1];

            final double sPlusV  = x0 * y0 + x1 * y1;
            final double sMinusV = x0 * y0 - x1 * y1;
            final double tPlusV  = x0 * y1 + x1 * y0;
            final double tMinusV = x0 * y1 - x1 * y0;

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
        final double leftReal  = eigenValues[leftStart];
        final double rightReal = eigenValues[rightStart];
        final double realShift = rightReal - leftReal;
        plan.eRealShift[entryIndex] = realShift;

        if (leftDim == 1 && rightDim == 1) {
            final double delta = leftReal - rightReal;
            plan.eInvDenom[entryIndex] = Math.abs(delta) < EIGEN_TOLERANCE ? 0.0 : 1.0 / delta;
            return;
        }

        if (leftDim == 1 && rightDim == 2) {
            final double rightImag = eigenValues[rightStart + stateCount];
            plan.eRightImag[entryIndex] = rightImag;
            final double denom = realShift * realShift + rightImag * rightImag;
            plan.eInvDenom[entryIndex] = denom < EIGEN_TOLERANCE ? 0.0 : 1.0 / denom;
            return;
        }

        if (leftDim == 2 && rightDim == 1) {
            final double leftImag = eigenValues[leftStart + stateCount];
            plan.eLeftImag[entryIndex] = leftImag;
            final double denom = realShift * realShift + leftImag * leftImag;
            plan.eInvDenom[entryIndex] = denom < EIGEN_TOLERANCE ? 0.0 : 1.0 / denom;
            return;
        }

        // leftDim == 2 && rightDim == 2
        final double leftImag  = eigenValues[leftStart  + stateCount];
        final double rightImag = eigenValues[rightStart + stateCount];
        final double plusImagShift  = leftImag + rightImag;
        final double minusImagShift = rightImag - leftImag;
        final double plusDenom  = realShift * realShift + plusImagShift  * plusImagShift;
        final double minusDenom = realShift * realShift + minusImagShift * minusImagShift;
        plan.eLeftImag[entryIndex]       = leftImag;
        plan.eRightImag[entryIndex]      = rightImag;
        plan.ePlusImagShift[entryIndex]  = plusImagShift;
        plan.eMinusImagShift[entryIndex] = minusImagShift;
        plan.eInvPlusDenom[entryIndex]   = plusDenom  < EIGEN_TOLERANCE ? 0.0 : 1.0 / plusDenom;
        plan.eInvMinusDenom[entryIndex]  = minusDenom < EIGEN_TOLERANCE ? 0.0 : 1.0 / minusDenom;
    }

    // 1×1 scalar block pair.
    private static void fillScalarCoefficient(double[] c, int off,
                                              double expLeft, double expRight,
                                              double inverseDelta, double time) {
        c[off] = inverseDelta == 0.0
                ? time * expLeft
                : (expLeft - expRight) * inverseDelta;
    }

    // 1×2 block pair (scalar left, complex right).
    private static void fillOneByTwoCoefficients(double[] c, int off,
                                                  double expLeft, double expRight,
                                                  double cosRight, double sinRight,
                                                  double realShift, double rightImag,
                                                  double inverseDenom, double time) {
        final double expShift = expRight / expLeft;
        final double ic, is;
        if (inverseDenom == 0.0) {
            ic = time; is = 0.0;
        } else {
            ic = (expShift * (realShift * cosRight + rightImag * sinRight) - realShift) * inverseDenom;
            is = (expShift * (realShift * sinRight - rightImag * cosRight) + rightImag) * inverseDenom;
        }
        c[off]     =  expLeft * ic;
        c[off + 1] =  expLeft * is;
        c[off + 2] = -expLeft * is;
        c[off + 3] =  expLeft * ic;
    }

    // 2×1 block pair (complex left, scalar right).
    private static void fillTwoByOneCoefficients(double[] c, int off,
                                                  double expLeft, double expRight,
                                                  double cosLeft, double sinLeft,
                                                  double realShift, double leftImag,
                                                  double inverseDenom, double time) {
        final double expShift = expRight / expLeft;
        final double ic, is;
        if (inverseDenom == 0.0) {
            ic = time; is = 0.0;
        } else {
            ic = (expShift * (realShift * cosLeft + leftImag * sinLeft) - realShift) * inverseDenom;
            is = (expShift * (realShift * sinLeft - leftImag * cosLeft) + leftImag) * inverseDenom;
        }
        final double l00 =  expLeft * cosLeft;
        final double l01 = -expLeft * sinLeft;
        final double l10 =  expLeft * sinLeft;
        final double l11 =  expLeft * cosLeft;
        c[off]     = l00 * ic + l01 * (-is);
        c[off + 1] = l00 * is + l01 * ic;
        c[off + 2] = l10 * ic + l11 * (-is);
        c[off + 3] = l10 * is + l11 * ic;
    }

    // 2×2 block pair (complex left, complex right).
    private static void fillTwoByTwoCoefficients(double[] c, int off,
                                                  double expLeft, double expRight,
                                                  double cosLeft, double sinLeft,
                                                  double cosRight, double sinRight,
                                                  double realShift,
                                                  double plusImagShift, double minusImagShift,
                                                  double inversePlusDenom, double inverseMinusDenom,
                                                  double time) {
        final double expShift = expRight / expLeft;

        final double cosPlusImag  = cosLeft * cosRight - sinLeft * sinRight;
        final double sinPlusImag  = sinLeft * cosRight + cosLeft * sinRight;
        final double cosMinusImag = cosRight * cosLeft + sinRight * sinLeft;
        final double sinMinusImag = sinRight * cosLeft - cosRight * sinLeft;

        final double pi0, pi1;
        if (inversePlusDenom == 0.0) {
            pi0 = time; pi1 = 0.0;
        } else {
            final double expCos = expShift * cosPlusImag;
            final double expSin = expShift * sinPlusImag;
            pi0 = (realShift * (expCos - 1.0) + plusImagShift * expSin) * inversePlusDenom;
            pi1 = (realShift * expSin - plusImagShift * (expCos - 1.0)) * inversePlusDenom;
        }

        final double mi0, mi1;
        if (inverseMinusDenom == 0.0) {
            mi0 = time; mi1 = 0.0;
        } else {
            final double expCos = expShift * cosMinusImag;
            final double expSin = expShift * sinMinusImag;
            mi0 = (realShift * (expCos - 1.0) + minusImagShift * expSin) * inverseMinusDenom;
            mi1 = (realShift * expSin - minusImagShift * (expCos - 1.0)) * inverseMinusDenom;
        }

        final double enR =  expLeft * cosLeft;
        final double enI = -expLeft * sinLeft;
        final double epI =  expLeft * sinLeft;

        final double plR = enR * pi0 - enI * pi1;
        final double plI = enR * pi1 + enI * pi0;
        final double miR = enR * mi0 - epI * mi1;
        final double miI = enR * mi1 + epI * mi0;

        final double u0 =  miR * 0.5,  v0 = -miI * 0.5,  p0 =  plR * 0.5,  q0 = -plI * 0.5;
        c[off]      = u0 + p0;  c[off + 4]  = v0 + q0;  c[off + 8]  = -v0 + q0;  c[off + 12] = u0 - p0;

        final double u1 =  miI * 0.5,  v1 =  miR * 0.5,  p1 =  plI * 0.5,  q1 =  plR * 0.5;
        c[off + 1]  = u1 + p1;  c[off + 5]  = v1 + q1;  c[off + 9]  = -v1 + q1;  c[off + 13] = u1 - p1;

        final double u2 = -miI * 0.5,  v2 = -miR * 0.5,  p2 =  plI * 0.5,  q2 =  plR * 0.5;
        c[off + 2]  = u2 + p2;  c[off + 6]  = v2 + q2;  c[off + 10] = -v2 + q2;  c[off + 14] = u2 - p2;

        final double u3 =  miR * 0.5,  v3 = -miI * 0.5,  p3 = -plR * 0.5,  q3 =  plI * 0.5;
        c[off + 3]  = u3 + p3;  c[off + 7]  = v3 + q3;  c[off + 11] = -v3 + q3;  c[off + 15] = u3 - p3;
    }
}

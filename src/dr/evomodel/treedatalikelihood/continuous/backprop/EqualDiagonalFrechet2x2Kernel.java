package dr.evomodel.treedatalikelihood.continuous.backprop;

import org.apache.commons.math.util.FastMath;

// Shared numeric utilities for equal-diagonal 2x2/2x2 Fréchet coefficient computation.
// All math is copied verbatim from BlockDiagonalFrechetExactPlan; that file is NOT modified.
final class EqualDiagonalFrechet2x2Kernel {

    private static final double EPS = 1.0e-12;
    private static final double SMALL_ROOT_SERIES_CUTOFF = 1.0e-3;
    private static final int SMALL_ROOT_SERIES_ORDER = 4;
    private static final double MOMENT_SERIES_CUTOFF = 1.0e-7;
    private static final double MOMENT_SERIES_CUTOFF_SQUARED =
            MOMENT_SERIES_CUTOFF * MOMENT_SERIES_CUTOFF;
    private static final int MOMENT_SERIES_TERMS = 30;
    private static final int MAX_SMALL_ROOT_DERIVATIVE_ORDER = 2 * SMALL_ROOT_SERIES_ORDER + 1;
    private static final int MAX_SMALL_ROOT_MOMENT_POWER = 2 * MAX_SMALL_ROOT_DERIVATIVE_ORDER;
    private static final int MAX_FACTORIAL_INDEX = 2 * SMALL_ROOT_SERIES_ORDER + 1;
    private static final double[] INVERSE_FACTORIALS = buildInverseFactorials(MAX_FACTORIAL_INDEX);
    private static final double[][] BINOMIAL = buildBinomialTable(MAX_SMALL_ROOT_DERIVATIVE_ORDER);

    private EqualDiagonalFrechet2x2Kernel() {
    }

    // -------------------------------------------------------------------------
    // Public types
    // -------------------------------------------------------------------------

    static final class PreparedCoefficients {
        double alpha;
        double beta;
        double gamma;
        double eta;
    }

    static final class Workspace {
        final PreparedCoefficients prepared = new PreparedCoefficients();
        final ComplexWorkspace complex = new ComplexWorkspace();
        final RealSeriesWorkspace realSeries = new RealSeriesWorkspace();
        final double[] shifted = new double[4];
        final double[] phi = new double[4];
    }

    // -------------------------------------------------------------------------
    // Main public API
    // -------------------------------------------------------------------------

    // Compute alpha/beta/gamma/eta into out.prepared and store them in workspace.prepared.
    // a = -t*leftDiagonal, b = -t*rightDiagonal (already scaled)
    // scaledLeftProduct = t^2 * leftUpper * leftLower, rAbs = |t| * sqrt(|leftUpper*leftLower|)
    static void fillPrepared(final double a,
                             final double b,
                             final double scaledLeftProduct,
                             final double scaledRightProduct,
                             final boolean leftProductNonNegative,
                             final boolean rightProductNonNegative,
                             final double rAbs,
                             final double qAbs,
                             final PreparedCoefficients out,
                             final Workspace workspace) {
        final double rRe = leftProductNonNegative ? rAbs : 0.0;
        final double rIm = leftProductNonNegative ? 0.0 : rAbs;
        final double qRe = rightProductNonNegative ? qAbs : 0.0;
        final double qIm = rightProductNonNegative ? 0.0 : qAbs;

        if (rAbs < SMALL_ROOT_SERIES_CUTOFF && qAbs < SMALL_ROOT_SERIES_CUTOFF) {
            fillCoefficientsBivariateSmallRootSeries(a, b, scaledLeftProduct, scaledRightProduct, out, workspace.realSeries);
        } else if (rAbs < SMALL_ROOT_SERIES_CUTOFF) {
            fillCoefficientsSmallLeftRootSeries(a, b, qRe, qIm, scaledLeftProduct, out, workspace.complex);
        } else if (qAbs < SMALL_ROOT_SERIES_CUTOFF) {
            fillCoefficientsSmallRightRootSeries(a, b, rRe, rIm, scaledRightProduct, out, workspace.complex);
        } else if (leftProductNonNegative && rightProductNonNegative) {
            fillCoefficientsDistinctReal(a, b, rRe, qRe, out);
        } else if (!leftProductNonNegative && rightProductNonNegative) {
            fillCoefficientsDistinctLeftImagRightReal(a, b, rIm, qRe, out, workspace.complex);
        } else if (leftProductNonNegative) {
            fillCoefficientsDistinctLeftRealRightImag(a, b, rRe, qIm, out, workspace.complex);
        } else {
            fillCoefficientsDistinctBothImag(a, b, rIm, qIm, out, workspace.complex);
        }
    }

    // Fill 16 coefficients at out[offset..offset+15] from pre-computed alpha/beta/gamma/eta.
    // The layout is the unrolled form of fillTwoByTwoCoefficientFromPrepared for all 4 basis vectors.
    static void fillCoefficientMatrix(final PreparedCoefficients prepared,
                                      final double leftUpper,
                                      final double leftLower,
                                      final double rightUpper,
                                      final double rightLower,
                                      final double[] out,
                                      final int offset) {
        final double alpha = prepared.alpha;
        final double beta  = prepared.beta;
        final double gamma = prepared.gamma;
        final double eta   = prepared.eta;

        // basis 0: [1,0,0,0]
        out[offset + 0]  = alpha;
        out[offset + 1]  = gamma * rightLower;
        out[offset + 2]  = beta  * leftUpper;
        out[offset + 3]  = eta   * leftUpper * rightLower;
        // basis 1: [0,1,0,0]
        out[offset + 4]  = gamma * rightUpper;
        out[offset + 5]  = alpha;
        out[offset + 6]  = eta   * leftUpper * rightUpper;
        out[offset + 7]  = beta  * leftUpper;
        // basis 2: [0,0,1,0]
        out[offset + 8]  = beta  * leftLower;
        out[offset + 9]  = eta   * leftLower * rightLower;
        out[offset + 10] = alpha;
        out[offset + 11] = gamma * rightLower;
        // basis 3: [0,0,0,1]
        out[offset + 12] = eta   * leftLower * rightUpper;
        out[offset + 13] = beta  * leftLower;
        out[offset + 14] = gamma * rightUpper;
        out[offset + 15] = alpha;
    }

    static double dividedDifferenceExp(final double left, final double right) {
        if (Math.abs(left - right) < EPS) {
            return FastMath.exp(left);
        }
        return (FastMath.exp(left) - FastMath.exp(right)) / (left - right);
    }

    // Coefficients for a 1x2 block pair where the left is a scalar.
    // shifted = [rightDiag-leftScalar, rightUpper, rightLower, rightDiag-leftScalar]
    // Note the phi index swap (phi[2]/phi[1]) matching the existing fillOneByTwoCoefficient.
    static void fillOneByTwoCoefficients(final double leftScalar,
                                         final double rightDiag,
                                         final double rightUpper,
                                         final double rightLower,
                                         final double[] out,
                                         final int offset,
                                         final Workspace workspace) {
        final double[] shifted = workspace.shifted;
        final double[] phi = workspace.phi;
        shifted[0] = rightDiag - leftScalar;
        shifted[1] = rightUpper;
        shifted[2] = rightLower;
        shifted[3] = rightDiag - leftScalar;  // equal-diagonal: both diagonals are the same
        fillPhi2x2(shifted, phi);
        final double scale = FastMath.exp(leftScalar);
        out[offset]     = scale * phi[0];
        out[offset + 1] = scale * phi[2];  // transposed to match fillOneByTwoCoefficient
        out[offset + 2] = scale * phi[1];
        out[offset + 3] = scale * phi[3];
    }

    // Coefficients for a 2x1 block pair where the right is a scalar.
    // shifted = [leftDiag-rightScalar, leftUpper, leftLower, leftDiag-rightScalar]
    static void fillTwoByOneCoefficients(final double leftDiag,
                                          final double leftUpper,
                                          final double leftLower,
                                          final double rightScalar,
                                          final double[] out,
                                          final int offset,
                                          final Workspace workspace) {
        final double[] shifted = workspace.shifted;
        final double[] phi = workspace.phi;
        shifted[0] = leftDiag - rightScalar;
        shifted[1] = leftUpper;
        shifted[2] = leftLower;
        shifted[3] = leftDiag - rightScalar;  // equal-diagonal: both diagonals are the same
        fillPhi2x2(shifted, phi);
        final double scale = FastMath.exp(rightScalar);
        out[offset]     = scale * phi[0];
        out[offset + 1] = scale * phi[1];
        out[offset + 2] = scale * phi[2];
        out[offset + 3] = scale * phi[3];
    }

    // -------------------------------------------------------------------------
    // Private inner types (equivalent to those in BlockDiagonalFrechetExactPlan)
    // -------------------------------------------------------------------------

    private static final class MutableComplex {
        double re;
        double im;
    }

    private static final class ComplexWorkspace {
        final MutableComplex c0 = new MutableComplex();
        final MutableComplex c1 = new MutableComplex();
        final MutableComplex c2 = new MutableComplex();
        final MutableComplex c3 = new MutableComplex();
        final MutableComplex c4 = new MutableComplex();
        final MutableComplex c5 = new MutableComplex();
        final MutableComplex c6 = new MutableComplex();
        final MutableComplex c7 = new MutableComplex();
        boolean momentCacheValid;
        double momentCacheCRe;
        double momentCacheCIm;
        int momentCacheMaxPower;
        final double[] momentRe = new double[MAX_SMALL_ROOT_MOMENT_POWER + 1];
        final double[] momentIm = new double[MAX_SMALL_ROOT_MOMENT_POWER + 1];
    }

    private static final class RealSeriesWorkspace {
        final double[] moments = new double[MAX_SMALL_ROOT_MOMENT_POWER + 1];
        final double[] derivatives =
                new double[(MAX_SMALL_ROOT_DERIVATIVE_ORDER + 1) * (MAX_SMALL_ROOT_DERIVATIVE_ORDER + 1)];
    }

    // -------------------------------------------------------------------------
    // Private helpers (copied verbatim from BlockDiagonalFrechetExactPlan)
    // -------------------------------------------------------------------------

    private static void fillCoefficientsBivariateSmallRootSeries(final double a,
                                                                  final double b,
                                                                  final double rSquared,
                                                                  final double qSquared,
                                                                  final PreparedCoefficients out,
                                                                  final RealSeriesWorkspace workspace) {
        fillRealMomentTable(b - a, workspace.moments);
        fillRealDerivativeTable(FastMath.exp(a), workspace.moments, workspace.derivatives);

        double alpha = 0.0;
        double beta = 0.0;
        double gamma = 0.0;
        double eta = 0.0;

        double rPower = 1.0;
        for (int i = 0; i <= SMALL_ROOT_SERIES_ORDER; ++i) {
            final double rEven = rPower * INVERSE_FACTORIALS[2 * i];
            final double rOdd = rPower * INVERSE_FACTORIALS[2 * i + 1];

            double qPower = 1.0;
            for (int j = 0; j <= SMALL_ROOT_SERIES_ORDER; ++j) {
                final double qEven = qPower * INVERSE_FACTORIALS[2 * j];
                final double qOdd = qPower * INVERSE_FACTORIALS[2 * j + 1];

                alpha += rEven * qEven * workspace.derivatives[derivativeIndex(2 * i, 2 * j)];
                beta += rOdd * qEven * workspace.derivatives[derivativeIndex(2 * i + 1, 2 * j)];
                gamma += rEven * qOdd * workspace.derivatives[derivativeIndex(2 * i, 2 * j + 1)];
                eta += rOdd * qOdd * workspace.derivatives[derivativeIndex(2 * i + 1, 2 * j + 1)];

                qPower *= qSquared;
            }
            rPower *= rSquared;
        }

        out.alpha = alpha;
        out.beta = beta;
        out.gamma = gamma;
        out.eta = eta;
    }

    private static void fillCoefficientsSmallLeftRootSeries(final double a,
                                                             final double b,
                                                             final double qRe,
                                                             final double qIm,
                                                             final double rSquared,
                                                             final PreparedCoefficients out,
                                                             final ComplexWorkspace workspace) {
        double alpha = 0.0;
        double beta = 0.0;
        double gamma = 0.0;
        double eta = 0.0;

        final MutableComplex plus = workspace.c0;
        final MutableComplex minus = workspace.c1;
        final MutableComplex scratch = workspace.c2;

        double rPower = 1.0;
        for (int i = 0; i <= SMALL_ROOT_SERIES_ORDER; ++i) {
            fillMixedDerivative(a, 0.0, b + qRe, qIm, 2 * i, 0, plus, workspace);
            fillMixedDerivative(a, 0.0, b - qRe, -qIm, 2 * i, 0, minus, workspace);
            final double evenAverageRe = 0.5 * (plus.re + minus.re);
            divideComplex(plus.re - minus.re, plus.im - minus.im, qRe, qIm, scratch);
            final double oddQEvenR = 0.5 * scratch.re;

            fillMixedDerivative(a, 0.0, b + qRe, qIm, 2 * i + 1, 0, plus, workspace);
            fillMixedDerivative(a, 0.0, b - qRe, -qIm, 2 * i + 1, 0, minus, workspace);
            final double oddAverageRe = 0.5 * (plus.re + minus.re);
            divideComplex(plus.re - minus.re, plus.im - minus.im, qRe, qIm, scratch);
            final double oddQOddR = 0.5 * scratch.re;

            alpha += rPower * INVERSE_FACTORIALS[2 * i] * evenAverageRe;
            beta += rPower * INVERSE_FACTORIALS[2 * i + 1] * oddAverageRe;
            gamma += rPower * INVERSE_FACTORIALS[2 * i] * oddQEvenR;
            eta += rPower * INVERSE_FACTORIALS[2 * i + 1] * oddQOddR;
            rPower *= rSquared;
        }

        out.alpha = alpha;
        out.beta = beta;
        out.gamma = gamma;
        out.eta = eta;
    }

    private static void fillCoefficientsSmallRightRootSeries(final double a,
                                                              final double b,
                                                              final double rRe,
                                                              final double rIm,
                                                              final double qSquared,
                                                              final PreparedCoefficients out,
                                                              final ComplexWorkspace workspace) {
        double alpha = 0.0;
        double beta = 0.0;
        double gamma = 0.0;
        double eta = 0.0;

        final MutableComplex plus = workspace.c0;
        final MutableComplex minus = workspace.c1;
        final MutableComplex scratch = workspace.c2;

        double qPower = 1.0;
        for (int j = 0; j <= SMALL_ROOT_SERIES_ORDER; ++j) {
            fillMixedDerivative(a + rRe, rIm, b, 0.0, 0, 2 * j, plus, workspace);
            fillMixedDerivative(a - rRe, -rIm, b, 0.0, 0, 2 * j, minus, workspace);
            final double evenAverageRe = 0.5 * (plus.re + minus.re);
            divideComplex(plus.re - minus.re, plus.im - minus.im, rRe, rIm, scratch);
            final double oddREvenQ = 0.5 * scratch.re;

            fillMixedDerivative(a + rRe, rIm, b, 0.0, 0, 2 * j + 1, plus, workspace);
            fillMixedDerivative(a - rRe, -rIm, b, 0.0, 0, 2 * j + 1, minus, workspace);
            final double oddAverageRe = 0.5 * (plus.re + minus.re);
            divideComplex(plus.re - minus.re, plus.im - minus.im, rRe, rIm, scratch);
            final double oddROddQ = 0.5 * scratch.re;

            alpha += qPower * INVERSE_FACTORIALS[2 * j] * evenAverageRe;
            beta += qPower * INVERSE_FACTORIALS[2 * j] * oddREvenQ;
            gamma += qPower * INVERSE_FACTORIALS[2 * j + 1] * oddAverageRe;
            eta += qPower * INVERSE_FACTORIALS[2 * j + 1] * oddROddQ;
            qPower *= qSquared;
        }

        out.alpha = alpha;
        out.beta = beta;
        out.gamma = gamma;
        out.eta = eta;
    }

    private static void fillCoefficientsDistinctReal(final double a,
                                                      final double b,
                                                      final double r,
                                                      final double q,
                                                      final PreparedCoefficients out) {
        // Compute each of the 4 unique exponentials once to avoid the 8 exp calls
        // that result from calling dividedDifferenceExp(a±r, b±q) independently.
        final double eAR = FastMath.exp(a + r);
        final double eAM = FastMath.exp(a - r);
        final double eBQ = FastMath.exp(b + q);
        final double eBM = FastMath.exp(b - q);

        final double dPP = (a + r) - (b + q);
        final double dPM = (a + r) - (b - q);
        final double dMP = (a - r) - (b + q);
        final double dMM = (a - r) - (b - q);

        final double phiPP = Math.abs(dPP) < EPS ? eAR : (eAR - eBQ) / dPP;
        final double phiPM = Math.abs(dPM) < EPS ? eAR : (eAR - eBM) / dPM;
        final double phiMP = Math.abs(dMP) < EPS ? eAM : (eAM - eBQ) / dMP;
        final double phiMM = Math.abs(dMM) < EPS ? eAM : (eAM - eBM) / dMM;

        out.alpha = 0.25 * (phiPP + phiPM + phiMP + phiMM);
        out.beta  = 0.25 * (phiPP + phiPM - phiMP - phiMM) / r;
        out.gamma = 0.25 * (phiPP - phiPM + phiMP - phiMM) / q;
        out.eta   = 0.25 * (phiPP - phiPM - phiMP + phiMM) / (r * q);
    }

    private static void fillCoefficientsDistinctLeftImagRightReal(final double a,
                                                                   final double b,
                                                                   final double r,
                                                                   final double q,
                                                                   final PreparedCoefficients out,
                                                                   final ComplexWorkspace workspace) {
        final MutableComplex phiPlus = workspace.c0;
        final MutableComplex phiMinus = workspace.c1;

        fillPhiIntegral(a, r, b + q, 0.0, phiPlus, workspace);
        fillPhiIntegral(a, r, b - q, 0.0, phiMinus, workspace);

        out.alpha = 0.5 * (phiPlus.re + phiMinus.re);
        out.beta = 0.5 * (phiPlus.im + phiMinus.im) / r;
        out.gamma = 0.5 * (phiPlus.re - phiMinus.re) / q;
        out.eta = 0.5 * (phiPlus.im - phiMinus.im) / (r * q);
    }

    private static void fillCoefficientsDistinctLeftRealRightImag(final double a,
                                                                   final double b,
                                                                   final double r,
                                                                   final double q,
                                                                   final PreparedCoefficients out,
                                                                   final ComplexWorkspace workspace) {
        final MutableComplex phiPlus = workspace.c0;
        final MutableComplex phiMinus = workspace.c1;

        fillPhiIntegral(a + r, 0.0, b, q, phiPlus, workspace);
        fillPhiIntegral(a - r, 0.0, b, q, phiMinus, workspace);

        out.alpha = 0.5 * (phiPlus.re + phiMinus.re);
        out.beta = 0.5 * (phiPlus.re - phiMinus.re) / r;
        out.gamma = 0.5 * (phiPlus.im + phiMinus.im) / q;
        out.eta = 0.5 * (phiPlus.im - phiMinus.im) / (r * q);
    }

    private static void fillCoefficientsDistinctBothImag(final double a,
                                                          final double b,
                                                          final double r,
                                                          final double q,
                                                          final PreparedCoefficients out,
                                                          final ComplexWorkspace workspace) {
        final MutableComplex phiSameSign = workspace.c0;
        final MutableComplex phiOppositeSign = workspace.c1;

        fillPhiIntegral(a, r, b, q, phiSameSign, workspace);
        fillPhiIntegral(a, r, b, -q, phiOppositeSign, workspace);

        out.alpha = 0.5 * (phiSameSign.re + phiOppositeSign.re);
        out.beta = 0.5 * (phiSameSign.im + phiOppositeSign.im) / r;
        out.gamma = 0.5 * (phiSameSign.im - phiOppositeSign.im) / q;
        out.eta = -0.5 * (phiSameSign.re - phiOppositeSign.re) / (r * q);
    }

    private static void fillCoefficientsDistinct(final double a,
                                                  final double b,
                                                  final double rRe,
                                                  final double rIm,
                                                  final double qRe,
                                                  final double qIm,
                                                  final PreparedCoefficients out,
                                                  final ComplexWorkspace workspace) {
        final MutableComplex phiPP = workspace.c0;
        final MutableComplex phiPM = workspace.c1;
        final MutableComplex phiMP = workspace.c2;
        final MutableComplex phiMM = workspace.c3;
        final MutableComplex scratch = workspace.c4;

        fillPhiIntegral(a + rRe, rIm, b + qRe, qIm, phiPP, workspace);
        fillPhiIntegral(a + rRe, rIm, b - qRe, -qIm, phiPM, workspace);
        fillPhiIntegral(a - rRe, -rIm, b + qRe, qIm, phiMP, workspace);
        fillPhiIntegral(a - rRe, -rIm, b - qRe, -qIm, phiMM, workspace);

        final double alpha = 0.25 * (phiPP.re + phiPM.re + phiMP.re + phiMM.re);

        divideComplex(
                phiPP.re + phiPM.re - phiMP.re - phiMM.re,
                phiPP.im + phiPM.im - phiMP.im - phiMM.im,
                rRe, rIm, scratch);
        final double beta = 0.25 * scratch.re;

        divideComplex(
                phiPP.re - phiPM.re + phiMP.re - phiMM.re,
                phiPP.im - phiPM.im + phiMP.im - phiMM.im,
                qRe, qIm, scratch);
        final double gamma = 0.25 * scratch.re;

        multiplyComplex(rRe, rIm, qRe, qIm, scratch);
        divideComplex(
                phiPP.re - phiPM.re - phiMP.re + phiMM.re,
                phiPP.im - phiPM.im - phiMP.im + phiMM.im,
                scratch.re, scratch.im, scratch);
        out.alpha = alpha;
        out.beta = beta;
        out.gamma = gamma;
        out.eta = 0.25 * scratch.re;
    }

    private static void fillPhiIntegral(final double lambdaRe,
                                         final double lambdaIm,
                                         final double muRe,
                                         final double muIm,
                                         final MutableComplex out,
                                         final ComplexWorkspace workspace) {
        fillIntegralMoment(muRe - lambdaRe, muIm - lambdaIm, 0, workspace.c5, workspace);
        final double scale = FastMath.exp(lambdaRe);
        final double cos = Math.cos(lambdaIm);
        final double sin = Math.sin(lambdaIm);
        out.re = scale * (cos * workspace.c5.re - sin * workspace.c5.im);
        out.im = scale * (sin * workspace.c5.re + cos * workspace.c5.im);
    }

    private static void fillIntegralMoment(final double cRe,
                                            final double cIm,
                                            final int power,
                                            final MutableComplex out,
                                            final ComplexWorkspace workspace) {
        if (cRe * cRe + cIm * cIm < MOMENT_SERIES_CUTOFF_SQUARED) {
            fillIntegrateMomentSeries(cRe, cIm, power, out);
            return;
        }

        final double scale = FastMath.exp(cRe);
        final double cos = Math.cos(cIm);
        final double sin = Math.sin(cIm);
        final double expRe = scale * cos;
        final double expIm = scale * sin;

        divideComplex(expRe - 1.0, expIm, cRe, cIm, out);
        if (power == 0) {
            return;
        }

        final MutableComplex expOverC = workspace.c6;
        final MutableComplex divided = workspace.c7;
        divideComplex(expRe, expIm, cRe, cIm, expOverC);
        for (int n = 1; n <= power; ++n) {
            divideComplex(n * out.re, n * out.im, cRe, cIm, divided);
            out.re = expOverC.re - divided.re;
            out.im = expOverC.im - divided.im;
        }
    }

    private static void fillIntegrateMomentSeries(final double cRe,
                                                   final double cIm,
                                                   final int power,
                                                   final MutableComplex out) {
        double totalRe = 0.0;
        double totalIm = 0.0;
        double cPowerRe = 1.0;
        double cPowerIm = 0.0;
        double factorial = 1.0;
        for (int n = 0; n < MOMENT_SERIES_TERMS; ++n) {
            if (n > 0) {
                factorial *= n;
                final double nextRe = cPowerRe * cRe - cPowerIm * cIm;
                final double nextIm = cPowerRe * cIm + cPowerIm * cRe;
                cPowerRe = nextRe;
                cPowerIm = nextIm;
            }
            final double termScale = 1.0 / (factorial * (n + power + 1.0));
            totalRe += termScale * cPowerRe;
            totalIm += termScale * cPowerIm;
        }
        out.re = totalRe;
        out.im = totalIm;
    }

    private static void fillComplexMomentTable(final double cRe,
                                                final double cIm,
                                                final int maxPower,
                                                final ComplexWorkspace workspace) {
        if (workspace.momentCacheValid
                && Double.doubleToLongBits(cRe) == Double.doubleToLongBits(workspace.momentCacheCRe)
                && Double.doubleToLongBits(cIm) == Double.doubleToLongBits(workspace.momentCacheCIm)
                && maxPower <= workspace.momentCacheMaxPower) {
            return;
        }

        if (cRe * cRe + cIm * cIm < MOMENT_SERIES_CUTOFF_SQUARED) {
            fillComplexMomentTableSeries(cRe, cIm, maxPower, workspace);
            markComplexMomentCache(cRe, cIm, maxPower, workspace);
            return;
        }

        final double scale = FastMath.exp(cRe);
        final double cos = Math.cos(cIm);
        final double sin = Math.sin(cIm);

        final double expRe = scale * cos;
        final double expIm = scale * sin;

        final double denom = cRe * cRe + cIm * cIm;

        workspace.momentRe[0] = ((expRe - 1.0) * cRe + expIm * cIm) / denom;
        workspace.momentIm[0] = (expIm * cRe - (expRe - 1.0) * cIm) / denom;

        if (maxPower > 0) {
            final double expOverCRe = (expRe * cRe + expIm * cIm) / denom;
            final double expOverCIm = (expIm * cRe - expRe * cIm) / denom;

            for (int power = 1; power <= maxPower; ++power) {
                final double previousRe = workspace.momentRe[power - 1];
                final double previousIm = workspace.momentIm[power - 1];

                final double numeratorRe = power * previousRe;
                final double numeratorIm = power * previousIm;

                final double dividedRe = (numeratorRe * cRe + numeratorIm * cIm) / denom;
                final double dividedIm = (numeratorIm * cRe - numeratorRe * cIm) / denom;

                workspace.momentRe[power] = expOverCRe - dividedRe;
                workspace.momentIm[power] = expOverCIm - dividedIm;
            }
        }

        markComplexMomentCache(cRe, cIm, maxPower, workspace);
    }

    private static void markComplexMomentCache(final double cRe,
                                                final double cIm,
                                                final int maxPower,
                                                final ComplexWorkspace workspace) {
        workspace.momentCacheValid = true;
        workspace.momentCacheCRe = cRe;
        workspace.momentCacheCIm = cIm;
        workspace.momentCacheMaxPower = maxPower;
    }

    private static void fillComplexMomentTableSeries(final double cRe,
                                                      final double cIm,
                                                      final int maxPower,
                                                      final ComplexWorkspace workspace) {
        for (int power = 0; power <= maxPower; ++power) {
            workspace.momentRe[power] = 0.0;
            workspace.momentIm[power] = 0.0;
        }

        double cPowerRe = 1.0;
        double cPowerIm = 0.0;
        double factorial = 1.0;

        for (int n = 0; n < MOMENT_SERIES_TERMS; ++n) {
            if (n > 0) {
                factorial *= n;

                final double nextRe = cPowerRe * cRe - cPowerIm * cIm;
                final double nextIm = cPowerRe * cIm + cPowerIm * cRe;
                cPowerRe = nextRe;
                cPowerIm = nextIm;
            }

            final double baseRe = cPowerRe / factorial;
            final double baseIm = cPowerIm / factorial;

            for (int power = 0; power <= maxPower; ++power) {
                final double denominator = n + power + 1.0;
                workspace.momentRe[power] += baseRe / denominator;
                workspace.momentIm[power] += baseIm / denominator;
            }
        }
    }

    private static void fillMixedDerivative(final double lambdaRe,
                                             final double lambdaIm,
                                             final double muRe,
                                             final double muIm,
                                             final int leftOrder,
                                             final int rightOrder,
                                             final MutableComplex out,
                                             final ComplexWorkspace workspace) {
        final double cRe = muRe - lambdaRe;
        final double cIm = muIm - lambdaIm;
        final int maxMomentPower = rightOrder + leftOrder;

        fillComplexMomentTable(cRe, cIm, maxMomentPower, workspace);

        double totalRe = 0.0;
        double totalIm = 0.0;

        for (int i = 0; i <= leftOrder; ++i) {
            final double coefficient =
                    ((i & 1) == 0 ? 1.0 : -1.0) * BINOMIAL[leftOrder][i];
            final int momentIndex = rightOrder + i;
            totalRe += coefficient * workspace.momentRe[momentIndex];
            totalIm += coefficient * workspace.momentIm[momentIndex];
        }

        final double scale = FastMath.exp(lambdaRe);
        final double cos = Math.cos(lambdaIm);
        final double sin = Math.sin(lambdaIm);

        out.re = scale * (cos * totalRe - sin * totalIm);
        out.im = scale * (sin * totalRe + cos * totalIm);
    }

    private static void fillRealDerivativeTable(final double expLambda,
                                                 final double[] moments,
                                                 final double[] derivatives) {
        for (int leftOrder = 0; leftOrder <= MAX_SMALL_ROOT_DERIVATIVE_ORDER; ++leftOrder) {
            final int derivativeOffset = leftOrder * (MAX_SMALL_ROOT_DERIVATIVE_ORDER + 1);
            for (int rightOrder = 0; rightOrder <= MAX_SMALL_ROOT_DERIVATIVE_ORDER; ++rightOrder) {
                double total = 0.0;
                for (int i = 0; i <= leftOrder; ++i) {
                    final double coefficient = ((i & 1) == 0 ? 1.0 : -1.0) * BINOMIAL[leftOrder][i];
                    total += coefficient * moments[rightOrder + i];
                }
                derivatives[derivativeOffset + rightOrder] = expLambda * total;
            }
        }
    }

    private static void fillRealMomentTable(final double c,
                                             final double[] moments) {
        if (Math.abs(c) < MOMENT_SERIES_CUTOFF) {
            fillRealMomentTableSeries(c, moments);
            return;
        }

        final double expC = FastMath.exp(c);
        final double expOverC = expC / c;
        moments[0] = (expC - 1.0) / c;
        for (int power = 1; power <= MAX_SMALL_ROOT_MOMENT_POWER; ++power) {
            moments[power] = expOverC - (power * moments[power - 1]) / c;
        }
    }

    private static void fillRealMomentTableSeries(final double c,
                                                   final double[] moments) {
        for (int power = 0; power <= MAX_SMALL_ROOT_MOMENT_POWER; ++power) {
            moments[power] = 0.0;
        }

        double cPower = 1.0;
        double factorial = 1.0;
        for (int n = 0; n < MOMENT_SERIES_TERMS; ++n) {
            if (n > 0) {
                factorial *= n;
                cPower *= c;
            }
            final double base = cPower / factorial;
            for (int power = 0; power <= MAX_SMALL_ROOT_MOMENT_POWER; ++power) {
                moments[power] += base / (n + power + 1.0);
            }
        }
    }

    private static int derivativeIndex(final int leftOrder,
                                       final int rightOrder) {
        return leftOrder * (MAX_SMALL_ROOT_DERIVATIVE_ORDER + 1) + rightOrder;
    }

    private static void fillPhi2x2(final double[] matrix,
                                    final double[] out) {
        final double a = matrix[0];
        final double b = matrix[1];
        final double c = matrix[2];
        final double d = matrix[3];

        final double tau = 0.5 * (a + d);
        final double diff = 0.5 * (a - d);
        final double delta2 = diff * diff + b * c;

        final double i0;
        final double i1;

        if (Math.abs(delta2) < EPS) {
            i0 = scalarExpIntegral(tau);
            i1 = weightedExpIntegral(tau);
        } else if (delta2 > 0.0) {
            final double delta = Math.sqrt(delta2);
            final double plus = scalarExpIntegral(tau + delta);
            final double minus = scalarExpIntegral(tau - delta);
            i0 = 0.5 * (plus + minus);
            i1 = 0.5 * (plus - minus) / delta;
        } else {
            final double omega = Math.sqrt(-delta2);
            final double denom = tau * tau + omega * omega;
            final double expTau = FastMath.exp(tau);
            final double cos = Math.cos(omega);
            final double sin = Math.sin(omega);
            i0 = (expTau * (tau * cos + omega * sin) - tau) / denom;
            i1 = (expTau * (tau * sin - omega * cos) + omega) / (omega * denom);
        }

        final double g1 = i1;
        final double g0 = i0 - tau * g1;

        out[0] = g0 + g1 * a;
        out[1] = g1 * b;
        out[2] = g1 * c;
        out[3] = g0 + g1 * d;
    }

    private static double scalarExpIntegral(final double value) {
        if (Math.abs(value) < EPS) {
            return 1.0;
        }
        return Math.expm1(value) / value;
    }

    private static double weightedExpIntegral(final double value) {
        if (Math.abs(value) < 1.0e-8) {
            return 0.5 + value / 3.0 + value * value / 8.0;
        }
        return (FastMath.exp(value) * (value - 1.0) + 1.0) / (value * value);
    }

    private static void divideComplex(final double numeratorRe,
                                       final double numeratorIm,
                                       final double denominatorRe,
                                       final double denominatorIm,
                                       final MutableComplex out) {
        final double denom = denominatorRe * denominatorRe + denominatorIm * denominatorIm;
        out.re = (numeratorRe * denominatorRe + numeratorIm * denominatorIm) / denom;
        out.im = (numeratorIm * denominatorRe - numeratorRe * denominatorIm) / denom;
    }

    private static void multiplyComplex(final double leftRe,
                                         final double leftIm,
                                         final double rightRe,
                                         final double rightIm,
                                         final MutableComplex out) {
        out.re = leftRe * rightRe - leftIm * rightIm;
        out.im = leftRe * rightIm + leftIm * rightRe;
    }

    private static double[] buildInverseFactorials(final int maxN) {
        final double[] inverse = new double[maxN + 1];
        double factorial = 1.0;
        inverse[0] = 1.0;
        for (int n = 1; n <= maxN; ++n) {
            factorial *= n;
            inverse[n] = 1.0 / factorial;
        }
        return inverse;
    }

    private static double[][] buildBinomialTable(final int maxN) {
        final double[][] binomial = new double[maxN + 1][maxN + 1];
        for (int n = 0; n <= maxN; ++n) {
            binomial[n][0] = 1.0;
            binomial[n][n] = 1.0;
            for (int k = 1; k < n; ++k) {
                binomial[n][k] = binomial[n - 1][k - 1] + binomial[n - 1][k];
            }
        }
        return binomial;
    }
}

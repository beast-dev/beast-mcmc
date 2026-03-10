package dr.evomodel.speciation;

import dr.evolution.tree.*;
import dr.inference.model.*;
import dr.math.FastFourierTransform;
import dr.xml.Reportable;

import java.util.Arrays;

/**
 * @author Frederik M. Andersen
 *
 * Time- and Age-dependent skyline-Weibull birth-death model.
 *
 * Birth and death rates are products of a piecewise-constant scale parameter and a Weibull age-dependent hazard:
 *
 * lambda(a, t) = lambdaScale(t) * lambdaAgeHazard(a)
 * mu(a, t)     = muScale(t)     * muAgeHazard(a)
 *
 * Two solver modes:
 * - FFT-Picard (default): Local convolution approach solved with FFT and Picard iteration
 * - Direct quadrature: anti-diagonal recursion with implicit quadratic solve
 */
public class AgeDependentSkylineBirthDeathModel extends AbstractModelLikelihood implements Reportable {
    private final Tree tree;
    private final Parameter times;
    private final Parameter birthScale;
    private final Parameter birthShape;
    private final Parameter deathScale;
    private final Parameter deathShape;

    private final int numEpochs;
    private final int numSteps;
    private final double h;

    private final double epsPicard;
    private final int maxIterPicard;

    private final double[] p0;
    private final double[] storedP0;
    private final double[] S;
    private final double[] storedS;

    private double logLikelihood = 0.0;
    private double storedLogLikelihood = 0.0;

    // Cached rate and hazard arrays
    private final double[] bScale;
    private final double[] dScale;
    private final double[] bHaz;
    private final double[] dHaz;
    private final double[] bCumHaz;
    private final double[] dCumHaz;

    // Compute flags
    private boolean ratesDirty = true;
    private boolean likelihoodKnown = false;
    private boolean storedLikelihoodKnown = false;

    // Per-epoch work arrays
    private final int[] epochIdx;
    private final double[] dScaleDelta;
    private final double[] bScaleDelta;

    // Solver selection
    private final boolean useDirectQuadrature;

    // Direct quadrature arrays (allocated only when useDirectQuadrature = true)
    private double[] birthRateAtGrid;
    private double[] deathRateAtGrid;
    private boolean[] isChangepoint;
    private double[] leftBirthRate;
    private double[] leftDeathRate;
    private double[] antiDiagBuf;
    private double[] branchDensBuffer;
    private boolean gridValid;
    private int sComputedUpTo;

    /**
     * useDirectQuadrature: if true, use direct quadrature solver; if false, use FFT-Picard
     */
    public AgeDependentSkylineBirthDeathModel(String name,
                                              Tree tree,
                                              Parameter times,
                                              Parameter birthScale,
                                              Parameter birthShape,
                                              Parameter deathScale,
                                              Parameter deathShape,
                                              int numSteps,
                                              double epsPicard,
                                              int maxIterPicard,
                                              boolean useDirectQuadrature) {
        super(name);

        this.tree = tree;
        this.times = times;
        this.birthScale = birthScale;
        this.birthShape = birthShape;
        this.deathScale = deathScale;
        this.deathShape = deathShape;

        this.numEpochs = times.getDimension();
        this.numSteps = numSteps;
        double maxTime = times.getValue(numEpochs - 1);
        this.h = maxTime / numSteps;

        this.epsPicard = epsPicard;
        this.maxIterPicard = maxIterPicard;
        this.useDirectQuadrature = useDirectQuadrature;

        addVariable(times);
        addVariable(birthScale);
        addVariable(birthShape);
        addVariable(deathScale);
        addVariable(deathShape);

        this.p0 = new double[numSteps + 1];
        this.storedP0 = new double[numSteps + 1];
        this.S = new double[numSteps + 1];
        this.storedS = new double[numSteps + 1];

        this.bScale = new double[numEpochs];
        this.dScale = new double[numEpochs];

        this.bHaz = new double[numSteps + 1];
        this.dHaz = new double[numSteps + 1];
        this.bCumHaz = new double[numSteps + 1];
        this.dCumHaz = new double[numSteps + 1];

        this.epochIdx = new int[numEpochs + 1];
        this.dScaleDelta = new double[numEpochs];
        this.bScaleDelta = new double[numEpochs];

        if (useDirectQuadrature) {
            int N = numSteps;
            this.birthRateAtGrid = new double[N + 1];
            this.deathRateAtGrid = new double[N + 1];
            this.isChangepoint = new boolean[N + 1];
            this.leftBirthRate = new double[N + 1];
            this.leftDeathRate = new double[N + 1];
            this.antiDiagBuf = new double[N + 1];
            this.branchDensBuffer = new double[N + 1];
            this.gridValid = false;
        }
    }

    /**
     * Original constructor, defaults to FFT-Picard solver.
     */
    public AgeDependentSkylineBirthDeathModel(String name,
                                              Tree tree,
                                              Parameter times,
                                              Parameter birthScale,
                                              Parameter birthShape,
                                              Parameter deathScale,
                                              Parameter deathShape,
                                              int numSteps,
                                              double epsPicard,
                                              int maxIterPicard) {
        this(name, tree, times, birthScale, birthShape, deathScale, deathShape,
             numSteps, epsPicard, maxIterPicard, false);
    }

    /**
     * Cache epoch indices, piecewise constant scales and age-hazard arrays for Skyline-Weibull rates:
     * r(t, a) = r(t) * h(a) with piecewise constant r(t) and Weibull hazard h(a) = shape * a^(shape - 1).
     */
    private void cacheSkylineWeibull() {
        if (!ratesDirty) return;

        // Epoch scales and indices
        epochIdx[0] = 0;
        bScaleDelta[0] = 0.0;
        dScaleDelta[0] = 0.0;
        for (int k = 0; k < numEpochs; k++) {
            dScale[k] = deathScale.getParameterValue(k);
            bScale[k] = birthScale.getParameterValue(k);
            epochIdx[k + 1] = closestIdx(times.getParameterValue(k));

            if (k == 0) continue;
            bScaleDelta[k] = bScale[k] - bScale[k - 1];
            dScaleDelta[k] = dScale[k] - dScale[k - 1];
        }

        // Age-hazard arrays
        double bShape = birthShape.getParameterValue(0);
        double dShape = deathShape.getParameterValue(0);
        for (int d = 0; d <= numSteps; d++) {
            double age = d * h;

            if (age == 0.0) {
                bHaz[d]    = (bShape == 1.0) ? 1.0 : 0.0;
                dHaz[d]    = (dShape == 1.0) ? 1.0 : 0.0;
                bCumHaz[d] = 0.0;
                dCumHaz[d] = 0.0;
            } else {
                bHaz[d]    = bShape * Math.pow(age, bShape - 1.0);
                dHaz[d]    = dShape * Math.pow(age, dShape - 1.0);
                bCumHaz[d] = Math.pow(age, bShape);
                dCumHaz[d] = Math.pow(age, dShape);
            }
        }

        ratesDirty = false;
    }

    // =========================================================================
    // FFT-Picard solver
    // =========================================================================

    /**
     * Compute extinction probability p0 and branch survival probability S (FFT path)
     */
    private void computeP0AndS() {
        cacheSkylineWeibull();

        // History integrals
        double[] p0History = new double[numSteps + 1];
        double[] SHistory  = new double[numSteps + 1];
        double[] RHistory = new double[numSteps + 1];

        for (int k = 0; k < numEpochs; k++) {
            int startIdx = epochIdx[k];
            int endIdx = epochIdx[k + 1];
            int epochLen = endIdx - startIdx;
            if (k == numEpochs - 1) epochLen++; // Include endpoint in last interval

            double bCurr = bScale[k];
            double dCurr = dScale[k];

            // Epoch boundary history correction
            if (k > 0) {
                for (int g = startIdx; g <= numSteps; g++) {
                    int d = g - startIdx;
                    double deltaR = dScaleDelta[k] * dCumHaz[d] + bScaleDelta[k] * bCumHaz[d];
                    double factor = Math.exp(-deltaR);
                    p0History[g] *= factor;
                    SHistory[g] *= factor;
                    RHistory[g] += deltaR;
                }
            }

            // exp(-R) and mu kernel prefix sum
            double[] localExpmR = new double[numSteps + 1];
            double[] muKernel = new double[numSteps + 1];
            double[] muKernelCum = new double[numSteps + 1];
            for (int d = 0; d <= numSteps; d++) {
                localExpmR[d] = Math.exp(-dCurr * dCumHaz[d] - bCurr * bCumHaz[d]);
                muKernel[d] = dCurr * dHaz[d] * localExpmR[d];
                muKernelCum[d] = (d == 0) ? muKernel[0] : muKernelCum[d - 1] + muKernel[d];
            }

            // Picard kernel: lambda(a) * exp(-R(a))
            double[] kernelFFT = buildPicardKernel(bCurr, dCurr, epochLen);
            int picardPad = kernelFFT.length - 2;
            double K0 = kernelFFT[picardPad]; // stored by buildPicardKernel

            // =======================
            // SOLVE p0 FOR EPOCH k
            // =======================

            // Boundary function p0B = dLocalInt + p0History
            double[] p0B = new double[epochLen];

            // Local death integral via prefix sum (trapezoidal rule)
            p0B[0] = 0.0;
            for (int i = 1; i < epochLen; i++) {
                p0B[i] = h * (muKernelCum[i] - 0.5 * muKernel[0] - 0.5 * muKernel[i]);
            }

            // Add history to boundary function
            if (k > 0) {
                for (int i = 0; i < epochLen; i++) {
                    p0B[i] += p0History[startIdx + i];
                }
            }

            // Picard iteration for p0
            double[] guessCurr = new double[epochLen];
            System.arraycopy(p0B, 0, guessCurr, 0, epochLen);

            double[] fBuf = new double[picardPad + 2];
            double[] guessNext = new double[epochLen];

            boolean converged = false;
            int iter = 0;
            while (!converged && iter < maxIterPicard) {
                Arrays.fill(fBuf, 0.0);
                for (int i = 0; i < epochLen; i++) {
                    fBuf[i] = guessCurr[i] * guessCurr[i];
                }
                fBuf[0] *= 0.5; // Left-endpoint trapezoidal correction
                FastFourierTransform.rfft(fBuf, picardPad, false);
                convolve(kernelFFT, fBuf, picardPad);

                double maxDiff = 0.0;
                for (int i = 0; i < epochLen; i++) {
                    // Trapezoidal endpoint correction
                    guessNext[i] = p0B[i] + fBuf[i] * h - 0.5 * h * K0 * guessCurr[i] * guessCurr[i];
                    maxDiff = Math.max(maxDiff, Math.abs(guessNext[i] - guessCurr[i]));
                }

                double[] tmp = guessCurr;
                guessCurr = guessNext;
                guessNext = tmp;

                if (maxDiff < epsPicard) converged = true;
                iter++;
            }

            System.arraycopy(guessCurr, 0, p0, startIdx, epochLen);

            // =======================
            // SOLVE S FOR EPOCH k
            // =======================

            // Boundary function SB = exp(-R) + SHistory
            double[] SB = new double[epochLen];

            for (int i = 0; i < epochLen; i++) {
                int a = startIdx + i;
                SB[i] = Math.exp(- RHistory[a] - dScale[0] * dCumHaz[a] - bScale[0] * bCumHaz[a]);

                if (k > 0) SB[i] += SHistory[startIdx + i];
            }

            // Picard iteration for S
            System.arraycopy(SB, 0, guessCurr, 0, epochLen);

            converged = false;
            iter = 0;
            while (!converged && iter < maxIterPicard) {
                Arrays.fill(fBuf, 0.0);
                for (int i = 0; i < epochLen; i++) {
                    fBuf[i] = 2.0 * p0[startIdx + i] * guessCurr[i];
                }
                fBuf[0] *= 0.5; // Left-endpoint trapezoidal correction
                FastFourierTransform.rfft(fBuf, picardPad, false);
                convolve(kernelFFT, fBuf, picardPad);

                double maxDiff = 0.0;
                for (int i = 0; i < epochLen; i++) {
                    double fi = 2.0 * p0[startIdx + i] * guessCurr[i];
                    guessNext[i] = SB[i] + fBuf[i] * h - 0.5 * h * K0 * fi;
                    maxDiff = Math.max(maxDiff, Math.abs(guessNext[i] - guessCurr[i]));
                }

                double[] tmp = guessCurr;
                guessCurr = guessNext;
                guessNext    = tmp;

                if (maxDiff < epsPicard) converged = true;
                iter++;
            }

            System.arraycopy(guessCurr, 0, S, startIdx, epochLen);

            // ==================
            // COMPUTE HISTORY
            // ==================
            if (k < numEpochs - 1) {
                // Death contribution to p0 history via trapezoidal sum
                for (int g = startIdx; g <= numSteps; g++) {
                    int pos = g - startIdx;
                    int lo = Math.max(0, pos - epochLen + 1);
                    double windowSum = muKernelCum[pos] - (lo > 0 ? muKernelCum[lo - 1] : 0.0);
                    double trapCorr = (startIdx == 0) ? 0.5 * muKernel[pos] : 0.0;
                    p0History[g] += h * windowSum - h * trapCorr;
                }

                // Birth kernel - Transformed
                int historyPad = nextPow2(numSteps + epochLen);
                double[] lambdaKernel = new double[historyPad + 2];
                for (int d = 0; d <= numSteps; d++) {
                    lambdaKernel[d] = bCurr * bHaz[d] * localExpmR[d];
                }
                FastFourierTransform.rfft(lambdaKernel, historyPad, false);

                // Birth contribution to p0 history signal - Transformed
                double[] p0Signal = new double[historyPad + 2];
                for (int idx = 0; idx < epochLen; idx++) {
                    int pastIdx = startIdx + idx;
                    double trapWeight = (pastIdx == 0) ? h * 0.5 : h;
                    p0Signal[idx] = p0[pastIdx] * p0[pastIdx] * trapWeight;
                }
                FastFourierTransform.rfft(p0Signal, historyPad, false);

                // Convolve and add to history
                convolve(lambdaKernel, p0Signal, historyPad);
                for (int pos = 0; pos < historyPad; pos++) {
                    int g = startIdx + pos;
                    if (g > numSteps) break;
                    p0History[g] += p0Signal[pos];
                }

                // S history signal
                double[] SSignal = new double[historyPad + 2];
                for (int idx = 0; idx < epochLen; idx++) {
                    int pastIdx = startIdx + idx;
                    double trapWeight = (pastIdx == 0) ? h * 0.5 : h;
                    SSignal[idx] = 2.0 * p0[pastIdx] * S[pastIdx] * trapWeight;
                }
                FastFourierTransform.rfft(SSignal, historyPad, false);

                // Convolve and add to history
                convolve(lambdaKernel, SSignal, historyPad);
                for (int pos = 0; pos < historyPad; pos++) {
                    int g = startIdx + pos;
                    if (g > numSteps) break;
                    SHistory[g] += SSignal[pos];
                }
            }
        }
    }


    /**
     * Computes branch length density g_t(ell) for a branch with node height t and length ell (FFT path)
     */
    public double branchLenDens(double t, double ell) {
        // Find starting epoch containing t
        int tIdx = closestIdx(t);

        // Find index first index after ell such that we can interpolate
        // Cap so we never read past the last computed index (numSteps - tIdx)
        int maxEll = Math.min((int) Math.ceil(ell / h) + 1, numSteps - tIdx);
        double[] gt = new double[maxEll + 1];
        int kStart = 0;
        for (int m = 1; m < numEpochs; m++) {
            if (tIdx >= epochIdx[m]) kStart = m;
        }

        // History integrals
        double[] gHistory = new double[maxEll + 1];
        double[] RHistory = new double[maxEll + 1];

        double p0t = interpolate(p0, t);

        for (int k = kStart; k < numEpochs; k++) {
            // Local ell range in this epoch
            int epochEndIdx = epochIdx[k + 1];
            if (epochEndIdx > numSteps) epochEndIdx = numSteps;

            int ellStart = (k == kStart) ? 0 : epochIdx[k] - tIdx;
            int ellEnd;
            if (k < numEpochs - 1) {
                ellEnd = epochEndIdx - tIdx - 1;
            } else {
                ellEnd = numSteps - tIdx;
            }
            ellEnd = Math.min(ellEnd, maxEll);
            if (ellEnd < ellStart || ellStart > maxEll) break;

            int ellLen = ellEnd - ellStart + 1;
            if (ellLen <= 0) continue;

            // Epoch boundary history correction
            if (k > kStart) {
                int ellBound = epochIdx[k] - tIdx;
                for (int l = ellBound; l <= maxEll; l++) {
                    int d = l - ellBound;
                    d = Math.min(d, numSteps);
                    double deltaR = dScaleDelta[k] * dCumHaz[d] + bScaleDelta[k] * bCumHaz[d];
                    gHistory[l] *= Math.exp(-deltaR);
                    RHistory[l] += deltaR;
                }
            }

            // Boundary function gB = lambda exp(-R) (1-p0)^2
            double[] gB = new double[ellLen];
            for (int i = 0; i < ellLen; i++) {
                int l = ellStart + i;
                gB[i] = bScale[kStart] * bHaz[l]
                        * Math.exp(-dScale[kStart] * dCumHaz[l]
                                  - bScale[kStart] * bCumHaz[l]
                                  - RHistory[l])
                        * (1.0 - p0t) * (1.0 - p0t);

                if (k > kStart) gB[i] += gHistory[ellStart + i];
            }

            // Picard kernel: lambda(a) * exp(-R(a))
            double bCurr = bScale[k];
            double dCurr = dScale[k];

            double[] kernelReal = buildPicardKernel(bCurr, dCurr, ellLen);
            int picardPad = kernelReal.length - 2;
            double gtK0 = kernelReal[picardPad];

            // Picard iteration for g_t
            double[] currentGuess = new double[ellLen];
            System.arraycopy(gB, 0, currentGuess, 0, ellLen);

            double[] fBuf = new double[picardPad + 2];
            double[] nextGuess = new double[ellLen];

            boolean converged = false;
            int iter = 0;
            while (!converged && iter < maxIterPicard) {
                Arrays.fill(fBuf, 0.0);
                for (int i = 0; i < ellLen; i++) {
                    fBuf[i] = 2.0 * p0[tIdx + ellStart + i] * currentGuess[i];
                }
                fBuf[0] *= 0.5; // Left-endpoint trapezoidal correction
                FastFourierTransform.rfft(fBuf, picardPad, false);
                convolve(kernelReal, fBuf, picardPad);

                double maxDiff = 0.0;
                for (int i = 0; i < ellLen; i++) {
                    double fi = 2.0 * p0[tIdx + ellStart + i] * currentGuess[i];
                    nextGuess[i] = gB[i] + fBuf[i] * h - 0.5 * h * gtK0 * fi;
                    maxDiff = Math.max(maxDiff, Math.abs(nextGuess[i] - currentGuess[i]));
                }

                double[] tmp = currentGuess;
                currentGuess = nextGuess;
                nextGuess    = tmp;

                if (maxDiff < epsPicard) converged = true;
                iter++;
            }

            System.arraycopy(currentGuess, 0, gt, ellStart, ellLen);

            // Convolve and add to history
            if (k < numEpochs - 1 && ellEnd < maxEll) {
                int historyLen = ellLen + maxEll - ellEnd;
                int historyPad = nextPow2(historyLen + ellLen - 1);

                double[] lambdaKernel = new double[historyPad + 2];
                for (int d = 0; d < historyLen && d <= numSteps; d++) {
                    lambdaKernel[d] = bCurr * bHaz[d] * Math.exp(-dCurr * dCumHaz[d] - bCurr * bCumHaz[d]);
                }
                FastFourierTransform.rfft(lambdaKernel, historyPad, false);

                double[] sigLambdaGt = new double[historyPad + 2];
                for (int idx = 0; idx < ellLen; idx++) {
                    int l = ellStart + idx;
                    double trapWeight = (l == 0) ? h * 0.5 : h;
                    sigLambdaGt[idx] = 2.0 * p0[tIdx + l] * gt[l] * trapWeight;
                }
                FastFourierTransform.rfft(sigLambdaGt, historyPad, false);
                convolve(lambdaKernel, sigLambdaGt, historyPad);

                for (int pos = 0; pos < historyPad; pos++) {
                    int l = ellStart + pos;
                    if (l > maxEll) break;
                    gHistory[l] += sigLambdaGt[pos];
                }
            }
        }

        return interpolate(gt, ell);
    }

    /**
     * FFT-Picard likelihood computation.
     */
    private double calculateLogLikelihoodFFT() {
        computeP0AndS();

        double originTime = times.getParameterValue(times.getDimension() - 1);
        double logL = 0.0;

        // Internal nodes, computes branch length densities
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            NodeRef node = tree.getInternalNode(i);
            double nodeHeight = tree.getNodeHeight(node);

            double branchLength;
            if (tree.isRoot(node)) {
                branchLength = originTime - nodeHeight;
            } else {
                branchLength = tree.getBranchLength(node);
            }
            branchLength = Math.max(branchLength, h);

            double gt = branchLenDens(nodeHeight, branchLength);

            double parentHeight = nodeHeight + branchLength;
            double p0Parent = interpolate(p0, parentHeight);

            logL += Math.log(gt) - Math.log(1.0 - p0Parent);
        }

        // External nodes, evaluates survival probabilities
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            NodeRef node = tree.getExternalNode(i);
            double branchLength = Math.max(tree.getBranchLength(node), h);
            double parentHeight = tree.getNodeHeight(node) + branchLength;

            double sVal = interpolate(S, branchLength);
            double p0Parent = interpolate(p0, parentHeight);

            logL += Math.log(sVal) - Math.log(1.0 - p0Parent);
        }

        if (Double.isNaN(logL) || logL == Double.POSITIVE_INFINITY) {
            return Double.NEGATIVE_INFINITY;
        }

        return logL;
    }

    // =========================================================================
    // Direct quadrature solver
    // =========================================================================

    /**
     * Build grid arrays for direct quadrature: map epoch boundaries to grid points,
     * fill rate arrays, flag changepoints.
     */
    private void buildGridDirect() {
        int N = numSteps;

        for (int i = 0; i <= N; i++) {
            // Find which epoch this grid point belongs to
            int epoch = 0;
            for (int k = 1; k < numEpochs; k++) {
                if (i >= epochIdx[k]) {
                    epoch = k;
                }
            }
            int birthEpoch = Math.min(epoch, numEpochs - 1);
            int deathEpoch = Math.min(epoch, deathScale.getDimension() - 1);

            birthRateAtGrid[i] = bScale[birthEpoch];
            deathRateAtGrid[i] = dScale[deathEpoch];

            // Detect changepoint
            if (i > 0) {
                int prevEpoch = 0;
                for (int k = 1; k < numEpochs; k++) {
                    if (i - 1 >= epochIdx[k]) {
                        prevEpoch = k;
                    }
                }
                int prevBirthEpoch = Math.min(prevEpoch, numEpochs - 1);
                int prevDeathEpoch = Math.min(prevEpoch, deathScale.getDimension() - 1);

                if (prevBirthEpoch != birthEpoch || prevDeathEpoch != deathEpoch) {
                    isChangepoint[i] = true;
                    leftBirthRate[i] = bScale[prevBirthEpoch];
                    leftDeathRate[i] = dScale[prevDeathEpoch];
                } else {
                    isChangepoint[i] = false;
                }
            } else {
                isChangepoint[i] = false;
            }
        }
    }

    /**
     * Compute expmR values along the anti-diagonal i + j = C.
     * buf[i] = exp(-sum of (b(s)*deltaBirthH + d(s)*deltaDeathH) from s=i to C-1)
     */
    private void computeAntiDiagonal(int C, double[] buf) {
        double logExpmR = 0.0;
        buf[C] = 1.0;

        for (int i = C - 1; i >= 0; i--) {
            int j = C - i;
            double bRate = getEffectiveBirthRate(i);
            double dRate = getEffectiveDeathRate(i);
            double deltaBH = bCumHaz[j] - bCumHaz[j - 1];
            double deltaDH = dCumHaz[j] - dCumHaz[j - 1];

            logExpmR -= bRate * deltaBH + dRate * deltaDH;
            buf[i] = Math.exp(logExpmR);
        }
    }

    // Vanilla: no changepoint averaging, just use the grid rate directly
    private double getEffectiveBirthRate(int i) {
        return birthRateAtGrid[i];
    }

    private double getEffectiveDeathRate(int i) {
        return deathRateAtGrid[i];
    }

    // Changepoint-averaged versions (commented out):
    // private double getEffectiveBirthRate(int i) {
    //     if (isChangepoint[i]) {
    //         return (leftBirthRate[i] + birthRateAtGrid[i]) / 2.0;
    //     }
    //     return birthRateAtGrid[i];
    // }
    //
    // private double getEffectiveDeathRate(int i) {
    //     if (isChangepoint[i]) {
    //         return (leftDeathRate[i] + deathRateAtGrid[i]) / 2.0;
    //     }
    //     return deathRateAtGrid[i];
    // }

    /**
     * Full birth rate at grid point (i, j): b(i) * h_b(j*h)
     */
    private double getBirthRate(int i, int j) {
        return getEffectiveBirthRate(i) * bHaz[j];
    }

    /**
     * Full death rate at grid point (i, j): d(i) * h_d(j*h)
     */
    private double getDeathRate(int i, int j) {
        return getEffectiveDeathRate(i) * dHaz[j];
    }

    /**
     * Vanilla: birth rate at age 0, just use point evaluation h_b(0).
     */
    private double getBirthRateAvgZero(int i) {
        return getBirthRate(i, 0);
    }

    /**
     * Vanilla: death rate at age 0, just use point evaluation h_d(0).
     */
    private double getDeathRateAvgZero(int i) {
        return getDeathRate(i, 0);
    }

    // Integral-averaged age-0 versions (commented out):
    // private double getBirthRateAvgZero(int i) {
    //     double deltaBH = bCumHaz[1] - bCumHaz[0]; // = bCumHaz[1] since bCumHaz[0]=0
    //     return getEffectiveBirthRate(i) * deltaBH / h;
    // }
    //
    // private double getDeathRateAvgZero(int i) {
    //     double deltaDH = dCumHaz[1] - dCumHaz[0];
    //     return getEffectiveDeathRate(i) * deltaDH / h;
    // }

    /**
     * Compute extinction probabilities using direct quadrature.
     */
    private void computeP0Direct() {
        p0[0] = 0.0;
        S[0] = 1.0;
        sComputedUpTo = 0;

        for (int m = 1; m <= numSteps; m++) {
            computeAntiDiagonal(m, antiDiagBuf);

            double trap_sum_mu = 0.0;
            double trap_sum_lam = 0.0;

            // i=0 endpoint (half weight)
            trap_sum_mu += getDeathRate(0, m) * antiDiagBuf[0] / 2.0;

            // Interior points
            for (int i = 1; i < m; i++) {
                int j = m - i;
                double expmR_ij = antiDiagBuf[i];
                trap_sum_mu += getDeathRate(i, j) * expmR_ij;
                trap_sum_lam += getBirthRate(i, j) * expmR_ij * p0[i] * p0[i];
            }

            // i=m endpoint (half weight for mu only)
            double expmR_m0 = antiDiagBuf[m];
            trap_sum_mu += getDeathRateAvgZero(m) * expmR_m0 / 2.0;

            trap_sum_mu *= h;
            trap_sum_lam *= h;

            // Implicit quadratic for p0[m]
            double a = h * getBirthRateAvgZero(m) * expmR_m0 / 2.0;
            double c = trap_sum_mu + trap_sum_lam;

            if (a == 0.0) {
                p0[m] = c;
            } else {
                double det = 1.0 - 4.0 * a * c;
                if (det < 0.0) {
                    p0[m] = 1.0 / (2.0 * a);
                } else {
                    p0[m] = (1.0 - Math.sqrt(det)) / (2.0 * a);
                }
            }
        }

        gridValid = true;
    }

    /**
     * Compute branch survival probability up to index idx (direct quadrature).
     */
    private void computeSUpToDirect(int idx) {
        if (idx <= sComputedUpTo) return;

        for (int i = sComputedUpTo + 1; i <= idx; i++) {
            computeAntiDiagonal(i, antiDiagBuf);

            double expmR_0i = antiDiagBuf[0];
            double trap_sum = 0.0;

            for (int j = 1; j < i; j++) {
                int k = i - j;
                trap_sum += getBirthRate(j, k) * antiDiagBuf[j] * p0[j] * S[j];
            }

            S[i] = (2.0 * h * trap_sum + expmR_0i) /
                    (1.0 - h * getBirthRateAvgZero(i) * antiDiagBuf[i] * p0[i]);
        }

        sComputedUpTo = idx;
    }

    /**
     * Compute branch length density for a node at height k with branch length l (direct quadrature).
     */
    private double branchLengthDensityDirect(int k, int l) {
        double p0_k = p0[k];
        double one_minus_p0_k = 1.0 - p0_k;

        branchDensBuffer[0] = getBirthRateAvgZero(k) * one_minus_p0_k * one_minus_p0_k;

        for (int m = 1; m <= l; m++) {
            int C = k + m;
            computeAntiDiagonal(C, antiDiagBuf);

            double birthRate_k_m = getBirthRate(k, m);
            double expmR_k_m = antiDiagBuf[k];
            double trap_sum = birthRate_k_m * expmR_k_m * p0_k * branchDensBuffer[0] / 2.0;

            for (int r = 1; r < m; r++) {
                int i = k + r;
                int j = m - r;
                trap_sum += getBirthRate(i, j) * antiDiagBuf[i] * p0[i] * branchDensBuffer[r];
            }

            int i_end = k + m;
            double num = 2.0 * h * trap_sum + birthRate_k_m * expmR_k_m * one_minus_p0_k * one_minus_p0_k;
            double denom = 1.0 - h * getBirthRateAvgZero(i_end) * antiDiagBuf[i_end] * p0[i_end];

            branchDensBuffer[m] = num / denom;
        }

        return Math.log(branchDensBuffer[l]) - Math.log(1.0 - p0[k + l]);
    }

    /**
     * Direct quadrature likelihood computation.
     */
    private double calculateLogLikelihoodDirect() {
        cacheSkylineWeibull();

        if (!gridValid) {
            buildGridDirect();
            computeP0Direct();
        }

        double logL = 0.0;

        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            NodeRef node = tree.getInternalNode(i);

            if (tree.isRoot(node)) {
                int k = (int) Math.round(tree.getNodeHeight(node) / h);
                logL += branchLengthDensityDirect(k, numSteps - k);
                continue;
            }

            int k = (int) Math.round(tree.getNodeHeight(node) / h);
            int l = Math.max(1, (int) Math.round(tree.getBranchLength(node) / h));
            logL += branchLengthDensityDirect(k, l);
        }

        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            NodeRef node = tree.getExternalNode(i);

            int l = Math.max(1, (int) Math.round(tree.getBranchLength(node) / h));
            computeSUpToDirect(l);
            logL += Math.log(S[l]) - Math.log(1.0 - p0[l]);
        }

        if (Double.isNaN(logL) || logL == Double.POSITIVE_INFINITY) {
            return Double.NEGATIVE_INFINITY;
        }

        return logL;
    }

    // =======
    // Utils
    // =======

    /**
     * Use appropriate solver.
     */
    private double calculateLogLikelihood() {
        if (useDirectQuadrature) {
            return calculateLogLikelihoodDirect();
        } else {
            return calculateLogLikelihoodFFT();
        }
    }

    /**
     * Build Picard kernel and transform it. Store K0 for trapezoidal correction.
     */
    private double[] buildPicardKernel(double bScale, double dScale, int M) {
        int picardPad = nextPow2(2 * M - 1);
        double[] kernel = new double[picardPad + 2];
        for (int i = 0; i < M; i++) {
            kernel[i] = bScale * bHaz[i] * Math.exp(-dScale * dCumHaz[i] - bScale * bCumHaz[i]);
        }
        double K0 = kernel[0];
        FastFourierTransform.rfft(kernel, picardPad, false);
        kernel[picardPad] = K0;
        return kernel;
    }

    /**
     * Pointwise complex multiply transformed kernel and signal, and inverse-transform.
     */
    private static void convolve(double[] kernel, double[] signal, int padM) {
        signal[0] = kernel[0] * signal[0];
        signal[1] = kernel[1] * signal[1];
        for (int i = 1; i < padM / 2; i++) {
            double kr = kernel[2 * i],  ki = kernel[2 * i + 1];
            double sr = signal[2 * i],  si = signal[2 * i + 1];
            signal[2 * i]     = kr * sr - ki * si;
            signal[2 * i + 1] = kr * si + ki * sr;
        }
        FastFourierTransform.rfft(signal, padM, true);
    }

    /**
     * Find next power of two
     */
    private static int nextPow2(int n) {
        if (n <= 1) return 2;
        int p = Integer.highestOneBit(n - 1) << 1;
        return Math.max(p, 2);
    }

    /**
     * Find the closest index in grid
     */
    private int closestIdx(double t) {
        int idx = (int) Math.round(t / h);
        if (idx < 0) idx = 0;
        if (idx > numSteps) idx = numSteps;

        return idx;
    }

    /**
     * Interpolates a grid array at a continuous time value
     */
    private double interpolate(double[] grid, double t) {
        double idx = t / h;
        int lo = (int) Math.floor(idx);
        if (lo < 0) lo = 0;
        if (lo >= grid.length - 1) return grid[grid.length - 1];
        double frac = idx - lo;
        return (1.0 - frac) * grid[lo] + frac * grid[lo + 1];
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        makeDirty();
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        makeDirty();
    }

    protected void storeState() {
        System.arraycopy(p0, 0, storedP0, 0, p0.length);
        System.arraycopy(S, 0, storedS, 0, S.length);
        storedLikelihoodKnown = likelihoodKnown;
        storedLogLikelihood = logLikelihood;
    }

    protected void restoreState() {
        System.arraycopy(storedP0, 0, p0, 0, p0.length);
        System.arraycopy(storedS, 0, S, 0, S.length);
        likelihoodKnown = storedLikelihoodKnown;
        logLikelihood = storedLogLikelihood;

        ratesDirty = true;
        if (useDirectQuadrature) {
            gridValid = false;
        }
    }

    protected void acceptState() {
    }

    public Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    public void makeDirty() {
        ratesDirty = true;
        likelihoodKnown = false;
        if (useDirectQuadrature) {
            gridValid = false;
        }
    }

    /**
     * Returns the full p0 grid array only for testing purpose
     */
    public double[] getP0Array() {
        if (useDirectQuadrature) {
            cacheSkylineWeibull();
            if (!gridValid) {
                buildGridDirect();
                computeP0Direct();
            }
        } else {
            computeP0AndS();
        }
        return p0.clone();
    }

    /**
     * Returns the full S grid array only for testing purpose
     */
    public double[] getSArray() {
        computeP0AndS();
        return S.clone();
    }

    public String getReport() {
        return "logLikelihood: " + getLogLikelihood();
    }

    public String toString() {
        return Double.toString(getLogLikelihood());
    }
}

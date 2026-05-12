package dr.evomodel.speciation;

import dr.evolution.tree.*;
import dr.evomodel.tree.TreeChangedEvent;
import dr.inference.model.*;
import dr.math.FastFourierTransform;
import dr.xml.Reportable;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.*;

/**
 * @author Frederik M. Andersen
 *
 * Time- and Age-dependent skyline birth-death model with linear-exponential age hazard.
 *
 * Birth and death rates are products of a piecewise-constant scale parameter and a
 * linear-exponential age-dependent hazard:
 *
 * lambda(a, t) = lambdaScale(t) * h_b(a)
 * mu(a, t)     = muScale(t)     * h_d(a)
 *
 * where h(a) = (1 + b*a) * exp(-gamma*a), with shape parameters [b, gamma].
 *
 * Two solver modes:
 * - FFT-Picard (default): Local convolution approach solved with FFT and Picard iteration
 * - Direct quadrature: anti-diagonal recursion with implicit quadratic solve
 *
 * Per-node caching: when only the tree changes (node height proposals), only affected nodes
 * and their direct children are recomputed. The p0/S grids (tree-independent) are reused.
 */
public class AgeDependentBirthDeathIEModelOld extends AbstractModelLikelihood implements Reportable {
    private final Tree tree;
    private final Parameter times;
    private final Parameter birthScale;
    private final Parameter birthShape;
    private final Parameter deathScale;
    private final Parameter deathShape;

    private final int numEpochs;
    private final int maxSteps;
    private int numSteps;
    private double h;

    private final double epsPicard;
    private final int maxIterPicard;

    private final double[] p0;
    private final double[] S;

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

    // Likelihood variant: when true, drop the root-branch term (from root height to origin).
    // Implicitly condition on non-extinction of both subtrees at root height, which falls out
    // automatically from the per-node -log(1 - p0(parentHeight)) terms of root's two children.
    private final boolean excludeRootBranch;

    // Direct quadrature arrays (allocated only when useDirectQuadrature = true)
    private double[] birthRateAtGrid;
    private double[] deathRateAtGrid;
    private double[] antiDiagBuf;
    private double[] branchDensBuffer;
    private boolean gridValid;
    private int sComputedUpTo;

    // Per-node caching
    private double[] cachedNodeLogL;
    private double[] storedCachedNodeLogL;
    private boolean[] nodeValid;
    private boolean[] storedNodeValid;

    // Cached p0/S for tree-only proposals
    private double[] storedP0;
    private double[] storedS;

    private boolean parametersDirty = true;
    private boolean storedParametersDirty = true;

    // Parallel computation
    private final int numThreads;
    private final ExecutorService threadPool;

    public AgeDependentBirthDeathIEModelOld(String name,
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
        this(name, tree, times, birthScale, birthShape, deathScale, deathShape,
             numSteps, epsPicard, maxIterPicard, useDirectQuadrature, 1, false);
    }

    public AgeDependentBirthDeathIEModelOld(String name,
                                            Tree tree,
                                            Parameter times,
                                            Parameter birthScale,
                                            Parameter birthShape,
                                            Parameter deathScale,
                                            Parameter deathShape,
                                            int numSteps,
                                            double epsPicard,
                                            int maxIterPicard,
                                            boolean useDirectQuadrature,
                                            int numThreads) {
        this(name, tree, times, birthScale, birthShape, deathScale, deathShape,
             numSteps, epsPicard, maxIterPicard, useDirectQuadrature, numThreads, false);
    }

    public AgeDependentBirthDeathIEModelOld(String name,
                                            Tree tree,
                                            Parameter times,
                                            Parameter birthScale,
                                            Parameter birthShape,
                                            Parameter deathScale,
                                            Parameter deathShape,
                                            int numSteps,
                                            double epsPicard,
                                            int maxIterPicard,
                                            boolean useDirectQuadrature,
                                            int numThreads,
                                            boolean excludeRootBranch) {
        super(name);

        this.tree = tree;
        this.times = times;
        this.birthScale = birthScale;
        this.birthShape = birthShape;
        this.deathScale = deathScale;
        this.deathShape = deathShape;

        this.numEpochs = times.getDimension();
        this.maxSteps = numSteps;
        this.numSteps = numSteps;

        this.epsPicard = epsPicard;
        this.maxIterPicard = maxIterPicard;
        this.useDirectQuadrature = useDirectQuadrature;
        this.excludeRootBranch = excludeRootBranch;
        this.numThreads = numThreads;
        this.threadPool = (numThreads > 1) ? Executors.newFixedThreadPool(numThreads) : null;

        if (tree instanceof Model) {
            addModel((Model) tree);
        }

        addVariable(times);
        addVariable(birthScale);
        addVariable(birthShape);
        addVariable(deathScale);
        addVariable(deathShape);

        this.p0 = new double[maxSteps + 1];
        this.S = new double[maxSteps + 1];

        this.bScale = new double[numEpochs];
        this.dScale = new double[numEpochs];

        this.bHaz = new double[maxSteps + 1];
        this.dHaz = new double[maxSteps + 1];
        this.bCumHaz = new double[maxSteps + 1];
        this.dCumHaz = new double[maxSteps + 1];

        this.epochIdx = new int[numEpochs + 1];
        this.dScaleDelta = new double[numEpochs];
        this.bScaleDelta = new double[numEpochs];

        if (useDirectQuadrature) {
            this.birthRateAtGrid = new double[maxSteps + 1];
            this.deathRateAtGrid = new double[maxSteps + 1];
            this.antiDiagBuf = new double[maxSteps + 1];
            this.branchDensBuffer = new double[maxSteps + 1];
            this.gridValid = false;
        }

        // Per-node caching
        int totalNodes = tree.getNodeCount();
        this.cachedNodeLogL = new double[totalNodes];
        this.storedCachedNodeLogL = new double[totalNodes];
        this.nodeValid = new boolean[totalNodes];
        this.storedNodeValid = new boolean[totalNodes];
        this.storedP0 = new double[maxSteps + 1];
        this.storedS = new double[maxSteps + 1];
    }

    public AgeDependentBirthDeathIEModelOld(String name,
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

    // =========================================================================
    // Node invalidation
    // =========================================================================

    private void invalidateNode(NodeRef node) {
        nodeValid[node.getNumber()] = false;
        for (int i = 0; i < tree.getChildCount(node); i++) {
            nodeValid[tree.getChild(node, i).getNumber()] = false;
        }
    }

    private void invalidateAllNodes() {
        Arrays.fill(nodeValid, false);
    }

    /**
     * Cache epoch indices, piecewise constant scales and age-hazard arrays.
     * Linear-exponential hazard: h(a) = (1 + b*a) * exp(-gamma*a)
     * Cumulative hazard H(a) = integral_0^a h(s) ds
     */
    private void cacheRates() {
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

        // Age-hazard arrays: linear-exponential h(a) = (1 + b*a) * exp(-gamma*a)
        double bB = birthShape.getParameterValue(0);
        double bGamma = birthShape.getParameterValue(1);
        double dB = deathShape.getParameterValue(0);
        double dGamma = deathShape.getParameterValue(1);

        for (int d = 0; d <= numSteps; d++) {
            double age = d * h;
            double expBG = Math.exp(-bGamma * age);
            double expDG = Math.exp(-dGamma * age);

            bHaz[d] = (1.0 + bB * age) * expBG;
            dHaz[d] = (1.0 + dB * age) * expDG;

            bCumHaz[d] = linExpCumHaz(age, bB, bGamma);
            dCumHaz[d] = linExpCumHaz(age, dB, dGamma);
        }

        ratesDirty = false;
    }

    /**
     * Cumulative hazard H(a) = integral_0^a (1 + b*s) * exp(-gamma*s) ds
     * gamma > 0: H(a) = (1/gamma + b/gamma^2)(1 - exp(-gamma*a)) - (b*a/gamma) * exp(-gamma*a)
     * gamma = 0: H(a) = a + b*a^2/2
     */
    private static double linExpCumHaz(double a, double b, double gamma) {
        if (gamma == 0.0) {
            return a + b * a * a / 2.0;
        }
        double emga = Math.exp(-gamma * a);
        double invG = 1.0 / gamma;
        return (invG + b * invG * invG) * (1.0 - emga) - b * a * invG * emga;
    }

    // =========================================================================
    // FFT-Picard solver
    // =========================================================================

    /**
     * Compute extinction probability p0 and branch survival probability S (FFT path)
     */
    private void computeP0AndS() {
        cacheRates();

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
                fBuf[0] *= 0.5; // Left-endpoint trapezoidal weight
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
                fBuf[0] *= 0.5; // Left-endpoint trapezoidal weight
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
                fBuf[0] *= 0.5; // Left-endpoint trapezoidal weight
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
     * FFT-Picard likelihood computation (non-cached path)
     */
    private double calculateLogLikelihoodFFT() {
        computeP0AndS();

        double originTime = times.getParameterValue(times.getDimension() - 1);
        double logL = 0.0;

        if (numThreads > 1) {
            // Collect internal node tasks
            List<double[]> tasks = new ArrayList<>();
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
                tasks.add(new double[]{nodeHeight, branchLength});
            }

            List<Future<Double>> futures = new ArrayList<>(tasks.size());
            for (double[] task : tasks) {
                final double nodeHeight = task[0];
                final double branchLength = task[1];
                futures.add(threadPool.submit(() -> {
                    double gt = branchLenDens(nodeHeight, branchLength);
                    double parentHeight = nodeHeight + branchLength;
                    double p0Parent = interpolate(p0, parentHeight);
                    return Math.log(gt) - Math.log(1.0 - p0Parent);
                }));
            }

            for (Future<Double> f : futures) {
                try {
                    logL += f.get();
                } catch (InterruptedException | ExecutionException e) {
                    throw new RuntimeException("Parallel node computation failed", e);
                }
            }
        } else {
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
    // Per-node contribution methods for node caching
    // =========================================================================

    /**
     * Compute per-node log-likelihood contributions (FFT path).
     * Only computes contributions for nodes where nodeValid[nodeNum] is false.
     * When numThreads > 1, internal node contributions are computed in parallel.
     */
    private void computeNodeContributionsFFT(double[] contributions) {
        double originTime = times.getParameterValue(times.getDimension() - 1);

        if (numThreads > 1) {
            // Collect invalidated internal nodes
            List<Future<?>> futures = new ArrayList<>();
            for (int i = 0; i < tree.getInternalNodeCount(); i++) {
                NodeRef node = tree.getInternalNode(i);
                final int nodeNum = node.getNumber();
                if (nodeValid[nodeNum]) continue;

                if (excludeRootBranch && tree.isRoot(node)) {
                    contributions[nodeNum] = 0.0;
                    continue;
                }

                final double nodeHeight = tree.getNodeHeight(node);
                double bl;
                if (tree.isRoot(node)) {
                    bl = originTime - nodeHeight;
                } else {
                    bl = tree.getBranchLength(node);
                }
                final double branchLength = Math.max(bl, h);

                futures.add(threadPool.submit(() -> {
                    double gt = branchLenDens(nodeHeight, branchLength);
                    double parentHeight = nodeHeight + branchLength;
                    double p0Parent = interpolate(p0, parentHeight);
                    contributions[nodeNum] = Math.log(gt) - Math.log(1.0 - p0Parent);
                }));
            }

            // Wait for all tasks
            for (Future<?> f : futures) {
                try {
                    f.get();
                } catch (InterruptedException | ExecutionException e) {
                    throw new RuntimeException("Parallel node computation failed", e);
                }
            }
        } else {
            for (int i = 0; i < tree.getInternalNodeCount(); i++) {
                NodeRef node = tree.getInternalNode(i);
                int nodeNum = node.getNumber();
                if (nodeValid[nodeNum]) continue;

                if (excludeRootBranch && tree.isRoot(node)) {
                    contributions[nodeNum] = 0.0;
                    continue;
                }

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

                contributions[nodeNum] = Math.log(gt) - Math.log(1.0 - p0Parent);
            }
        }

        // External nodes are cheap (just interpolation), no need to parallelize
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            NodeRef node = tree.getExternalNode(i);
            int nodeNum = node.getNumber();
            if (nodeValid[nodeNum]) continue;

            double branchLength = Math.max(tree.getBranchLength(node), h);
            double parentHeight = tree.getNodeHeight(node) + branchLength;

            double sVal = interpolate(S, branchLength);
            double p0Parent = interpolate(p0, parentHeight);

            contributions[nodeNum] = Math.log(sVal) - Math.log(1.0 - p0Parent);
        }
    }

    /**
     * Compute per-node log-likelihood contributions (direct quadrature path).
     * Only computes contributions for nodes where nodeValid[nodeNum] is false.
     */
    private void computeNodeContributionsDirect(double[] contributions) {
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            NodeRef node = tree.getInternalNode(i);
            int nodeNum = node.getNumber();
            if (nodeValid[nodeNum]) continue;

            if (excludeRootBranch && tree.isRoot(node)) {
                contributions[nodeNum] = 0.0;
                continue;
            }

            int k = (int) Math.round(tree.getNodeHeight(node) / h);
            int l;
            if (tree.isRoot(node)) {
                l = numSteps - k;
            } else {
                l = Math.max(1, (int) Math.round(tree.getBranchLength(node) / h));
            }

            contributions[nodeNum] = branchLengthDensityDirect(k, l);
        }

        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            NodeRef node = tree.getExternalNode(i);
            int nodeNum = node.getNumber();
            if (nodeValid[nodeNum]) continue;

            int l = Math.max(1, (int) Math.round(tree.getBranchLength(node) / h));
            contributions[nodeNum] = Math.log(S[l]) - Math.log(1.0 - p0[l]);
        }
    }

    // =========================================================================
    // Direct quadrature solver
    // =========================================================================

    /**
     * Build grid arrays for direct quadrature: map epoch boundaries to grid points,
     * fill rate arrays.
     */
    private void buildGridDirect() {
        for (int i = 0; i <= numSteps; i++) {
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

    private double getEffectiveBirthRate(int i) {
        return birthRateAtGrid[i];
    }

    private double getEffectiveDeathRate(int i) {
        return deathRateAtGrid[i];
    }

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

            // i=m endpoint
            double expmR_m0 = antiDiagBuf[m];
            trap_sum_mu += getDeathRate(m, 0) * expmR_m0 / 2.0;

            trap_sum_mu *= h;
            trap_sum_lam *= h;

            // Implicit quadratic for p0[m]
            double a = h * getBirthRate(m, 0) * expmR_m0 / 2.0;
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
                    (1.0 - h * getBirthRate(i, 0) * antiDiagBuf[i] * p0[i]);
        }

        sComputedUpTo = idx;
    }

    /**
     * Compute branch length density for a node at height k with branch length l (direct quadrature).
     */
    private double branchLengthDensityDirect(int k, int l) {
        double p0_k = p0[k];
        double one_minus_p0_k = 1.0 - p0_k;

        branchDensBuffer[0] = getBirthRate(k, 0) * one_minus_p0_k * one_minus_p0_k;

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
            double denom = 1.0 - h * getBirthRate(i_end, 0) * antiDiagBuf[i_end] * p0[i_end];

            branchDensBuffer[m] = num / denom;
        }

        return Math.log(branchDensBuffer[l]) - Math.log(1.0 - p0[k + l]);
    }

    /**
     * Direct quadrature likelihood computation (non-cached).
     */
    private double calculateLogLikelihoodDirect() {
        cacheRates();

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
     * Single solve at maxSteps resolution, with optional per-node caching
     * for tree-only proposals.
     */
    private double calculateLogLikelihood() {
        double T = times.getParameterValue(numEpochs - 1);
        this.numSteps = maxSteps;
        this.h = T / numSteps;
        ratesDirty = true;

        // Per-node cached path
        boolean needP0S = parametersDirty;
        if (needP0S) {
            invalidateAllNodes();
        }

        cacheRates();

        if (useDirectQuadrature) {
            if (needP0S) {
                gridValid = false;
                buildGridDirect();
                computeP0Direct();
                computeSUpToDirect(numSteps);
            }
            computeNodeContributionsDirect(cachedNodeLogL);
        } else {
            if (needP0S) {
                computeP0AndS();
            }
            computeNodeContributionsFFT(cachedNodeLogL);
        }

        // Mark all recomputed nodes as valid, sum contributions
        double logL = 0.0;
        int totalNodes = tree.getNodeCount();
        for (int i = 0; i < totalNodes; i++) {
            nodeValid[i] = true;
            logL += cachedNodeLogL[i];
        }

        parametersDirty = false;

        if (Double.isNaN(logL) || Double.isInfinite(logL)) {
            return Double.NEGATIVE_INFINITY;
        }
        return logL;
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
        if (model == tree) {
            if (object instanceof TreeChangedEvent) {
                TreeChangedEvent event = (TreeChangedEvent) object;
                if (event.getNode() != null) {
                    invalidateNode(event.getNode());
                } else {
                    invalidateAllNodes();
                }
            } else {
                invalidateAllNodes();
            }
            likelihoodKnown = false;
        } else {
            makeDirty();
        }
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        makeDirty();
    }

    protected void storeState() {
        System.arraycopy(cachedNodeLogL, 0, storedCachedNodeLogL, 0, cachedNodeLogL.length);
        System.arraycopy(nodeValid, 0, storedNodeValid, 0, nodeValid.length);
        System.arraycopy(p0, 0, storedP0, 0, numSteps + 1);
        System.arraycopy(S, 0, storedS, 0, numSteps + 1);
        storedParametersDirty = parametersDirty;
        storedLikelihoodKnown = likelihoodKnown;
        storedLogLikelihood = logLikelihood;
    }

    protected void restoreState() {
        double[] tmpD;
        boolean[] tmpB;

        tmpD = cachedNodeLogL; cachedNodeLogL = storedCachedNodeLogL; storedCachedNodeLogL = tmpD;
        tmpB = nodeValid; nodeValid = storedNodeValid; storedNodeValid = tmpB;
        tmpD = p0;
        System.arraycopy(storedP0, 0, p0, 0, maxSteps + 1);
        System.arraycopy(storedS, 0, S, 0, maxSteps + 1);
        parametersDirty = storedParametersDirty;

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
        parametersDirty = true;
        if (useDirectQuadrature) {
            gridValid = false;
        }
    }

    /**
     * Returns the full p0 grid array only for testing purpose
     */
    public double[] getP0Array() {
        this.numSteps = maxSteps;
        this.h = times.getParameterValue(numEpochs - 1) / numSteps;
        ratesDirty = true;
        if (useDirectQuadrature) {
            gridValid = false;
            cacheRates();
            buildGridDirect();
            computeP0Direct();
        } else {
            computeP0AndS();
        }
        return p0.clone();
    }

    /**
     * Returns the full S grid array only for testing purpose
     */
    public double[] getSArray() {
        this.numSteps = maxSteps;
        this.h = times.getParameterValue(numEpochs - 1) / numSteps;
        ratesDirty = true;
        computeP0AndS();
        return S.clone();
    }

    public String getReport() {
        getLogLikelihood();

        StringBuilder sb = new StringBuilder();
        sb.append("logLikelihood: ").append(logLikelihood).append("\n");

        return sb.toString();
    }

    public String toString() {
        return Double.toString(getLogLikelihood());
    }
}

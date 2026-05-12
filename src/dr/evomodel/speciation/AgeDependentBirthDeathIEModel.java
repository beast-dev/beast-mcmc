package dr.evomodel.speciation;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.tree.TreeChangedEvent;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.FastFourierTransform;
import dr.util.TaskPool;
import dr.xml.Reportable;

import java.util.Arrays;
import java.util.concurrent.atomic.AtomicInteger;

/**
 * @author Frederik M. Andersen
 *
 * Computes the tree likelihood for a skyline, time- and age-dependent birth-death process under
 * symmetric speciation (both daughter lineages reset to age 0 at the speciation event) via
 * integral-equation formulation on a regular time grid.
 *
 * Let p0(t) denote the probability that a lineage born at time t (age 0) leaves no sampled
 * descendants by the present, S(t) the probability that it survives to leave exactly one
 * sampled descendant, and g(t, ell) the density of branches of length ell starting at time t.
 * lambda(t, a) and mu(t, a) are the (skyline, age-dependent) birth and death hazards. With
 * R(t, a) = int_0^a (lambda(t, s) + mu(t, s)) ds the integrated total hazard, the three
 * quantities satisfy Volterra integral equations of the second kind:
 *
 *     p0(t) = int_0^t mu(t-u, u) exp(-R(t-u, u)) du
 *             + int_0^t lambda(t-u, u) exp(-R(t-u, u)) p0(t-u)^2 du
 *
 *     S(t)  = exp(-R(0, t))
 *             + int_0^t lambda(t-u, u) exp(-R(t-u, u)) 2 p0(t-u) S(t-u) du
 *
 *     g(t, ell) = lambda(t, ell) exp(-R(t, ell)) (1 - p0(t))^2
 *             + int_0^ell lambda(t, u) exp(-R(t, u)) 2 p0(t+u) g(t+u, ell - u) du
 *
 * Conditioning:
 *   - on survival from the origin:
 *      S(t_origin) / (1 - p0(t_origin))
 *   - on survival of both mrca subtrees:
 *      prod_v g(t_v, ell_v) / (1 - p0(t_v + ell_v))   (per-internal-node terms)
 *      and    S(ell_tip)    / (1 - p0(parentHeight))  (per-tip terms)
 *
 * Currently implemented for birth and death rates of scaled linear-exponential form:
 *      rate(t, a) = scale(t) (1 + r * gamma * a) exp(-gamma * a)
 * where scale is piecewise constant, r >= 0 and gamma >= 0
 * r = 0 or gamma = 0 collapses to a constant-rate process.
 *
 * Implementation:
 *   - Each Volterra equation is discretized on a uniform grid of step h = originTime / numSteps
 *     and split into current epoch and history convolution equations. These are solved using
 *     Picard iterations, each iteration evaluatied via FFT.
 *   - Node caching ensures that branches not touched by a tree-change operator are not recomputed.
 *   - Internal branch solver can work in parallel.
 *   - An old direct anti-diagonal quadrature solver is retained as an alternative to the FFT-Picard
 *
 * Inputs:
 *   - tree                                   sampled tree
 *   - birthScale, birthShape                 lambda(t, a) parameters
 *   - deathScale, deathShape                 mu(t, a) parameters
 *   - epochTimes                             skyline epoch boundaries
 *   - originTime, numSteps                   integration domain and step count (h = originTime / numSteps)
 *   - epsPicard, maxIterPicard               Picard iteration tolerance and cap
 *   - useDirectQuadrature                    use anti-diagonal direct solver instead of FFT-Picard
 *   - excludeRootBranch                      condition on mrca instead of origin
 *   - numThreads                             worker count for parallel internal-node solves
 *
 * Output:
 *   - logLikelihood                          conditioned tree log-likelihood
 */
public class AgeDependentBirthDeathIEModel extends AbstractModelLikelihood implements Reportable {
    private final Tree tree;
    private final boolean excludeRootBranch;

    private final Parameter birthScale;
    private final Parameter birthShape;
    private final Parameter deathScale;
    private final Parameter deathShape;
    private final Parameter epochTimes;
    private final boolean constBirth;
    private final boolean constDeath;
    private final int numEpochs;
    private final int numBoundaries;

    private final double originTime;
    private final int numSteps;
    private final double h;
    private final double invH;

    private final boolean useDirectQuadrature;
    private final double epsPicard;
    private final int maxIterPicard;

    // Work arrays
    private double[] p0;
    private double[] storedP0;
    private double[] S;
    private double[] storedS;

    private final double[] bScale;
    private final double[] dScale;
    private final double[] bHaz;
    private final double[] dHaz;
    private final double[] bCumHaz;
    private final double[] dCumHaz;
    private final int[] epochIdx;
    private final double[] bScaleDelta;
    private final double[] dScaleDelta;

    // Per-epoch kernel caches
    private final double[][] localExpmR;
    private final double[][] lamKernelTime;
    private final double[][] muKernelTime;
    private final double[][] muKernelCum;

    // Direct quadrature buffers
    private double[] birthRateAtGrid;
    private double[] deathRateAtGrid;
    private double[] antiDiagBuf;
    private double[] branchDensBuffer;
    private int sComputedUpTo;

    // Per-node caching
    private final double[] cachedNodeLogL;
    private final double[] storedCachedNodeLogL;
    private final boolean[] nodeValid;
    private final boolean[] storedNodeValid;
    private final int[] modifiedNodes;
    private int modifiedNodeCount = 0;

    // Compute / dirty flags
    private boolean ratesDirty = true;
    private boolean parametersDirty = true;
    private boolean storedParametersDirty = true;
    private boolean rateStateDirty = false;
    private boolean storedStateDirty = true;
    private boolean likelihoodKnown = false;
    private boolean storedLikelihoodKnown = false;

    private double logLikelihood = 0.0;
    private double storedLogLikelihood = 0.0;

    // Sequential scratch (used by computeP0AndS, single-threaded)
    private final double[] seqP0History;
    private final double[] seqSHistory;
    private final double[] seqRHistory;
    private final double[] seqP0B;
    private final double[] seqSB;
    private final double[] seqGuessCurr;
    private final double[] seqGuessNext;
    private final double[] seqFftA;
    private final double[] seqFftB;

    // Per-worker scratch (used by branchLenDens, parallel)
    private final double[][] workerGt;
    private final double[][] workerGHistory;
    private final double[][] workerRHistory;
    private final double[][] workerGB;
    private final double[][] workerGuessCurr;
    private final double[][] workerGuessNext;
    private final double[][] workerFftA;
    private final double[][] workerFftB;

    // Tip-pass batching
    private final int numExternal;
    private final int[] tipNumsBuf;
    private final double[] tipBranchBuf;

    // Parallelization
    private final int numThreads;
    private final TaskPool taskPool;

    public AgeDependentBirthDeathIEModel(String name,
                                         Tree tree,
                                         Parameter epochTimes,
                                         double originTime,
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
        this.epochTimes = epochTimes;
        this.originTime = originTime;
        this.birthScale = birthScale;
        this.birthShape = birthShape;
        this.deathScale = deathScale;
        this.deathShape = deathShape;

        this.numBoundaries = (epochTimes != null) ? epochTimes.getDimension() : 0;
        this.numEpochs = numBoundaries + 1;
        this.numSteps = numSteps;
        this.h = originTime / numSteps;
        this.invH = 1.0 / h;

        this.epsPicard = epsPicard;
        this.maxIterPicard = maxIterPicard;
        this.useDirectQuadrature = useDirectQuadrature;
        this.excludeRootBranch = excludeRootBranch;
        this.numThreads = Math.max(1, numThreads);
        this.taskPool = new TaskPool(this.numThreads, this.numThreads);

        this.constBirth = birthScale.getDimension() == 1;
        this.constDeath = deathScale.getDimension() == 1;

        if (tree instanceof Model) {
            addModel((Model) tree);
        }

        if (epochTimes != null) {
            addVariable(epochTimes);
        }
        addVariable(birthScale);
        addVariable(birthShape);
        addVariable(deathScale);
        addVariable(deathShape);

        this.p0 = new double[numSteps + 1];
        this.S = new double[numSteps + 1];
        this.storedP0 = new double[numSteps + 1];
        this.storedS = new double[numSteps + 1];

        this.bScale = new double[numEpochs];
        this.dScale = new double[numEpochs];

        this.bHaz = new double[numSteps + 1];
        this.dHaz = new double[numSteps + 1];
        this.bCumHaz = new double[numSteps + 1];
        this.dCumHaz = new double[numSteps + 1];

        this.epochIdx = new int[numEpochs + 1];
        this.bScaleDelta = new double[numEpochs];
        this.dScaleDelta = new double[numEpochs];

        this.localExpmR = new double[numEpochs][numSteps + 1];
        this.lamKernelTime = new double[numEpochs][numSteps + 1];
        this.muKernelTime = new double[numEpochs][numSteps + 1];
        this.muKernelCum = new double[numEpochs][numSteps + 1];

        // FFT pad bound
        int maxFftPad = nextPow2(Math.max(2, 2 * numSteps));

        if (useDirectQuadrature) {
            this.birthRateAtGrid = new double[numSteps + 1];
            this.deathRateAtGrid = new double[numSteps + 1];
            this.antiDiagBuf = new double[numSteps + 1];
            this.branchDensBuffer = new double[numSteps + 1];
        }

        int totalNodes = tree.getNodeCount();
        this.cachedNodeLogL = new double[totalNodes];
        this.storedCachedNodeLogL = new double[totalNodes];
        this.nodeValid = new boolean[totalNodes];
        this.storedNodeValid = new boolean[totalNodes];
        this.modifiedNodes = new int[totalNodes];

        // Sequential scratch for computeP0AndS
        this.seqP0History = new double[numSteps + 1];
        this.seqSHistory = new double[numSteps + 1];
        this.seqRHistory = new double[numSteps + 1];
        this.seqP0B = new double[numSteps + 1];
        this.seqSB = new double[numSteps + 1];
        this.seqGuessCurr = new double[numSteps + 1];
        this.seqGuessNext = new double[numSteps + 1];
        this.seqFftA = new double[maxFftPad + 2];
        this.seqFftB = new double[maxFftPad + 2];

        // Per-worker scratch for branchLenDens
        this.workerGt = new double[this.numThreads][numSteps + 1];
        this.workerGHistory = new double[this.numThreads][numSteps + 1];
        this.workerRHistory = new double[this.numThreads][numSteps + 1];
        this.workerGB = new double[this.numThreads][numSteps + 1];
        this.workerGuessCurr = new double[this.numThreads][numSteps + 1];
        this.workerGuessNext = new double[this.numThreads][numSteps + 1];
        this.workerFftA = new double[this.numThreads][maxFftPad + 2];
        this.workerFftB = new double[this.numThreads][maxFftPad + 2];

        this.numExternal = tree.getExternalNodeCount();
        this.tipNumsBuf = new int[numExternal];
        this.tipBranchBuf = new double[numExternal];
    }


    /**
     * Cache rates and integration kernels
     */
    private void cacheRates() {
        if (!ratesDirty) return;

        epochIdx[0] = 0;
        bScaleDelta[0] = 0.0;
        dScaleDelta[0] = 0.0;
        for (int k = 0; k < numEpochs; k++) {
            dScale[k] = deathScale.getParameterValue(constDeath ? 0 : k);
            bScale[k] = birthScale.getParameterValue(constBirth ? 0 : k);
            epochIdx[k + 1] = (k < numBoundaries) ? closestIdx(epochTimes.getParameterValue(k)) : numSteps;

            if (k == 0) continue;
            bScaleDelta[k] = bScale[k] - bScale[k - 1];
            dScaleDelta[k] = dScale[k] - dScale[k - 1];
        }

        double bR = birthShape.getParameterValue(0);
        double bGamma = birthShape.getParameterValue(1);
        double dR = deathShape.getParameterValue(0);
        double dGamma = deathShape.getParameterValue(1);
        double bB = bR * bGamma;
        double dB = dR * dGamma;

        for (int d = 0; d <= numSteps; d++) {
            double age = d * h;
            bHaz[d] = (1.0 + bB * age) * Math.exp(-bGamma * age);
            dHaz[d] = (1.0 + dB * age) * Math.exp(-dGamma * age);
            bCumHaz[d] = linExpCumHaz(age, bB, bGamma);
            dCumHaz[d] = linExpCumHaz(age, dB, dGamma);
        }

        // Per-epoch integration kernels
        for (int k = 0; k < numEpochs; k++) {
            double bs = bScale[k];
            double ds = dScale[k];
            double[] lemr = localExpmR[k];
            double[] lkt = lamKernelTime[k];
            double[] mkt = muKernelTime[k];
            double[] mkc = muKernelCum[k];
            for (int d = 0; d <= numSteps; d++) {
                double v = Math.exp(-ds * dCumHaz[d] - bs * bCumHaz[d]);
                lemr[d] = v;
                lkt[d] = bs * bHaz[d] * v;
                mkt[d] = ds * dHaz[d] * v;
                mkc[d] = (d == 0) ? mkt[0] : mkc[d - 1] + mkt[d];
            }
        }

        ratesDirty = false;
    }

    private static double linExpCumHaz(double a, double b, double gamma) {
        if (gamma == 0.0) {
            return a + b * a * a / 2.0;
        }
        double emga = Math.exp(-gamma * a);
        double invG = 1.0 / gamma;
        return (invG + b * invG * invG) * (1.0 - emga) - b * a * invG * emga;
    }

    /**
     * Build and transform the epoch-k Picard kernel of length M. Returns padded array length.
     */
    private int buildLamKernelFFT(int k, int M, double[] scratch) {
        int padM = nextPow2(2 * M - 1);
        Arrays.fill(scratch, 0, padM + 2, 0.0);
        double[] lkt = lamKernelTime[k];
        System.arraycopy(lkt, 0, scratch, 0, M);
        double K0 = lkt[0];
        FastFourierTransform.rfft(scratch, padM, false);
        scratch[padM] = K0;
        return padM;
    }

    /**
     * Build and transform the cross-epoch kernel of length L.
     */
    private void buildLamKernelHistory(int k, int L, int pad, double[] scratch) {
        Arrays.fill(scratch, 0, pad + 2, 0.0);
        double[] lkt = lamKernelTime[k];
        System.arraycopy(lkt, 0, scratch, 0, L);
        FastFourierTransform.rfft(scratch, pad, false);
    }

    // ==============================================
    // p0 and S grid solve (FFT-Picard, sequential)
    // ==============================================

    private void computeP0AndS() {
        cacheRates();

        double[] p0History = seqP0History;
        double[] SHistory = seqSHistory;
        double[] RHistory = seqRHistory;
        Arrays.fill(p0History, 0, numSteps + 1, 0.0);
        Arrays.fill(SHistory, 0, numSteps + 1, 0.0);
        Arrays.fill(RHistory, 0, numSteps + 1, 0.0);

        double[] p0B = seqP0B;
        double[] SB = seqSB;
        double[] fBuf = seqFftA;
        double[] kernelFFT = seqFftB;

        for (int k = 0; k < numEpochs; k++) {
            int startIdx = epochIdx[k];
            int endIdx = epochIdx[k + 1];
            int epochLen = endIdx - startIdx;
            if (k == numEpochs - 1) epochLen++;

            // History correction (skipped if rates are constant across epochs)
            if (k > 0 && (!constBirth || !constDeath) && (bScaleDelta[k] != 0.0 || dScaleDelta[k] != 0.0)) {
                for (int g = startIdx; g <= numSteps; g++) {
                    int d = g - startIdx;
                    double deltaR = dScaleDelta[k] * dCumHaz[d] + bScaleDelta[k] * bCumHaz[d];
                    double factor = Math.exp(-deltaR);
                    p0History[g] *= factor;
                    SHistory[g] *= factor;
                    RHistory[g] += deltaR;
                }
            }

            double[] muKernel = muKernelTime[k];
            double[] muCum = muKernelCum[k];

            int picardPad = buildLamKernelFFT(k, epochLen, kernelFFT);
            double K0 = kernelFFT[picardPad];

            // p0 boundary function
            p0B[0] = 0.0;
            for (int i = 1; i < epochLen; i++) {
                p0B[i] = h * (muCum[i] - 0.5 * muKernel[0] - 0.5 * muKernel[i]);
            }
            if (k > 0) {
                for (int i = 0; i < epochLen; i++) p0B[i] += p0History[startIdx + i];
            }

            // Picard iterations for p0
            double[] gc = seqGuessCurr;
            double[] gn = seqGuessNext;
            System.arraycopy(p0B, 0, gc, 0, epochLen);

            boolean converged = false;
            int iter = 0;
            while (!converged && iter < maxIterPicard) {
                Arrays.fill(fBuf, 0, picardPad + 2, 0.0);
                for (int i = 0; i < epochLen; i++) fBuf[i] = gc[i] * gc[i];
                fBuf[0] *= 0.5;
                FastFourierTransform.rfft(fBuf, picardPad, false);
                convolve(kernelFFT, fBuf, picardPad);

                double maxDiff = 0.0;
                for (int i = 0; i < epochLen; i++) {
                    gn[i] = p0B[i] + fBuf[i] * h - 0.5 * h * K0 * gc[i] * gc[i];
                    double diff = Math.abs(gn[i] - gc[i]);
                    if (diff > maxDiff) maxDiff = diff;
                }
                double[] tmp = gc; gc = gn; gn = tmp;
                if (maxDiff < epsPicard) converged = true;
                iter++;
            }
            System.arraycopy(gc, 0, p0, startIdx, epochLen);

            // S boundary function
            for (int i = 0; i < epochLen; i++) {
                int a = startIdx + i;
                SB[i] = Math.exp(-RHistory[a] - dScale[0] * dCumHaz[a] - bScale[0] * bCumHaz[a]);
                if (k > 0) SB[i] += SHistory[startIdx + i];
            }

            // Picard iteration for S
            System.arraycopy(SB, 0, gc, 0, epochLen);

            converged = false;
            iter = 0;
            while (!converged && iter < maxIterPicard) {
                Arrays.fill(fBuf, 0, picardPad + 2, 0.0);
                for (int i = 0; i < epochLen; i++) fBuf[i] = 2.0 * p0[startIdx + i] * gc[i];
                fBuf[0] *= 0.5;
                FastFourierTransform.rfft(fBuf, picardPad, false);
                convolve(kernelFFT, fBuf, picardPad);

                double maxDiff = 0.0;
                for (int i = 0; i < epochLen; i++) {
                    double fi = 2.0 * p0[startIdx + i] * gc[i];
                    gn[i] = SB[i] + fBuf[i] * h - 0.5 * h * K0 * fi;
                    double diff = Math.abs(gn[i] - gc[i]);
                    if (diff > maxDiff) maxDiff = diff;
                }
                double[] tmp = gc; gc = gn; gn = tmp;
                if (maxDiff < epsPicard) converged = true;
                iter++;
            }
            System.arraycopy(gc, 0, S, startIdx, epochLen);

            // History updates for next epochs
            if (k < numEpochs - 1) {
                // Death contribution to p0 history
                for (int g = startIdx; g <= numSteps; g++) {
                    int pos = g - startIdx;
                    int lo = Math.max(0, pos - epochLen + 1);
                    double windowSum = muCum[pos] - (lo > 0 ? muCum[lo - 1] : 0.0);
                    double trapCorr = (startIdx == 0) ? 0.5 * muKernel[pos] : 0.0;
                    p0History[g] += h * windowSum - h * trapCorr;
                }

                int historyPad = nextPow2(numSteps + epochLen);
                buildLamKernelHistory(k, numSteps + 1, historyPad, kernelFFT);

                // p0 birth-history signal
                Arrays.fill(fBuf, 0, historyPad + 2, 0.0);
                for (int idx = 0; idx < epochLen; idx++) {
                    int pastIdx = startIdx + idx;
                    double tw = (pastIdx == 0) ? h * 0.5 : h;
                    fBuf[idx] = p0[pastIdx] * p0[pastIdx] * tw;
                }
                FastFourierTransform.rfft(fBuf, historyPad, false);
                convolve(kernelFFT, fBuf, historyPad);
                for (int pos = 0; pos < historyPad; pos++) {
                    int g = startIdx + pos;
                    if (g > numSteps) break;
                    p0History[g] += fBuf[pos];
                }

                // S birth-history signal
                Arrays.fill(fBuf, 0, historyPad + 2, 0.0);
                for (int idx = 0; idx < epochLen; idx++) {
                    int pastIdx = startIdx + idx;
                    double tw = (pastIdx == 0) ? h * 0.5 : h;
                    fBuf[idx] = 2.0 * p0[pastIdx] * S[pastIdx] * tw;
                }
                FastFourierTransform.rfft(fBuf, historyPad, false);
                convolve(kernelFFT, fBuf, historyPad);
                for (int pos = 0; pos < historyPad; pos++) {
                    int g = startIdx + pos;
                    if (g > numSteps) break;
                    SHistory[g] += fBuf[pos];
                }
            }
        }
    }

    // ===============================================
    // Branch length density (FFT-Picard, parallel)
    // ===============================================

    public double branchLenDens(int worker, double t, double ell) {
        int tIdx = closestIdx(t);
        int maxEll = Math.min((int) Math.ceil(ell * invH) + 1, numSteps - tIdx);

        double[] gt = workerGt[worker];
        double[] gHistory = workerGHistory[worker];
        double[] RHistory = workerRHistory[worker];
        double[] gB = workerGB[worker];
        double[] gc = workerGuessCurr[worker];
        double[] gn = workerGuessNext[worker];
        double[] fBuf = workerFftA[worker];
        double[] kernelFFT = workerFftB[worker];

        Arrays.fill(gt, 0, maxEll + 1, 0.0);
        Arrays.fill(gHistory, 0, maxEll + 1, 0.0);
        Arrays.fill(RHistory, 0, maxEll + 1, 0.0);

        int kStart = 0;
        for (int m = 1; m < numEpochs; m++) {
            if (tIdx >= epochIdx[m]) kStart = m;
        }

        double p0t = interpolate(p0, t);

        for (int k = kStart; k < numEpochs; k++) {
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

            if (k > kStart && (!constBirth || !constDeath) && (bScaleDelta[k] != 0.0 || dScaleDelta[k] != 0.0)) {
                int ellBound = epochIdx[k] - tIdx;
                for (int l = ellBound; l <= maxEll; l++) {
                    int d = Math.min(l - ellBound, numSteps);
                    double deltaR = dScaleDelta[k] * dCumHaz[d] + bScaleDelta[k] * bCumHaz[d];
                    gHistory[l] *= Math.exp(-deltaR);
                    RHistory[l] += deltaR;
                }
            }

            // Boundary function
            double oneMinusP0t = 1.0 - p0t;
            double oneMinusP0t2 = oneMinusP0t * oneMinusP0t;
            double bsK = bScale[kStart];
            double dsK = dScale[kStart];
            for (int i = 0; i < ellLen; i++) {
                int l = ellStart + i;
                gB[i] = bsK * bHaz[l]
                        * Math.exp(-dsK * dCumHaz[l] - bsK * bCumHaz[l] - RHistory[l])
                        * oneMinusP0t2;
                if (k > kStart) gB[i] += gHistory[ellStart + i];
            }

            int picardPad = buildLamKernelFFT(k, ellLen, kernelFFT);
            double gtK0 = kernelFFT[picardPad];

            System.arraycopy(gB, 0, gc, 0, ellLen);

            boolean converged = false;
            int iter = 0;
            while (!converged && iter < maxIterPicard) {
                Arrays.fill(fBuf, 0, picardPad + 2, 0.0);
                for (int i = 0; i < ellLen; i++) {
                    fBuf[i] = 2.0 * p0[tIdx + ellStart + i] * gc[i];
                }
                fBuf[0] *= 0.5;
                FastFourierTransform.rfft(fBuf, picardPad, false);
                convolve(kernelFFT, fBuf, picardPad);

                double maxDiff = 0.0;
                for (int i = 0; i < ellLen; i++) {
                    double fi = 2.0 * p0[tIdx + ellStart + i] * gc[i];
                    gn[i] = gB[i] + fBuf[i] * h - 0.5 * h * gtK0 * fi;
                    double diff = Math.abs(gn[i] - gc[i]);
                    if (diff > maxDiff) maxDiff = diff;
                }
                double[] tmp = gc; gc = gn; gn = tmp;
                if (maxDiff < epsPicard) converged = true;
                iter++;
            }
            System.arraycopy(gc, 0, gt, ellStart, ellLen);

            if (k < numEpochs - 1 && ellEnd < maxEll) {
                int historyLen = ellLen + maxEll - ellEnd;
                int historyPad = nextPow2(historyLen + ellLen - 1);

                buildLamKernelHistory(k, historyLen, historyPad, kernelFFT);

                Arrays.fill(fBuf, 0, historyPad + 2, 0.0);
                for (int idx = 0; idx < ellLen; idx++) {
                    int l = ellStart + idx;
                    double tw = (l == 0) ? h * 0.5 : h;
                    fBuf[idx] = 2.0 * p0[tIdx + l] * gt[l] * tw;
                }
                FastFourierTransform.rfft(fBuf, historyPad, false);
                convolve(kernelFFT, fBuf, historyPad);
                for (int pos = 0; pos < historyPad; pos++) {
                    int l = ellStart + pos;
                    if (l > maxEll) break;
                    gHistory[l] += fBuf[pos];
                }
            }
        }

        // Boundary-aware interpolation: gt is only populated up to index
        // (numSteps - tIdx). Floating-point ell can land between this last
        // populated slot and the unpopulated next slot, so clamp the index.
        int gtUpper = numSteps - tIdx;
        double idx = ell * invH;
        if (idx >= gtUpper) return gt[gtUpper];
        int lo = (int) Math.floor(idx);
        double frac = idx - lo;
        return (1.0 - frac) * gt[lo] + frac * gt[lo + 1];
    }

    // =======================
    // Per-node contributions
    // =======================

    private void computeNodeContributionsFFT() {
        int internalCount = tree.getInternalNodeCount();

        // Collect invalidated internal nodes into a flat list
        int[] dirty = modifiedNodes;
        int dirtyCount = 0;
        for (int i = 0; i < internalCount; i++) {
            NodeRef node = tree.getInternalNode(i);
            int nodeNum = node.getNumber();
            if (nodeValid[nodeNum]) continue;

            if (excludeRootBranch && tree.isRoot(node)) {
                cachedNodeLogL[nodeNum] = 0.0;
                dirty[dirtyCount++] = nodeNum;
                continue;
            }
            dirty[dirtyCount++] = nodeNum;
        }
        modifiedNodeCount = dirtyCount;

        if (dirtyCount > 0) {
            // Indices of nodes that need a Picard solve
            final int[] internalToSolve = new int[dirtyCount];
            final double[] internalHeights = new double[dirtyCount];
            final double[] internalBranchLens = new double[dirtyCount];
            int solveCount = 0;
            for (int i = 0; i < dirtyCount; i++) {
                int nodeNum = dirty[i];
                NodeRef node = tree.getNode(nodeNum);
                if (excludeRootBranch && tree.isRoot(node)) continue;
                double nodeHeight = tree.getNodeHeight(node);
                double bl = tree.isRoot(node) ? originTime - nodeHeight : tree.getBranchLength(node);
                internalToSolve[solveCount] = nodeNum;
                internalHeights[solveCount] = nodeHeight;
                internalBranchLens[solveCount] = Math.max(bl, h);
                solveCount++;
            }
            final int totalSolves = solveCount;

            if (totalSolves == 1 || numThreads == 1) {
                for (int i = 0; i < totalSolves; i++) {
                    int nodeNum = internalToSolve[i];
                    double nodeHeight = internalHeights[i];
                    double bl = internalBranchLens[i];
                    double gt = branchLenDens(0, nodeHeight, bl);
                    double p0Parent = interpolate(p0, nodeHeight + bl);
                    cachedNodeLogL[nodeNum] = Math.log(gt) - Math.log(1.0 - p0Parent);
                }
            } else {
                final AtomicInteger next = new AtomicInteger(0);
                taskPool.fork((task, thread) -> {
                    int idx;
                    while ((idx = next.getAndIncrement()) < totalSolves) {
                        int nodeNum = internalToSolve[idx];
                        double nodeHeight = internalHeights[idx];
                        double bl = internalBranchLens[idx];
                        double gt = branchLenDens(thread, nodeHeight, bl);
                        double p0Parent = interpolate(p0, nodeHeight + bl);
                        cachedNodeLogL[nodeNum] = Math.log(gt) - Math.log(1.0 - p0Parent);
                        }
                });
            }
        }

        // Batch tips by branch length, share the interpolation factors for equal-length tips (cherries / equal-date tips).
        int dirtyTips = 0;
        for (int i = 0; i < numExternal; i++) {
            NodeRef node = tree.getExternalNode(i);
            int nodeNum = node.getNumber();
            if (nodeValid[nodeNum]) continue;
            tipNumsBuf[dirtyTips] = nodeNum;
            tipBranchBuf[dirtyTips] = Math.max(tree.getBranchLength(node), h);
            dirtyTips++;
        }
        if (dirtyTips > 0) {
            sortByKey(tipNumsBuf, tipBranchBuf, dirtyTips);
            double prevBl = -1.0;
            double sVal = 0.0, p0Parent = 0.0, contrib = 0.0;
            for (int i = 0; i < dirtyTips; i++) {
                int nodeNum = tipNumsBuf[i];
                double bl = tipBranchBuf[i];
                if (bl != prevBl) {
                    sVal = interpolate(S, bl);
                    p0Parent = interpolate(p0, tree.getNodeHeight(tree.getNode(nodeNum)) + bl);
                    contrib = Math.log(sVal) - Math.log(1.0 - p0Parent);
                    prevBl = bl;
                } else {
                    // Same branch length, but parentHeight = nodeHeight + bl varies.
                    // p0 depends on parentHeight, so recompute. Only S(bl) is shared.
                    p0Parent = interpolate(p0, tree.getNodeHeight(tree.getNode(nodeNum)) + bl);
                    contrib = Math.log(sVal) - Math.log(1.0 - p0Parent);
                }
                cachedNodeLogL[nodeNum] = contrib;
                modifiedNodes[modifiedNodeCount++] = nodeNum;
            }
        }
    }

    private void computeNodeContributionsDirect() {
        int dirtyCount = 0;
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            NodeRef node = tree.getInternalNode(i);
            int nodeNum = node.getNumber();
            if (nodeValid[nodeNum]) continue;

            if (excludeRootBranch && tree.isRoot(node)) {
                cachedNodeLogL[nodeNum] = 0.0;
                modifiedNodes[dirtyCount++] = nodeNum;
                continue;
            }

            int k = (int) Math.round(tree.getNodeHeight(node) * invH);
            int l;
            if (tree.isRoot(node)) {
                l = numSteps - k;
            } else {
                l = Math.max(1, (int) Math.round(tree.getBranchLength(node) * invH));
            }
            cachedNodeLogL[nodeNum] = branchLengthDensityDirect(k, l);
            modifiedNodes[dirtyCount++] = nodeNum;
        }
        for (int i = 0; i < numExternal; i++) {
            NodeRef node = tree.getExternalNode(i);
            int nodeNum = node.getNumber();
            if (nodeValid[nodeNum]) continue;
            int l = Math.max(1, (int) Math.round(tree.getBranchLength(node) * invH));
            cachedNodeLogL[nodeNum] = Math.log(S[l]) - Math.log(1.0 - p0[l]);
            modifiedNodes[dirtyCount++] = nodeNum;
        }
        modifiedNodeCount = dirtyCount;
    }

    // ======================
    // Direct quadrature TODO: Should probably just be deleted...
    // ======================

    private void buildGridDirect() {
        for (int i = 0; i <= numSteps; i++) {
            int epoch = 0;
            for (int kk = 1; kk < numEpochs; kk++) {
                if (i >= epochIdx[kk]) epoch = kk;
            }
            int birthEpoch = Math.min(epoch, numEpochs - 1);
            int deathEpoch = Math.min(epoch, deathScale.getDimension() - 1);
            birthRateAtGrid[i] = bScale[birthEpoch];
            deathRateAtGrid[i] = dScale[deathEpoch];
        }
    }

    private void computeAntiDiagonal(int C, double[] buf) {
        double logExpmR = 0.0;
        buf[C] = 1.0;
        for (int i = C - 1; i >= 0; i--) {
            int j = C - i;
            double bRate = birthRateAtGrid[i];
            double dRate = deathRateAtGrid[i];
            double deltaBH = bCumHaz[j] - bCumHaz[j - 1];
            double deltaDH = dCumHaz[j] - dCumHaz[j - 1];
            logExpmR -= bRate * deltaBH + dRate * deltaDH;
            buf[i] = Math.exp(logExpmR);
        }
    }

    private double getBirthRate(int i, int j) { return birthRateAtGrid[i] * bHaz[j]; }
    private double getDeathRate(int i, int j) { return deathRateAtGrid[i] * dHaz[j]; }

    private void computeP0Direct() {
        p0[0] = 0.0;
        S[0] = 1.0;
        sComputedUpTo = 0;

        for (int m = 1; m <= numSteps; m++) {
            computeAntiDiagonal(m, antiDiagBuf);
            double tsm = 0.0, tsl = 0.0;
            tsm += getDeathRate(0, m) * antiDiagBuf[0] / 2.0;
            for (int i = 1; i < m; i++) {
                int j = m - i;
                double e = antiDiagBuf[i];
                tsm += getDeathRate(i, j) * e;
                tsl += getBirthRate(i, j) * e * p0[i] * p0[i];
            }
            double expmR_m0 = antiDiagBuf[m];
            tsm += getDeathRate(m, 0) * expmR_m0 / 2.0;
            tsm *= h;
            tsl *= h;

            double a = h * getBirthRate(m, 0) * expmR_m0 / 2.0;
            double c = tsm + tsl;
            if (a == 0.0) {
                p0[m] = c;
            } else {
                double det = 1.0 - 4.0 * a * c;
                p0[m] = det < 0.0 ? 1.0 / (2.0 * a) : (1.0 - Math.sqrt(det)) / (2.0 * a);
            }
        }
    }

    private void computeSUpToDirect(int idx) {
        if (idx <= sComputedUpTo) return;
        for (int i = sComputedUpTo + 1; i <= idx; i++) {
            computeAntiDiagonal(i, antiDiagBuf);
            double expmR_0i = antiDiagBuf[0];
            double trap_sum = 0.0;
            for (int j = 1; j < i; j++) {
                int kk = i - j;
                trap_sum += getBirthRate(j, kk) * antiDiagBuf[j] * p0[j] * S[j];
            }
            S[i] = (2.0 * h * trap_sum + expmR_0i)
                    / (1.0 - h * getBirthRate(i, 0) * antiDiagBuf[i] * p0[i]);
        }
        sComputedUpTo = idx;
    }

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

    // ============
    // Likelihood
    // ============

    private double calculateLogLikelihood() {
        boolean needP0S = parametersDirty;
        if (needP0S) {
            ratesDirty = true;
            invalidateAllNodes();
        }

        cacheRates();

        if (useDirectQuadrature) {
            if (needP0S) {
                buildGridDirect();
                computeP0Direct();
                computeSUpToDirect(numSteps);
                rateStateDirty = true;
            }
            computeNodeContributionsDirect();
        } else {
            if (needP0S) {
                computeP0AndS();
                rateStateDirty = true;
            }
            computeNodeContributionsFFT();
        }

        // Sum contributions and mark all nodes valid
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

    // ========
    // Utils
    // ========

    private static void convolve(double[] kernel, double[] signal, int padM) {
        signal[0] = kernel[0] * signal[0];
        signal[1] = kernel[1] * signal[1];
        for (int i = 1; i < padM / 2; i++) {
            double kr = kernel[2 * i],     ki = kernel[2 * i + 1];
            double sr = signal[2 * i],     si = signal[2 * i + 1];
            signal[2 * i]     = kr * sr - ki * si;
            signal[2 * i + 1] = kr * si + ki * sr;
        }
        FastFourierTransform.rfft(signal, padM, true);
    }

    private static int nextPow2(int n) {
        if (n <= 1) return 2;
        int p = Integer.highestOneBit(n - 1) << 1;
        return Math.max(p, 2);
    }

    private int closestIdx(double t) {
        int idx = (int) Math.round(t * invH);
        if (idx < 0) idx = 0;
        if (idx > numSteps) idx = numSteps;
        return idx;
    }

    private double interpolate(double[] grid, double t) {
        double idx = t * invH;
        int lo = (int) Math.floor(idx);
        if (lo < 0) lo = 0;
        if (lo >= grid.length - 1) return grid[grid.length - 1];
        double frac = idx - lo;
        return (1.0 - frac) * grid[lo] + frac * grid[lo + 1];
    }

    private static void sortByKey(int[] values, double[] keys, int n) {
        for (int k = 1; k < n; k++) {
            int v = values[k];
            double kk = keys[k];
            int m = k - 1;
            while (m >= 0 && keys[m] > kk) {
                values[m + 1] = values[m];
                keys[m + 1] = keys[m];
                m--;
            }
            values[m + 1] = v;
            keys[m + 1] = kk;
        }
    }

    // ==================
    // Node invalidation
    // ==================

    private void invalidateNode(NodeRef node) {
        nodeValid[node.getNumber()] = false;
        for (int i = 0; i < tree.getChildCount(node); i++) {
            nodeValid[tree.getChild(node, i).getNumber()] = false;
        }
    }

    private void invalidateAllNodes() {
        Arrays.fill(nodeValid, false);
    }

    // ======================
    // MCMC state management
    // ======================

    public Model getModel() { return this; }

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
        if (useDirectQuadrature) ;
        invalidateAllNodes();
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
        storedLogLikelihood = logLikelihood;
        storedLikelihoodKnown = likelihoodKnown;
        storedParametersDirty = parametersDirty;

        if (storedStateDirty || rateStateDirty) {
            System.arraycopy(p0, 0, storedP0, 0, numSteps + 1);
            System.arraycopy(S, 0, storedS, 0, numSteps + 1);
            rateStateDirty = false;
        }

        if (storedStateDirty) {
            System.arraycopy(cachedNodeLogL, 0, storedCachedNodeLogL, 0, cachedNodeLogL.length);
            storedStateDirty = false;
        } else {
            for (int i = 0; i < modifiedNodeCount; i++) {
                int n = modifiedNodes[i];
                storedCachedNodeLogL[n] = cachedNodeLogL[n];
            }
        }
        modifiedNodeCount = 0;

        System.arraycopy(nodeValid, 0, storedNodeValid, 0, nodeValid.length);
    }

    protected void restoreState() {
        logLikelihood = storedLogLikelihood;
        likelihoodKnown = storedLikelihoodKnown;
        parametersDirty = storedParametersDirty;

        if (rateStateDirty) {
            // p0/S were rewritten this step; swap back to the snapshot.
            double[] tmp;
            tmp = p0; p0 = storedP0; storedP0 = tmp;
            tmp = S;  S  = storedS;  storedS  = tmp;
            // The kernel caches were also rebuilt for the rejected state; mark
            // for full reinitialization on the next accepted store.
            storedStateDirty = true;
            // Force kernel recompute on next solve since they reflect rejected rates
            ratesDirty = true;
            if (useDirectQuadrature) ;
        }

        if (storedStateDirty) {
            // Full restore on next accept will reinitialize stored snapshot
            System.arraycopy(storedCachedNodeLogL, 0, cachedNodeLogL, 0, cachedNodeLogL.length);
        } else {
            for (int i = 0; i < modifiedNodeCount; i++) {
                int n = modifiedNodes[i];
                cachedNodeLogL[n] = storedCachedNodeLogL[n];
            }
        }
        System.arraycopy(storedNodeValid, 0, nodeValid, 0, nodeValid.length);

        rateStateDirty = false;
        modifiedNodeCount = 0;
    }

    protected void acceptState() { }

    public String getReport() {
        getLogLikelihood();
        return "logLikelihood: " + logLikelihood
                + " (numSteps=" + numSteps
                + ", numThreads=" + numThreads
                + ", p0Origin=" + p0[numSteps] + ")\n";
    }

    public String toString() {
        return Double.toString(getLogLikelihood());
    }
}

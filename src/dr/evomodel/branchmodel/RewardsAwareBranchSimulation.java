/*
 *
 * Copyright (c) 2002-2024 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

package dr.evomodel.branchmodel;

import dr.math.Poisson;
import dr.math.distributions.PoissonDistribution;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

/**
 * Efficient reward-conditioned CTMC branch path simulator.
 *
 * This class is a Java translation of the uniformization + SIR + constrained
 * dwell-time sampler in {@code paper/assets/dependent-ctmc/simul.py}.  Expensive
 * endpoint bridge quantities are cached for a fixed {@code (Q, rewardRates, xT, T)}
 * so repeated samples with different starts or reward proportions can reuse the
 * same endpoint messages and reward feasibility bounds.
 *
 * The returned path is represented by merged CTMC states and dwell times.  Virtual
 * self-transitions introduced by uniformization are removed.
 *
 * @author Filippo Monti
 */
public class RewardsAwareBranchSimulation {

    public static final double REWARD_TOL = 1.0e-10;
    public static final double AFFINE_TOL = 1.0e-8;
    private static final double LOG_ZERO = Double.NEGATIVE_INFINITY;
    private static final int DEFAULT_EXACT_TRIANGULATION_MAX_COMPONENTS = 12;

    public enum Component {
        CONTINUOUS,
        ATOM,
        AUTO
    }

    public enum DwellSampler {
        AUTO,
        HIT_RUN,
        EXACT_TWO_REWARD,
        EXACT_TRIANGULATION
    }

    public enum ResamplingMethod {
        SYSTEMATIC,
        STRATIFIED,
        MULTINOMIAL
    }

    public static final class Options {
        public int nParticles = 200;
        public int nHitRun = 300;
        public int nBurnin = 100;
        public int nDwellSamples = 1;
        public int hitRunThin = 1;
        public int maxAdaptiveHitRunBurnin = -1;
        public int maxParticleMultiplier = 4;
        public int exactTriangulationMaxComponents = DEFAULT_EXACT_TRIANGULATION_MAX_COMPONENTS;

        public double essThreshold = 0.1;
        public double rewardDensityTemperature = 0.1;

        public boolean rewardAwareProposal = true;
        public boolean rewardDensityAwareProposal = false;
        public boolean adaptiveHitRunBurnin = false;
        public boolean autoTuneSir = true;
        public boolean validatePath = true;

        public Component component = Component.CONTINUOUS;
        public DwellSampler dwellSampler = DwellSampler.AUTO;
        public ResamplingMethod resamplingMethod = ResamplingMethod.SYSTEMATIC;

        public Options copy() {
            Options copy = new Options();
            copy.nParticles = nParticles;
            copy.nHitRun = nHitRun;
            copy.nBurnin = nBurnin;
            copy.nDwellSamples = nDwellSamples;
            copy.hitRunThin = hitRunThin;
            copy.maxAdaptiveHitRunBurnin = maxAdaptiveHitRunBurnin;
            copy.maxParticleMultiplier = maxParticleMultiplier;
            copy.exactTriangulationMaxComponents = exactTriangulationMaxComponents;
            copy.essThreshold = essThreshold;
            copy.rewardDensityTemperature = rewardDensityTemperature;
            copy.rewardAwareProposal = rewardAwareProposal;
            copy.rewardDensityAwareProposal = rewardDensityAwareProposal;
            copy.adaptiveHitRunBurnin = adaptiveHitRunBurnin;
            copy.autoTuneSir = autoTuneSir;
            copy.validatePath = validatePath;
            copy.component = component;
            copy.dwellSampler = dwellSampler;
            copy.resamplingMethod = resamplingMethod;
            return copy;
        }
    }

    public static final class PathSample {
        private final int[] states;
        private final double[] dwellTimes;
        private final double ess;
        private final Diagnostics diagnostics;

        private PathSample(int[] states, double[] dwellTimes, double ess, Diagnostics diagnostics) {
            this.states = states;
            this.dwellTimes = dwellTimes;
            this.ess = ess;
            this.diagnostics = diagnostics;
        }

        public int[] getStates() {
            return states.clone();
        }

        public double[] getDwellTimes() {
            return dwellTimes.clone();
        }

        public double getEffectiveSampleSizeFraction() {
            return ess;
        }

        public Diagnostics getDiagnostics() {
            return diagnostics;
        }

        public double getTotalTime() {
            double total = 0.0;
            for (double dwellTime : dwellTimes) {
                total += dwellTime;
            }
            return total;
        }

        public double getRewardProportion(double[] rewardRates) {
            if (rewardRates == null) {
                throw new IllegalArgumentException("rewardRates must be non-null");
            }
            if (states.length != dwellTimes.length) {
                throw new IllegalStateException("states and dwellTimes have different lengths");
            }
            double totalTime = getTotalTime();
            double reward = 0.0;
            for (int i = 0; i < states.length; i++) {
                reward += rewardRates[states[i]] * dwellTimes[i];
            }
            return reward / totalTime;
        }
    }

    public static final class Diagnostics {
        public String component;
        public boolean cachedPrecomputation;
        public boolean rewardAwareProposalRequested;
        public boolean rewardDensityAwareProposalRequested;
        public boolean rewardAwareProposalUsed;
        public boolean rewardDensityAwareProposalUsed;
        public double rewardDensityTemperature;
        public TruncationInfo truncationInfo;
        public SirDiagnostics sir;
        public DwellDiagnostics dwell;
        public int selectedParticle = -1;
        public int selectedN = -1;
        public int[] selectedSkeleton;
        public String atomType;
        public int[] rewardClassStates;
    }

    public static final class TruncationInfo {
        public final int nMax;
        public final double endpointMass;
        public final double tailBound;
        public final double relativeTailBound;

        private TruncationInfo(int nMax, double endpointMass, double tailBound) {
            this.nMax = nMax;
            this.endpointMass = endpointMass;
            this.tailBound = tailBound;
            this.relativeTailBound = endpointMass > 0.0 ? tailBound / endpointMass : Double.POSITIVE_INFINITY;
        }
    }

    public static final class SirDiagnostics {
        public double ess;
        public double essRaw;
        public double essThreshold;
        public int nParticles;
        public int nFiniteParticles;
        public double maxLogWeight;
        public boolean autoTuned;
        public final List<SirAttempt> attempts = new ArrayList<SirAttempt>();
        public ResamplingInfo resamplingInfo;
    }

    public static final class SirAttempt {
        public final double ess;
        public final int nParticles;
        public final int nFiniteParticles;
        public final boolean rewardAware;
        public final boolean rewardDensityAware;
        public final int batchSize;

        private SirAttempt(double ess, int nParticles, int nFiniteParticles,
                           boolean rewardAware, boolean rewardDensityAware, int batchSize) {
            this.ess = ess;
            this.nParticles = nParticles;
            this.nFiniteParticles = nFiniteParticles;
            this.rewardAware = rewardAware;
            this.rewardDensityAware = rewardDensityAware;
            this.batchSize = batchSize;
        }
    }

    public static final class ResamplingInfo {
        public final ResamplingMethod method;
        public final int chosenSlot;
        public final int uniqueSelectedParticles;

        private ResamplingInfo(ResamplingMethod method, int chosenSlot, int uniqueSelectedParticles) {
            this.method = method;
            this.chosenSlot = chosenSlot;
            this.uniqueSelectedParticles = uniqueSelectedParticles;
        }
    }

    public static final class DwellDiagnostics {
        public String sampler;
        public double constraintError;
        public double sumError;
        public double rewardError;
        public double maxRetainedConstraintError;
        public int skippedMoves;
        public int boundaryHits;
        public int totalMoves;
        public int acceptedMoves;
        public int effectiveBurnin;
        public int thin;
        public int nRetained;
        public int nVertices;
        public double initialMinComponent;
        public double[] selectedSample;
        public double[][] retainedSamples;
        public int retry;
    }

    private final double[][] qMatrix;
    private final double[] rewardRates;
    private final int stateCount;
    private final int xT;
    private final double branchLength;
    private final double epsTruncation;
    private final int maxTruncationN;

    private final double omega;
    private final double omegaT;
    private final double[][] transitionKernel;
    private final int[] nMaxByStart;
    private final double[][] endpointMessages;
    private final TruncationInfo[] truncationInfoByStart;
    private final boolean deterministicReward;

    private final double minReward;
    private final double maxReward;
    private final SkeletonScratch skeletonScratch;
    private final RewardDensityScratch rewardDensityScratch;

    private FutureBounds futureBounds;
    private FutureMoments futureMoments;

    public RewardsAwareBranchSimulation(double[][] qMatrix, double[] rewardRates,
                                        int xT, double branchLength) {
        this(qMatrix, rewardRates, xT, branchLength, 1.0e-10, -1);
    }

    public RewardsAwareBranchSimulation(double[][] qMatrix, double[] rewardRates,
                                        int xT, double branchLength,
                                        double epsTruncation,
                                        int maxTruncationN) {
        validateStaticInputs(qMatrix, rewardRates, xT, branchLength, epsTruncation);
        this.qMatrix = copyMatrix(qMatrix);
        this.rewardRates = rewardRates.clone();
        this.stateCount = qMatrix.length;
        this.xT = xT;
        this.branchLength = branchLength;
        this.epsTruncation = epsTruncation;
        this.maxTruncationN = maxTruncationN;

        Uniformization uniformization = uniformize(this.qMatrix);
        this.omega = uniformization.omega;
        this.transitionKernel = uniformization.kernel;
        this.omegaT = this.omega * this.branchLength;

        EndpointTruncation truncation = endpointAwareTruncationAllStarts(
                this.transitionKernel, this.xT, this.omegaT, this.epsTruncation, this.maxTruncationN);
        this.nMaxByStart = truncation.nMaxByStart;
        this.endpointMessages = truncation.endpointMessages;
        this.truncationInfoByStart = truncation.infoByStart;
        this.minReward = min(this.rewardRates);
        this.maxReward = max(this.rewardRates);
        this.deterministicReward = this.maxReward - this.minReward <= REWARD_TOL;
        this.skeletonScratch = new SkeletonScratch(this.stateCount);
        this.rewardDensityScratch = new RewardDensityScratch();
    }

    public double getUniformizationRate() {
        return omega;
    }

    public double[][] getTransitionKernel() {
        return copyMatrix(transitionKernel);
    }

    public int getCachedNMax() {
        return endpointMessages.length - 1;
    }

    public int getNMaxForStart(int x0) {
        checkStateIndex(x0, "x0");
        return nMaxByStart[x0];
    }

    public TruncationInfo getTruncationInfoForStart(int x0) {
        checkStateIndex(x0, "x0");
        return truncationInfoByStart[x0];
    }

    public PathSample sample(int x0, double rho, Random random) {
        return sample(x0, rho, new Options(), random);
    }

    public PathSample sample(int x0, double rho, Options options, Random random) {
        if (options == null) {
            options = new Options();
        } else {
            options = options.copy();
        }
        if (random == null) {
            random = new Random();
        }
        validateSamplingInputs(x0, rho, options);

        Diagnostics diagnostics = new Diagnostics();
        diagnostics.cachedPrecomputation = true;
        diagnostics.component = options.component.name().toLowerCase();
        diagnostics.rewardAwareProposalRequested = options.rewardAwareProposal;
        diagnostics.rewardDensityAwareProposalRequested = options.rewardDensityAwareProposal;
        diagnostics.rewardDensityTemperature = options.rewardDensityTemperature;

        boolean atomCompatible = rewardAtomCompatible(qMatrix, rewardRates, x0, xT, rho);
        boolean rhoOnBoundary = isRewardBoundary(rho);

        if (options.component == Component.ATOM) {
            if (!atomCompatible) {
                throw new IllegalArgumentException("The requested reward-class atom is incompatible with x0, xT, and rho.");
            }
            return sampleRewardAtomPath(x0, rho, random, "atom", diagnostics, options.validatePath);
        }

        if (deterministicReward) {
            if (Math.abs(rho - rewardRates[0]) > REWARD_TOL) {
                throw new IllegalArgumentException("All reward rates are equal, so the requested reward is impossible.");
            }
            return sampleRewardAtomPath(x0, rho, random, "deterministic_reward", diagnostics, options.validatePath);
        }

        if (options.component == Component.AUTO && atomCompatible && rhoOnBoundary) {
            return sampleRewardAtomPath(x0, rho, random, "atom", diagnostics, options.validatePath);
        }

        int nMax = nMaxByStart[x0];
        if (nMax < 0) {
            throw new IllegalArgumentException("No valid paths from x0=" + x0 + " to xT=" + xT + ".");
        }
        diagnostics.truncationInfo = truncationInfoByStart[x0];

        double[] qN = endpointConditionedCountDistribution(x0, nMax, endpointMessages, omegaT);
        FutureBounds activeFutureBounds = null;
        FutureMoments activeFutureMoments = null;
        boolean activeRewardAware = options.rewardAwareProposal;
        if (activeRewardAware) {
            activeFutureBounds = getFutureRewardBounds(nMax);
            if (options.rewardDensityAwareProposal) {
                activeFutureMoments = getFutureRewardMoments(nMax);
            }
        }

        int currentBatchSize = options.nParticles;
        int maxParticles = options.nParticles * options.maxParticleMultiplier;
        ParticleStore particles = new ParticleStore(maxParticles);
        WeightNormalization normalized = null;
        SirDiagnostics sirDiagnostics = new SirDiagnostics();
        sirDiagnostics.essThreshold = options.essThreshold;
        sirDiagnostics.autoTuned = options.autoTuneSir;

        while (true) {
            drawParticleBatch(x0, rho, qN, nMax, currentBatchSize,
                    activeRewardAware, activeFutureBounds, activeFutureMoments,
                    options, random, particles);

            try {
                normalized = normalizeLogWeights(particles);
                sirDiagnostics.attempts.add(new SirAttempt(
                        normalized.essFraction,
                        particles.size(),
                        normalized.nFiniteParticles,
                        activeRewardAware,
                        activeRewardAware && options.rewardDensityAwareProposal,
                        currentBatchSize));
            } catch (IllegalArgumentException e) {
                normalized = null;
                sirDiagnostics.attempts.add(new SirAttempt(
                        0.0,
                        particles.size(),
                        0,
                        activeRewardAware,
                        activeRewardAware && options.rewardDensityAwareProposal,
                        currentBatchSize));
            }

            double ess = normalized == null ? 0.0 : normalized.essFraction;
            if (!options.autoTuneSir || ess >= options.essThreshold) {
                break;
            }

            if (!activeRewardAware) {
                activeRewardAware = true;
                activeFutureBounds = getFutureRewardBounds(nMax);
                if (options.rewardDensityAwareProposal) {
                    activeFutureMoments = getFutureRewardMoments(nMax);
                }
                currentBatchSize = options.nParticles;
                continue;
            }

            int nextTotal = Math.min(particles.size() * 2, maxParticles);
            if (nextTotal > particles.size()) {
                currentBatchSize = nextTotal - particles.size();
                continue;
            }
            break;
        }

        if (normalized == null) {
            throw new IllegalArgumentException("All SIR weights are zero. The conditioning event may be incompatible or very unlikely.");
        }

        sirDiagnostics.ess = normalized.essFraction;
        sirDiagnostics.essRaw = normalized.essRaw;
        sirDiagnostics.nParticles = particles.size();
        sirDiagnostics.nFiniteParticles = normalized.nFiniteParticles;
        sirDiagnostics.maxLogWeight = normalized.maxLogWeight;
        diagnostics.sir = sirDiagnostics;
        diagnostics.rewardAwareProposalUsed = activeRewardAware;
        diagnostics.rewardDensityAwareProposalUsed = activeRewardAware && options.rewardDensityAwareProposal;

        SelectedParticle selected = selectResampledParticle(normalized.weights, random, options.resamplingMethod);
        sirDiagnostics.resamplingInfo = selected.resamplingInfo;
        diagnostics.selectedParticle = selected.index;

        int[] chosenSkeleton = particles.getSkeleton(selected.index);
        if (chosenSkeleton == null) {
            throw new IllegalStateException("Selected a null skeleton despite finite SIR weight.");
        }
        diagnostics.selectedSkeleton = chosenSkeleton.clone();
        diagnostics.selectedN = chosenSkeleton.length - 1;

        double[] dwellTimes;
        DwellDiagnostics selectedDwellDiagnostics = null;
        if (chosenSkeleton.length == 1) {
            dwellTimes = new double[]{branchLength};
        } else {
            double[] cChosen = rewardsForSkeleton(chosenSkeleton, rewardRates);
            ReconstructedPath reconstructed = null;
            for (int retry = 0; retry < 50 && reconstructed == null; retry++) {
                DwellResult dwellResult = sampleRewardConstrainedDwell(
                        cChosen, rho, options, random, true);
                for (int i = dwellResult.retainedSamples.length - 1; i >= 0; i--) {
                    double[] spacings = dwellResult.retainedSamples[i];
                    double[] candidateDwellTimes = scale(spacings, branchLength);
                    ReconstructedPath candidate = reconstructCtmcPath(chosenSkeleton, candidateDwellTimes);
                    if (candidate.states.length > 0 &&
                            candidate.states[0] == x0 &&
                            candidate.states[candidate.states.length - 1] == xT) {
                        dwellResult.diagnostics.selectedSample = spacings.clone();
                        dwellResult.diagnostics.retry = retry;
                        selectedDwellDiagnostics = dwellResult.diagnostics;
                        reconstructed = candidate;
                        break;
                    }
                }
            }
            if (reconstructed == null) {
                throw new IllegalStateException("Could not obtain positive endpoint dwell times for the selected skeleton.");
            }
            diagnostics.dwell = selectedDwellDiagnostics;
            if (options.validatePath) {
                validateConditionedPath(reconstructed.states, reconstructed.times, rewardRates, x0, xT, branchLength, rho);
            }
            return new PathSample(reconstructed.states, reconstructed.times, normalized.essFraction, diagnostics);
        }

        ReconstructedPath reconstructed = reconstructCtmcPath(chosenSkeleton, dwellTimes);
        if (options.validatePath) {
            validateConditionedPath(reconstructed.states, reconstructed.times, rewardRates, x0, xT, branchLength, rho);
        }
        return new PathSample(reconstructed.states, reconstructed.times, normalized.essFraction, diagnostics);
    }

    public static PathSample simulateConditional(double[][] qMatrix, double[] rewardRates,
                                                int x0, int xT, double branchLength,
                                                double rho, Options options, Random random) {
        RewardsAwareBranchSimulation simulation =
                new RewardsAwareBranchSimulation(qMatrix, rewardRates, xT, branchLength);
        return simulation.sample(x0, rho, options, random);
    }

    public static double dividedDifferenceDensity(double[] rewardRates, double rho) {
        double logDensity = dividedDifferenceLogDensity(rewardRates, rho);
        return Double.isFinite(logDensity) ? Math.exp(logDensity) : 0.0;
    }

    public static double dividedDifferenceLogDensity(double[] rewardRates, double rho) {
        return dividedDifferenceLogDensity(rewardRates, rho, new RewardDensityScratch());
    }

    private static double dividedDifferenceLogDensity(double[] rewardRates, double rho,
                                                      RewardDensityScratch scratch) {
        scratch.ensureKnotCount(rewardRates.length);
        System.arraycopy(rewardRates, 0, scratch.sortedKnots, 0, rewardRates.length);
        Arrays.sort(scratch.sortedKnots, 0, rewardRates.length);
        return bsplineRewardLogDensityFromSortedKnots(scratch.sortedKnots, rewardRates.length, rho, 1.0e-12, scratch);
    }

    private static double dividedDifferenceLogDensity(int[] skeleton, double[] rewardRates, double rho,
                                                      RewardDensityScratch scratch) {
        scratch.ensureKnotCount(skeleton.length);
        for (int i = 0; i < skeleton.length; i++) {
            scratch.sortedKnots[i] = rewardRates[skeleton[i]];
        }
        Arrays.sort(scratch.sortedKnots, 0, skeleton.length);
        return bsplineRewardLogDensityFromSortedKnots(scratch.sortedKnots, skeleton.length, rho, 1.0e-12, scratch);
    }

    private void drawParticleBatch(int x0, double rho, double[] qN, int nMax, int batchSize,
                                   boolean activeRewardAware, FutureBounds activeFutureBounds,
                                   FutureMoments activeFutureMoments, Options options, Random random,
                                   ParticleStore particles) {
        for (int i = 0; i < batchSize; i++) {
            int n = randomChoice(qN, random);
            if (n > nMax) {
                throw new IllegalStateException("Sampled count outside truncation range.");
            }

            if (n == 0) {
                if (x0 == xT) {
                    int[] skeleton = new int[]{x0};
                    double logWeight = dividedDifferenceLogDensity(skeleton, rewardRates, rho, rewardDensityScratch);
                    particles.add(skeleton, logWeight);
                } else {
                    particles.add(null, LOG_ZERO);
                }
                continue;
            }

            try {
                SkeletonResult skeletonResult = sampleBridgeSkeleton(
                        n, x0, xT, transitionKernel, endpointMessages, random,
                        rewardRates, rho, activeRewardAware,
                        activeFutureBounds, activeFutureMoments,
                        activeRewardAware && options.rewardDensityAwareProposal,
                        options.rewardDensityTemperature,
                        skeletonScratch, minReward, maxReward);

                double logDensity = dividedDifferenceLogDensity(skeletonResult.skeleton, rewardRates, rho, rewardDensityScratch);
                double logWeight = LOG_ZERO;
                if (Double.isFinite(logDensity)) {
                    logWeight = logDensity + skeletonResult.logBridgeCorrection;
                }
                particles.add(skeletonResult.skeleton, logWeight);
            } catch (IllegalArgumentException e) {
                particles.add(null, LOG_ZERO);
            }
        }
    }

    private PathSample sampleRewardAtomPath(int x0, double rho, Random random,
                                           String componentLabel, Diagnostics diagnostics,
                                           boolean validatePath) {
        boolean[] atomMask = rewardAtomStateMask(rewardRates, rho, REWARD_TOL);
        if (!atomMask[x0] || !atomMask[xT]) {
            throw new IllegalArgumentException("The requested atom is incompatible with x0, xT, and rho.");
        }

        int[] atomStates = trueIndices(atomMask);
        int localX0 = indexOf(atomStates, x0);
        int localXT = indexOf(atomStates, xT);
        double[][] qAtom = subMatrix(qMatrix, atomStates);
        double atomOmega = maxAbsDiagonal(qAtom);

        diagnostics.component = componentLabel;
        diagnostics.atomType = "reward_class";
        diagnostics.rewardClassStates = atomStates.clone();
        diagnostics.rewardAwareProposalUsed = false;
        diagnostics.rewardDensityAwareProposalUsed = false;

        if (atomOmega == 0.0) {
            if (localX0 != localXT) {
                throw new IllegalArgumentException("No atomic path connects x0 to xT in a zero-rate chain.");
            }
            int[] states = new int[]{x0};
            double[] times = new double[]{branchLength};
            if (validatePath) {
                validateConditionedPath(states, times, rewardRates, x0, xT, branchLength, rho);
            }
            diagnostics.selectedN = 0;
            diagnostics.selectedSkeleton = states.clone();
            return new PathSample(states, times, 1.0, diagnostics);
        }

        Uniformization atomUniformization = uniformize(qAtom);
        double[][] pAtom = atomUniformization.kernel;
        clampTransitionKernel(pAtom);
        double atomOmegaT = atomOmega * branchLength;
        EndpointTruncation truncation = endpointAwareTruncation(
                pAtom, localX0, localXT, atomOmegaT, epsTruncation, maxTruncationN);
        diagnostics.truncationInfo = truncation.infoByStart[localX0];

        int nMax = truncation.endpointMessages.length - 1;
        double[] qN = endpointConditionedCountDistribution(localX0, nMax, truncation.endpointMessages, atomOmegaT);
        int n = randomChoice(qN, random);

        int[] localSkeleton;
        if (n == 0) {
            localSkeleton = new int[]{localX0};
        } else {
            localSkeleton = sampleBridgeSkeleton(
                    n, localX0, localXT, pAtom, truncation.endpointMessages, random,
                    null, Double.NaN, false, null, null, false, 1.0,
                    new SkeletonScratch(pAtom.length), Double.NaN, Double.NaN).skeleton;
        }

        int[] skeleton = new int[localSkeleton.length];
        for (int i = 0; i < localSkeleton.length; i++) {
            skeleton[i] = atomStates[localSkeleton[i]];
        }
        double[] dwellTimes = n == 0 ? new double[]{branchLength} : scale(sampleDirichlet(skeleton.length, random), branchLength);
        ReconstructedPath path = reconstructCtmcPath(skeleton, dwellTimes);

        if (validatePath) {
            validateConditionedPath(path.states, path.times, rewardRates, x0, xT, branchLength, rho);
        }

        diagnostics.selectedN = n;
        diagnostics.selectedSkeleton = skeleton.clone();
        return new PathSample(path.states, path.times, 1.0, diagnostics);
    }

    private FutureBounds getFutureRewardBounds(int nMax) {
        if (deterministicReward) {
            return null;
        }
        if (futureBounds == null) {
            futureBounds = precomputeFutureRewardBounds(transitionKernel, endpointMessages, rewardRates);
        }
        if (nMax == endpointMessages.length - 1) {
            return futureBounds;
        }
        return new FutureBounds(sliceRows(futureBounds.min, nMax + 1), sliceRows(futureBounds.max, nMax + 1));
    }

    private FutureMoments getFutureRewardMoments(int nMax) {
        if (deterministicReward) {
            return null;
        }
        if (futureMoments == null) {
            futureMoments = precomputeFutureRewardMoments(transitionKernel, endpointMessages, rewardRates);
        }
        if (nMax == endpointMessages.length - 1) {
            return futureMoments;
        }
        return new FutureMoments(sliceRows(futureMoments.mean, nMax + 1), sliceRows(futureMoments.variance, nMax + 1));
    }

    private static SkeletonResult sampleBridgeSkeleton(int n, int x0, int xT,
                                                       double[][] p, double[][] endpointMessages,
                                                       Random random, double[] rewardRates, double rho,
                                                       boolean rewardAware, FutureBounds futureBounds,
                                                       FutureMoments futureMoments,
                                                       boolean rewardDensityAware,
                                                       double rewardDensityTemperature,
                                                       SkeletonScratch scratch,
                                                       double globalMin, double globalMax) {
        if (n < 0) {
            throw new IllegalArgumentException("N must be non-negative.");
        }
        if (n == 0) {
            if (x0 != xT) {
                throw new IllegalArgumentException("A zero-step skeleton requires x0 == xT.");
            }
            return new SkeletonResult(new int[]{x0}, 0.0);
        }
        if (n == 1) {
            if (p[x0][xT] <= 0.0) {
                throw new IllegalArgumentException("No one-step bridge from state " + x0 + " to " + xT + ".");
            }
            return new SkeletonResult(new int[]{x0, xT}, 0.0);
        }

        boolean useRewardFilter = rewardAware && rewardRates != null && Double.isFinite(rho);
        if (rewardDensityTemperature <= 0.0) {
            throw new IllegalArgumentException("rewardDensityTemperature must be positive.");
        }

        int stateCount = p.length;
        if (scratch == null) {
            scratch = new SkeletonScratch(stateCount);
        } else {
            scratch.ensureStateCount(stateCount);
        }
        int[] skeleton = new int[n + 1];
        skeleton[0] = x0;
        double logBridgeCorrection = 0.0;

        double knownMin = 0.0;
        double knownMax = 0.0;
        double knownRewardSum = 0.0;
        double rewardSpan = 0.0;
        if (useRewardFilter) {
            knownMin = rewardRates[x0];
            knownMax = rewardRates[x0];
            knownRewardSum = rewardRates[x0];
            rewardSpan = Math.max(globalMax - globalMin, REWARD_TOL);
        }

        for (int k = 1; k <= n; k++) {
            int prev = skeleton[k - 1];
            int stepsLeft = n - k;
            double[] weights = scratch.weights;
            double total = 0.0;
            for (int state = 0; state < stateCount; state++) {
                double weight = p[prev][state] * endpointMessages[stepsLeft][state];
                weights[state] = weight;
                total += weight;
            }
            if (total <= 0.0) {
                throw new IllegalArgumentException("Zero bridge weight at step " + k + ".");
            }

            double[] baseProbs = scratch.baseProbs;
            for (int state = 0; state < stateCount; state++) {
                baseProbs[state] = weights[state] / total;
            }

            double[] proposalProbs = baseProbs;
            if (useRewardFilter) {
                double[] feasibleMin = scratch.feasibleMin;
                double[] feasibleMax = scratch.feasibleMax;
                boolean[] mask = scratch.mask;

                boolean anyFeasible = false;
                for (int state = 0; state < stateCount; state++) {
                    double candidateMin = Math.min(knownMin, rewardRates[state]);
                    double candidateMax = Math.max(knownMax, rewardRates[state]);
                    if (futureBounds != null) {
                        feasibleMin[state] = Math.min(candidateMin, futureBounds.min[stepsLeft][state]);
                        feasibleMax[state] = Math.max(candidateMax, futureBounds.max[stepsLeft][state]);
                    } else if (stepsLeft > 0) {
                        feasibleMin[state] = Math.min(candidateMin, globalMin);
                        feasibleMax[state] = Math.max(candidateMax, globalMax);
                    } else {
                        feasibleMin[state] = candidateMin;
                        feasibleMax[state] = candidateMax;
                    }
                    mask[state] = baseProbs[state] > 0.0 &&
                            rho >= feasibleMin[state] - REWARD_TOL &&
                            rho <= feasibleMax[state] + REWARD_TOL;
                    anyFeasible |= mask[state];
                }

                if (!anyFeasible) {
                    throw new IllegalArgumentException("Reward-aware bridge proposal found no feasible continuation.");
                }
                proposalProbs = rewardAwareProposal(baseProbs, mask, rho, knownRewardSum, rewardRates,
                        stepsLeft, n, feasibleMin, feasibleMax, rewardDensityAware,
                        rewardDensityTemperature, rewardSpan, futureMoments,
                        scratch.proposalProbs);
            }

            int state = randomChoice(proposalProbs, random);
            skeleton[k] = state;

            if (useRewardFilter) {
                double pBase = baseProbs[state];
                double pProposal = proposalProbs[state];
                if (pBase <= 0.0 || pProposal <= 0.0) {
                    throw new IllegalStateException("Invalid bridge proposal probability.");
                }
                logBridgeCorrection += Math.log(pBase) - Math.log(pProposal);
                double reward = rewardRates[state];
                knownMin = Math.min(knownMin, reward);
                knownMax = Math.max(knownMax, reward);
                knownRewardSum += reward;
            }
        }

        if (skeleton[skeleton.length - 1] != xT) {
            throw new IllegalStateException("Bridge skeleton ended at " + skeleton[skeleton.length - 1] + ", expected " + xT + ".");
        }
        return new SkeletonResult(skeleton, logBridgeCorrection);
    }

    private static double[] rewardAwareProposal(double[] baseProbs, boolean[] mask, double rho,
                                                double knownRewardSum, double[] rewardRates,
                                                int stepsLeft, int n, double[] feasibleMin,
                                                double[] feasibleMax, boolean rewardDensityAware,
                                                double rewardDensityTemperature,
                                                double rewardSpan, FutureMoments futureMoments,
                                                double[] proposal) {
        Arrays.fill(proposal, 0, baseProbs.length, 0.0);
        double proposalMass = 0.0;

        if (rewardDensityAware) {
            for (int state = 0; state < baseProbs.length; state++) {
                if (!mask[state]) {
                    continue;
                }
                double expectedReward;
                double scale;
                if (futureMoments != null) {
                    expectedReward = (knownRewardSum + rewardRates[state] + futureMoments.mean[stepsLeft][state]) / (n + 1.0);
                    double varianceFloor = 0.15 * rewardSpan / Math.sqrt(n + 1.0);
                    varianceFloor *= varianceFloor;
                    double expectedVariance = futureMoments.variance[stepsLeft][state] / ((n + 1.0) * (n + 1.0));
                    scale = Math.sqrt(Math.max(expectedVariance, varianceFloor));
                } else {
                    double feasibleWidth = Math.max(feasibleMax[state] - feasibleMin[state], REWARD_TOL);
                    expectedReward = 0.5 * (feasibleMin[state] + feasibleMax[state]);
                    scale = Math.max(feasibleWidth, 0.05 * rewardSpan);
                }
                double z = (rho - expectedReward) / scale;
                double densityProxy = Math.exp(-0.5 * z * z) / scale;
                if (rewardDensityTemperature != 1.0) {
                    densityProxy = Math.pow(densityProxy, rewardDensityTemperature);
                }
                if (Double.isFinite(densityProxy) && densityProxy > 0.0) {
                    proposal[state] = baseProbs[state] * densityProxy;
                    proposalMass += proposal[state];
                }
            }
            if (proposalMass <= 0.0) {
                Arrays.fill(proposal, 0, baseProbs.length, 0.0);
                proposalMass = 0.0;
                for (int state = 0; state < baseProbs.length; state++) {
                    if (mask[state]) {
                        proposal[state] = baseProbs[state];
                        proposalMass += proposal[state];
                    }
                }
            }
        } else {
            for (int state = 0; state < baseProbs.length; state++) {
                if (mask[state]) {
                    proposal[state] = baseProbs[state];
                    proposalMass += proposal[state];
                }
            }
        }

        if (proposalMass <= 0.0) {
            throw new IllegalArgumentException("Reward-aware bridge proposal has zero feasible mass.");
        }
        for (int state = 0; state < proposal.length; state++) {
            proposal[state] /= proposalMass;
        }
        return proposal;
    }

    private static DwellResult sampleRewardConstrainedDwell(double[] c, double rho,
                                                            Options options, Random random,
                                                            boolean returnDiagnostics) {
        if (options.dwellSampler == DwellSampler.AUTO || options.dwellSampler == DwellSampler.EXACT_TWO_REWARD) {
            TwoRewardPartition partition = twoRewardPartition(c, 1.0e-12);
            if (partition != null) {
                return exactTwoRewardDwellSample(c, partition, rho, random, options.nDwellSamples);
            }
            if (options.dwellSampler == DwellSampler.EXACT_TWO_REWARD) {
                throw new IllegalArgumentException("Selected skeleton does not have exactly two reward values.");
            }
        }

        if (options.dwellSampler == DwellSampler.AUTO &&
                (Math.abs(rho - min(c)) <= REWARD_TOL || Math.abs(rho - max(c)) <= REWARD_TOL)) {
            return rewardBoundaryDwellSample(c, rho, random, options.nDwellSamples);
        }

        if (options.dwellSampler == DwellSampler.EXACT_TRIANGULATION) {
            return exactTriangulationDwellSample(c, rho, random, options.nDwellSamples);
        }

        if (options.dwellSampler == DwellSampler.AUTO && c.length <= options.exactTriangulationMaxComponents) {
            try {
                return exactTriangulationDwellSample(c, rho, random, options.nDwellSamples);
            } catch (IllegalArgumentException ignored) {
                // Fall through to hit-and-run for higher-dimensional non-simplex slices.
            }
        }

        return hitAndRunSample(c, rho, options.nHitRun, options.nBurnin, random,
                options.nDwellSamples, options.hitRunThin, options.adaptiveHitRunBurnin,
                options.maxAdaptiveHitRunBurnin);
    }

    private static DwellResult rewardBoundaryDwellSample(double[] c, double rho,
                                                         Random random, int nSamples) {
        double cMin = min(c);
        double cMax = max(c);
        boolean useMin;
        String sampler;
        if (Math.abs(rho - cMin) <= REWARD_TOL) {
            useMin = true;
            sampler = "exact_min_reward_face";
        } else if (Math.abs(rho - cMax) <= REWARD_TOL) {
            useMin = false;
            sampler = "exact_max_reward_face";
        } else {
            throw new IllegalArgumentException("rho is not on a reward boundary.");
        }

        int active = 0;
        boolean[] mask = new boolean[c.length];
        for (int i = 0; i < c.length; i++) {
            mask[i] = Math.abs(c[i] - (useMin ? cMin : cMax)) <= REWARD_TOL;
            if (mask[i]) {
                active++;
            }
        }

        double[][] retained = new double[nSamples][c.length];
        for (int sample = 0; sample < nSamples; sample++) {
            double[] dirichlet = sampleDirichlet(active, random);
            int k = 0;
            for (int i = 0; i < c.length; i++) {
                if (mask[i]) {
                    retained[sample][i] = dirichlet[k++];
                }
            }
        }

        DwellDiagnostics diagnostics = dwellDiagnostics(c, rho, retained,
                0, active < c.length ? 1 : 0, 0, 0, 1, sampler);
        return new DwellResult(retained[retained.length - 1], retained, diagnostics);
    }

    private static DwellResult exactTwoRewardDwellSample(double[] c, TwoRewardPartition partition,
                                                         double rho, Random random, int nSamples) {
        if (nSamples <= 0) {
            throw new IllegalArgumentException("nSamples must be positive.");
        }
        if (partition == null) {
            throw new IllegalArgumentException("exactTwoRewardDwellSample requires exactly two reward values.");
        }
        if (rho < partition.low - REWARD_TOL || rho > partition.high + REWARD_TOL) {
            throw new IllegalArgumentException("Polytope is empty for the requested reward.");
        }

        double highTotal = (rho - partition.low) / (partition.high - partition.low);
        highTotal = clip(highTotal, 0.0, 1.0);
        double lowTotal = 1.0 - highTotal;

        int nLow = countTrue(partition.lowMask);
        int nHigh = countTrue(partition.highMask);
        double[][] retained = new double[nSamples][c.length];
        for (int sample = 0; sample < nSamples; sample++) {
            if (lowTotal > 0.0) {
                double[] low = sampleDirichlet(nLow, random);
                int k = 0;
                for (int i = 0; i < c.length; i++) {
                    if (partition.lowMask[i]) {
                        retained[sample][i] = lowTotal * low[k++];
                    }
                }
            }
            if (highTotal > 0.0) {
                double[] high = sampleDirichlet(nHigh, random);
                int k = 0;
                for (int i = 0; i < c.length; i++) {
                    if (partition.highMask[i]) {
                        retained[sample][i] = highTotal * high[k++];
                    }
                }
            }
        }

        int boundaryHits = highTotal <= REWARD_TOL || lowTotal <= REWARD_TOL ? 1 : 0;
        DwellDiagnostics diagnostics = dwellDiagnostics(c, rho, retained,
                0, boundaryHits, 0, 0, 1, "exact_two_reward");
        return new DwellResult(retained[retained.length - 1], retained, diagnostics);
    }

    private static DwellResult exactTriangulationDwellSample(double[] c, double rho,
                                                             Random random, int nSamples) {
        if (nSamples <= 0) {
            throw new IllegalArgumentException("nSamples must be positive.");
        }
        List<double[]> vertices = polytopeVertices(c, rho, REWARD_TOL);
        double[][] retained = new double[nSamples][c.length];

        if (vertices.size() == 1) {
            for (int i = 0; i < nSamples; i++) {
                retained[i] = vertices.get(0).clone();
            }
        } else {
            double[][] diffs = new double[vertices.size() - 1][c.length];
            double[] v0 = vertices.get(0);
            for (int i = 1; i < vertices.size(); i++) {
                for (int j = 0; j < c.length; j++) {
                    diffs[i - 1][j] = vertices.get(i)[j] - v0[j];
                }
            }
            int rank = rank(diffs, 1.0e-10);
            if (rank == 0) {
                for (int i = 0; i < nSamples; i++) {
                    retained[i] = v0.clone();
                }
            } else if (rank == 1) {
                int[] endpoints = farthestPair(vertices);
                double[] left = vertices.get(endpoints[0]);
                double[] right = vertices.get(endpoints[1]);
                for (int sample = 0; sample < nSamples; sample++) {
                    double t = random.nextDouble();
                    for (int j = 0; j < c.length; j++) {
                        retained[sample][j] = left[j] + t * (right[j] - left[j]);
                    }
                }
            } else if (vertices.size() == rank + 1) {
                for (int sample = 0; sample < nSamples; sample++) {
                    double[] barycentric = sampleDirichlet(vertices.size(), random);
                    for (int vertex = 0; vertex < vertices.size(); vertex++) {
                        double weight = barycentric[vertex];
                        double[] v = vertices.get(vertex);
                        for (int j = 0; j < c.length; j++) {
                            retained[sample][j] += weight * v[j];
                        }
                    }
                }
            } else {
                throw new IllegalArgumentException(
                        "Exact triangulation requires Delaunay triangulation for this slice; use AUTO or HIT_RUN.");
            }
        }

        cleanupSmallEntries(retained, 1.0e-14);
        DwellDiagnostics diagnostics = dwellDiagnostics(c, rho, retained,
                0, anyEntryAtMost(retained, 1.0e-12) ? 1 : 0, 0, 0, 1,
                "exact_triangulation");
        diagnostics.nVertices = vertices.size();
        return new DwellResult(retained[retained.length - 1], retained, diagnostics);
    }

    private static DwellResult hitAndRunSample(double[] c, double rho, int nIter, int nBurnin,
                                               Random random, int nSamples, int thin,
                                               boolean adaptiveBurnin,
                                               int maxAdaptiveBurnin) {
        if (nSamples <= 0) {
            throw new IllegalArgumentException("nSamples must be positive.");
        }
        if (thin <= 0) {
            throw new IllegalArgumentException("thin must be positive.");
        }
        if (nIter < 0 || nBurnin < 0) {
            throw new IllegalArgumentException("nIter and nBurnin must be non-negative.");
        }

        double cMin = min(c);
        double cMax = max(c);
        if (rho < cMin - 1.0e-12 || rho > cMax + 1.0e-12) {
            throw new IllegalArgumentException("Polytope is empty for the requested reward.");
        }

        if (cMax - cMin <= 1.0e-12) {
            if (Math.abs(rho - c[0]) > REWARD_TOL) {
                throw new IllegalArgumentException("Polytope is empty: all rewards are equal.");
            }
            double[][] retained = new double[nSamples][c.length];
            for (int i = 0; i < nSamples; i++) {
                retained[i] = sampleDirichlet(c.length, random);
            }
            DwellDiagnostics diagnostics = dwellDiagnostics(c, rho, retained,
                    0, 0, 0, 0, thin, "hit_and_run");
            return new DwellResult(retained[retained.length - 1], retained, diagnostics);
        }

        List<double[]> vertices = polytopeVertices(c, rho, REWARD_TOL);
        double[] u = meanVector(vertices);
        cleanupSmallEntries(u, 1.0e-14);
        double initialMinComponent = min(u);

        List<double[]> nullBasis = nullSpaceBasisOfRewardSlice(c);
        int nullDim = nullBasis.size();
        if (nullDim == 0) {
            double[][] retained = new double[nSamples][c.length];
            for (int i = 0; i < nSamples; i++) {
                retained[i] = u.clone();
            }
            DwellDiagnostics diagnostics = dwellDiagnostics(c, rho, retained,
                    0, anyAtMost(u, 1.0e-12) ? 1 : 0, 0, 0, thin, "hit_and_run");
            diagnostics.nVertices = vertices.size();
            diagnostics.initialMinComponent = initialMinComponent;
            return new DwellResult(u, retained, diagnostics);
        }

        int effectiveBurnin = nBurnin;
        if (adaptiveBurnin) {
            int dimensionBurnin = 10 * Math.max(1, nullDim) * c.length;
            effectiveBurnin = Math.max(effectiveBurnin, dimensionBurnin);
            if (maxAdaptiveBurnin >= 0) {
                effectiveBurnin = Math.min(effectiveBurnin, maxAdaptiveBurnin);
            }
        }

        HitRunState state = new HitRunState(u);
        int totalMoves = 0;
        for (int i = 0; i < effectiveBurnin; i++) {
            takeHitAndRunStep(state, nullBasis, random);
            totalMoves++;
        }

        double[][] retained = new double[nSamples][c.length];
        for (int sample = 0; sample < nSamples; sample++) {
            int stepsToNext = sample == 0 ? nIter : thin;
            for (int step = 0; step < stepsToNext; step++) {
                takeHitAndRunStep(state, nullBasis, random);
                totalMoves++;
            }
            retained[sample] = state.u.clone();
        }

        double[] finalSample = retained[retained.length - 1];
        if (Math.abs(sum(finalSample) - 1.0) > 1.0e-8 || Math.abs(dot(c, finalSample) - rho) > 1.0e-8) {
            throw new IllegalStateException("Hit-and-run sample drifted off the affine constraints.");
        }

        DwellDiagnostics diagnostics = dwellDiagnostics(c, rho, retained,
                state.skippedMoves, state.boundaryHits, totalMoves,
                effectiveBurnin, thin, "hit_and_run");
        diagnostics.nVertices = vertices.size();
        diagnostics.initialMinComponent = initialMinComponent;
        return new DwellResult(finalSample, retained, diagnostics);
    }

    private static void takeHitAndRunStep(HitRunState state, List<double[]> nullBasis, Random random) {
        double[] direction = randomNullDirection(nullBasis, random);
        double tLo = Double.NEGATIVE_INFINITY;
        double tHi = Double.POSITIVE_INFINITY;
        for (int i = 0; i < state.u.length; i++) {
            double d = direction[i];
            if (d > 1.0e-15) {
                tLo = Math.max(tLo, -state.u[i] / d);
            } else if (d < -1.0e-15) {
                tHi = Math.min(tHi, -state.u[i] / d);
            }
        }
        if (tHi <= tLo + 1.0e-14) {
            state.skippedMoves++;
            return;
        }
        double t = tLo + random.nextDouble() * (tHi - tLo);
        if (t - tLo <= 1.0e-12 || tHi - t <= 1.0e-12) {
            state.boundaryHits++;
        }
        for (int i = 0; i < state.u.length; i++) {
            state.u[i] += t * direction[i];
            if (Math.abs(state.u[i]) < 1.0e-14) {
                state.u[i] = 0.0;
            }
            if (state.u[i] < -1.0e-10) {
                throw new IllegalStateException("Hit-and-run left the feasible polytope.");
            }
            if (state.u[i] <= 1.0e-12) {
                state.boundaryHits++;
            }
        }
    }

    private static DwellDiagnostics dwellDiagnostics(double[] c, double rho, double[][] retained,
                                                     int skippedMoves, int boundaryHits, int totalMoves,
                                                     int burnin, int thin, String sampler) {
        DwellDiagnostics diagnostics = new DwellDiagnostics();
        diagnostics.sampler = sampler;
        diagnostics.retainedSamples = copyMatrix(retained);
        diagnostics.nRetained = retained.length;
        diagnostics.skippedMoves = skippedMoves;
        diagnostics.boundaryHits = boundaryHits;
        diagnostics.totalMoves = totalMoves;
        diagnostics.acceptedMoves = totalMoves - skippedMoves;
        diagnostics.effectiveBurnin = burnin;
        diagnostics.thin = thin;

        double maxConstraintError = 0.0;
        for (int i = 0; i < retained.length; i++) {
            double sumError = sum(retained[i]) - 1.0;
            double rewardError = dot(c, retained[i]) - rho;
            maxConstraintError = Math.max(maxConstraintError, Math.max(Math.abs(sumError), Math.abs(rewardError)));
        }
        double[] finalSample = retained[retained.length - 1];
        diagnostics.sumError = sum(finalSample) - 1.0;
        diagnostics.rewardError = dot(c, finalSample) - rho;
        diagnostics.constraintError = Math.max(Math.abs(diagnostics.sumError), Math.abs(diagnostics.rewardError));
        diagnostics.maxRetainedConstraintError = maxConstraintError;
        diagnostics.selectedSample = finalSample.clone();
        return diagnostics;
    }

    private static double bsplineRewardLogDensity(double[] rewardRates, double rho, double tol) {
        double[] knots = rewardRates.clone();
        Arrays.sort(knots);
        return bsplineRewardLogDensityFromSortedKnots(knots, knots.length, rho, tol, new RewardDensityScratch());
    }

    private static double bsplineRewardLogDensityFromSortedKnots(double[] knots, int knotCount,
                                                                 double rho, double tol,
                                                                 RewardDensityScratch scratch) {
        int n = knotCount - 1;
        if (n <= 0) {
            return LOG_ZERO;
        }
        double supportWidth = knots[knotCount - 1] - knots[0];
        if (supportWidth <= tol) {
            return LOG_ZERO;
        }
        if (rho <= knots[0] + tol || rho >= knots[knotCount - 1] - tol) {
            return LOG_ZERO;
        }

        double logBasis = bsplineBasisLogValue(knots, knotCount, rho, tol, scratch);
        if (!Double.isFinite(logBasis)) {
            return LOG_ZERO;
        }
        return Math.log(n) - Math.log(supportWidth) + logBasis;
    }

    private static double bsplineBasisLogValue(double[] knots, int knotCount, double x, double tol,
                                               RewardDensityScratch scratch) {
        int degree = knotCount - 2;
        if (degree < 0) {
            return LOG_ZERO;
        }

        int nIntervals = knotCount - 1;
        scratch.ensureBasisCount(nIntervals);
        double[] logs = scratch.basisLogsA;
        double[] next = scratch.basisLogsB;
        Arrays.fill(logs, 0, nIntervals, LOG_ZERO);
        for (int i = 0; i < nIntervals; i++) {
            double left = knots[i];
            double right = knots[i + 1];
            if (right - left <= tol) {
                continue;
            }
            if ((left <= x && x < right) || (i == nIntervals - 1 && left < x && x <= right)) {
                logs[i] = 0.0;
            }
        }

        int currentLength = nIntervals;
        for (int k = 1; k <= degree; k++) {
            int nRemaining = currentLength - 1;
            for (int i = 0; i < nRemaining; i++) {
                double denomLeft = knots[i + k] - knots[i];
                double denomRight = knots[i + k + 1] - knots[i + 1];
                double logLeft = LOG_ZERO;
                double logRight = LOG_ZERO;
                if (denomLeft > tol) {
                    double coeffLeft = Math.max((x - knots[i]) / denomLeft, 0.0);
                    if (coeffLeft > 0.0 && Double.isFinite(logs[i])) {
                        logLeft = Math.log(coeffLeft) + logs[i];
                    }
                }
                if (denomRight > tol) {
                    double coeffRight = Math.max((knots[i + k + 1] - x) / denomRight, 0.0);
                    if (coeffRight > 0.0 && Double.isFinite(logs[i + 1])) {
                        logRight = Math.log(coeffRight) + logs[i + 1];
                    }
                }
                next[i] = logAdd(logLeft, logRight);
            }
            double[] tmp = logs;
            logs = next;
            next = tmp;
            currentLength = nRemaining;
        }
        return currentLength == 0 ? LOG_ZERO : logs[0];
    }

    private static Uniformization uniformize(double[][] qMatrix) {
        double omega = maxAbsDiagonal(qMatrix);
        double[][] p = new double[qMatrix.length][qMatrix.length];
        if (omega == 0.0) {
            for (int i = 0; i < p.length; i++) {
                p[i][i] = 1.0;
            }
            return new Uniformization(omega, p);
        }
        for (int i = 0; i < qMatrix.length; i++) {
            for (int j = 0; j < qMatrix.length; j++) {
                p[i][j] = (i == j ? 1.0 : 0.0) + qMatrix[i][j] / omega;
            }
        }
        clampTransitionKernel(p);
        return new Uniformization(omega, p);
    }

    private static void clampTransitionKernel(double[][] p) {
        for (int i = 0; i < p.length; i++) {
            for (int j = 0; j < p[i].length; j++) {
                if (Math.abs(p[i][j]) < 1.0e-15) {
                    p[i][j] = 0.0;
                }
                if (p[i][j] < -1.0e-12) {
                    throw new IllegalArgumentException("Uniformized kernel has a negative entry.");
                }
                if (p[i][j] < 0.0) {
                    p[i][j] = 0.0;
                }
            }
        }
    }

    private static EndpointTruncation endpointAwareTruncation(double[][] p, int x0, int xT,
                                                              double omegaT, double eps, int maxN) {
        EndpointTruncation all = endpointAwareTruncationAllStarts(p, xT, omegaT, eps, maxN);
        if (all.nMaxByStart[x0] < 0) {
            throw new IllegalArgumentException("No valid endpoint-conditioned path for start state " + x0 + ".");
        }
        int nMax = all.nMaxByStart[x0];
        double[][] messages = sliceRows(all.endpointMessages, nMax + 1);
        TruncationInfo[] info = new TruncationInfo[p.length];
        info[x0] = all.infoByStart[x0];
        int[] nMaxByStart = new int[p.length];
        Arrays.fill(nMaxByStart, -1);
        nMaxByStart[x0] = nMax;
        return new EndpointTruncation(nMaxByStart, messages, info);
    }

    private static EndpointTruncation endpointAwareTruncationAllStarts(double[][] p, int xT,
                                                                       double omegaT, double eps, int maxN) {
        if (maxN < 0) {
            double scale = Math.sqrt(Math.max(omegaT, 1.0));
            maxN = Math.max(1000, (int) (omegaT + 50.0 * scale + 500.0));
        }

        int states = p.length;
        boolean[] reachable = statesThatCanReachEndpoint(p, xT);
        int[] nByStart = new int[states];
        Arrays.fill(nByStart, -1);
        TruncationInfo[] infoByStart = new TruncationInfo[states];
        double[] endpointMass = new double[states];
        List<double[]> messageList = new ArrayList<>();

        double[] h = new double[states];
        double[] hNext = new double[states];
        h[xT] = 1.0;
        int reachableCount = countTrue(reachable);
        if (reachableCount == 0) {
            throw new IllegalArgumentException("No state can reach the requested endpoint.");
        }

        for (int n = 0; n <= maxN; n++) {
            messageList.add(h.clone());
            double poissonPmf = Math.exp(logPoissonPmf(n, omegaT));
            double tailBound = poissonTailAfter(n, omegaT);

            for (int state = 0; state < states; state++) {
                endpointMass[state] += poissonPmf * h[state];
                if (reachable[state] && nByStart[state] < 0 &&
                        endpointMass[state] > 0.0 &&
                        tailBound <= eps * endpointMass[state]) {
                    nByStart[state] = n;
                    infoByStart[state] = new TruncationInfo(n, endpointMass[state], tailBound);
                }
            }

            boolean done = true;
            for (int state = 0; state < states; state++) {
                if (reachable[state] && nByStart[state] < 0) {
                    done = false;
                    break;
                }
            }
            if (done) {
                int globalNMax = 0;
                for (int state = 0; state < states; state++) {
                    if (reachable[state]) {
                        globalNMax = Math.max(globalNMax, nByStart[state]);
                    }
                }
                double[][] messages = messageList.subList(0, globalNMax + 1).toArray(new double[0][]);
                return new EndpointTruncation(nByStart, messages, infoByStart);
            }

            matrixVectorProduct(p, h, hNext);
            double[] tmp = h; h = hNext; hNext = tmp;
        }
        throw new IllegalArgumentException("Could not satisfy endpoint-aware truncation before maxN=" + maxN + ".");
    }

    private static boolean[] statesThatCanReachEndpoint(double[][] p, int xT) {
        int states = p.length;
        boolean[] reachable = new boolean[states];
        int[] stack = new int[states];
        int stackSize = 0;
        reachable[xT] = true;
        stack[stackSize++] = xT;

        while (stackSize > 0) {
            int target = stack[--stackSize];
            for (int predecessor = 0; predecessor < states; predecessor++) {
                if (!reachable[predecessor] && p[predecessor][target] > 0.0) {
                    reachable[predecessor] = true;
                    stack[stackSize++] = predecessor;
                }
            }
        }
        return reachable;
    }

    private static FutureBounds precomputeFutureRewardBounds(double[][] p, double[][] endpointMessages,
                                                             double[] rewardRates) {
        int nMax = endpointMessages.length - 1;
        int states = p.length;
        double[][] futureMin = new double[nMax + 1][states];
        double[][] futureMax = new double[nMax + 1][states];
        for (int state = 0; state < states; state++) {
            futureMin[0][state] = Double.POSITIVE_INFINITY;
            futureMax[0][state] = Double.NEGATIVE_INFINITY;
        }

        for (int m = 1; m <= nMax; m++) {
            for (int current = 0; current < states; current++) {
                double bestMin = Double.POSITIVE_INFINITY;
                double bestMax = Double.NEGATIVE_INFINITY;
                for (int next = 0; next < states; next++) {
                    if (p[current][next] > 0.0 && endpointMessages[m - 1][next] > 0.0) {
                        bestMin = Math.min(bestMin, Math.min(rewardRates[next], futureMin[m - 1][next]));
                        bestMax = Math.max(bestMax, Math.max(rewardRates[next], futureMax[m - 1][next]));
                    }
                }
                futureMin[m][current] = bestMin;
                futureMax[m][current] = bestMax;
            }
        }
        return new FutureBounds(futureMin, futureMax);
    }

    private static FutureMoments precomputeFutureRewardMoments(double[][] p, double[][] endpointMessages,
                                                               double[] rewardRates) {
        int nMax = endpointMessages.length - 1;
        int states = p.length;
        double[][] futureMean = new double[nMax + 1][states];
        double[][] futureSecond = new double[nMax + 1][states];
        double[][] futureVariance = new double[nMax + 1][states];

        for (int m = 1; m <= nMax; m++) {
            for (int current = 0; current < states; current++) {
                double denom = endpointMessages[m][current];
                if (denom <= 0.0) {
                    continue;
                }
                double mean = 0.0;
                double second = 0.0;
                for (int next = 0; next < states; next++) {
                    double bridgeProb = p[current][next] * endpointMessages[m - 1][next] / denom;
                    if (bridgeProb <= 0.0) {
                        continue;
                    }
                    double nextMean = rewardRates[next] + futureMean[m - 1][next];
                    double nextSecond = rewardRates[next] * rewardRates[next] +
                            2.0 * rewardRates[next] * futureMean[m - 1][next] +
                            futureSecond[m - 1][next];
                    mean += bridgeProb * nextMean;
                    second += bridgeProb * nextSecond;
                }
                futureMean[m][current] = mean;
                futureSecond[m][current] = second;
                futureVariance[m][current] = Math.max(second - mean * mean, 0.0);
            }
        }
        return new FutureMoments(futureMean, futureVariance);
    }

    private static double[] endpointConditionedCountDistribution(int x0, int nMax,
                                                                 double[][] endpointMessages,
                                                                 double omegaT) {
        double[] logQ = new double[nMax + 1];
        Arrays.fill(logQ, LOG_ZERO);
        double maxLog = LOG_ZERO;
        for (int n = 0; n <= nMax; n++) {
            double bridgeProb = endpointMessages[n][x0];
            if (bridgeProb > 0.0) {
                logQ[n] = logPoissonPmf(n, omegaT) + Math.log(bridgeProb);
                maxLog = Math.max(maxLog, logQ[n]);
            }
        }
        if (!Double.isFinite(maxLog)) {
            throw new IllegalArgumentException("No valid paths from x0 to xT under the truncated bridge proposal.");
        }
        double[] q = new double[nMax + 1];
        double sum = 0.0;
        for (int n = 0; n <= nMax; n++) {
            if (Double.isFinite(logQ[n])) {
                q[n] = Math.exp(logQ[n] - maxLog);
                sum += q[n];
            }
        }
        if (sum <= 0.0) {
            throw new IllegalArgumentException("No valid paths from x0 to xT in the given branch length.");
        }
        for (int n = 0; n <= nMax; n++) {
            q[n] /= sum;
        }
        return q;
    }

    private static WeightNormalization normalizeLogWeights(ParticleStore particles) {
        double maxLog = LOG_ZERO;
        int finite = 0;
        for (int i = 0; i < particles.size(); i++) {
            double logWeight = particles.getLogWeight(i);
            if (Double.isFinite(logWeight)) {
                finite++;
                maxLog = Math.max(maxLog, logWeight);
            }
        }
        if (finite == 0) {
            throw new IllegalArgumentException("All SIR weights are zero.");
        }

        double[] weights = new double[particles.size()];
        double sum = 0.0;
        for (int i = 0; i < particles.size(); i++) {
            double logWeight = particles.getLogWeight(i);
            if (Double.isFinite(logWeight)) {
                weights[i] = Math.exp(logWeight - maxLog);
                sum += weights[i];
            }
        }
        if (sum <= 0.0) {
            throw new IllegalArgumentException("All SIR weights underflowed to zero.");
        }
        double sumSquares = 0.0;
        for (int i = 0; i < weights.length; i++) {
            weights[i] /= sum;
            sumSquares += weights[i] * weights[i];
        }
        double essRaw = 1.0 / sumSquares;
        return new WeightNormalization(weights, essRaw / weights.length, essRaw, finite, maxLog);
    }

    private static int[] resampleIndices(double[] weights, int nSamples, Random random, ResamplingMethod method) {
        int[] indices = new int[nSamples];
        if (method == ResamplingMethod.MULTINOMIAL) {
            for (int i = 0; i < nSamples; i++) {
                indices[i] = randomChoice(weights, random);
            }
            return indices;
        }

        double[] cumulative = cumulative(weights);
        cumulative[cumulative.length - 1] = 1.0;
        double systematicOffset = method == ResamplingMethod.SYSTEMATIC ? random.nextDouble() : 0.0;
        for (int i = 0; i < nSamples; i++) {
            double position;
            if (method == ResamplingMethod.SYSTEMATIC) {
                position = (systematicOffset + i) / nSamples;
            } else if (method == ResamplingMethod.STRATIFIED) {
                position = (random.nextDouble() + i) / nSamples;
            } else {
                throw new IllegalArgumentException("Unknown resampling method.");
            }
            indices[i] = searchSortedRight(cumulative, position);
        }
        return indices;
    }

    private static SelectedParticle selectResampledParticle(double[] weights, Random random, ResamplingMethod method) {
        int[] indices = resampleIndices(weights, weights.length, random, method);
        int chosenSlot = random.nextInt(indices.length);
        int unique = countUnique(indices);
        return new SelectedParticle(indices[chosenSlot], new ResamplingInfo(method, chosenSlot, unique));
    }

    private static boolean[] rewardAtomStateMask(double[] rewardRates, double rho, double tol) {
        boolean[] mask = new boolean[rewardRates.length];
        for (int i = 0; i < rewardRates.length; i++) {
            mask[i] = Math.abs(rewardRates[i] - rho) <= tol;
        }
        return mask;
    }

    private static boolean rewardAtomCompatible(double[][] qMatrix, double[] rewardRates,
                                                int x0, int xT, double rho) {
        boolean[] mask = rewardAtomStateMask(rewardRates, rho, REWARD_TOL);
        if (!mask[x0] || !mask[xT]) {
            return false;
        }
        if (x0 == xT) {
            return true;
        }

        boolean[] seen = new boolean[qMatrix.length];
        int[] stack = new int[qMatrix.length];
        int stackSize = 0;
        seen[x0] = true;
        stack[stackSize++] = x0;
        while (stackSize > 0) {
            int state = stack[--stackSize];
            for (int next = 0; next < qMatrix.length; next++) {
                if (next != state && mask[next] && !seen[next] && qMatrix[state][next] > REWARD_TOL) {
                    if (next == xT) {
                        return true;
                    }
                    seen[next] = true;
                    stack[stackSize++] = next;
                }
            }
        }
        return false;
    }

    private boolean isRewardBoundary(double rho) {
        return Math.abs(rho - minReward) <= REWARD_TOL ||
                Math.abs(rho - maxReward) <= REWARD_TOL;
    }

    private static void validateConditionedPath(int[] states, double[] times, double[] rewardRates,
                                                int x0, int xT, double branchLength, double rho) {
        if (states.length == 0 || times.length == 0 || states.length != times.length) {
            throw new IllegalStateException("Sampled path is empty or malformed.");
        }
        if (states[0] != x0 || states[states.length - 1] != xT) {
            throw new IllegalStateException("Sampled path endpoints are " + states[0] + "->" +
                    states[states.length - 1] + ", expected " + x0 + "->" + xT + ".");
        }
        for (double time : times) {
            if (time < -AFFINE_TOL) {
                throw new IllegalStateException("Sampled path contains a negative dwell time.");
            }
        }
        double totalTime = sum(times);
        if (Math.abs(totalTime - branchLength) > AFFINE_TOL * Math.max(1.0, Math.abs(branchLength))) {
            throw new IllegalStateException("Sampled path has total time " + totalTime + ", expected " + branchLength + ".");
        }
        double reward = 0.0;
        for (int i = 0; i < states.length; i++) {
            reward += rewardRates[states[i]] * times[i];
        }
        double realizedRho = reward / branchLength;
        if (Math.abs(realizedRho - rho) > AFFINE_TOL * Math.max(1.0, Math.abs(rho))) {
            throw new IllegalStateException("Sampled path has reward " + realizedRho + ", expected " + rho + ".");
        }
    }

    private static ReconstructedPath reconstructCtmcPath(int[] skeleton, double[] dwellTimes) {
        int[] statesTmp = new int[skeleton.length];
        double[] timesTmp = new double[dwellTimes.length];
        int size = 0;
        for (int i = 0; i < skeleton.length; i++) {
            double dwell = dwellTimes[i];
            if (dwell <= 1.0e-12) {
                continue;
            }
            int state = skeleton[i];
            if (size > 0 && statesTmp[size - 1] == state) {
                timesTmp[size - 1] += dwell;
            } else {
                statesTmp[size] = state;
                timesTmp[size] = dwell;
                size++;
            }
        }
        return new ReconstructedPath(Arrays.copyOf(statesTmp, size), Arrays.copyOf(timesTmp, size));
    }

    private static TwoRewardPartition twoRewardPartition(double[] c, double tol) {
        double[] sorted = c.clone();
        Arrays.sort(sorted);
        double[] values = new double[2];
        int count = 0;
        for (double value : sorted) {
            if (count == 0 || Math.abs(value - values[count - 1]) > tol) {
                if (count == 2) {
                    return null;
                }
                values[count++] = value;
            }
        }
        if (count != 2) {
            return null;
        }
        boolean[] lowMask = new boolean[c.length];
        boolean[] highMask = new boolean[c.length];
        for (int i = 0; i < c.length; i++) {
            lowMask[i] = Math.abs(c[i] - values[0]) <= tol;
            highMask[i] = Math.abs(c[i] - values[1]) <= tol;
        }
        return new TwoRewardPartition(values[0], values[1], lowMask, highMask);
    }

    private static List<double[]> polytopeVertices(double[] c, double rho, double tol) {
        List<double[]> vertices = new ArrayList<double[]>();
        for (int i = 0; i < c.length; i++) {
            if (Math.abs(c[i] - rho) <= tol) {
                double[] v = new double[c.length];
                v[i] = 1.0;
                addUniqueVertex(vertices, v, 1.0e-12);
            }
        }
        for (int i = 0; i < c.length; i++) {
            for (int j = i + 1; j < c.length; j++) {
                double denom = c[i] - c[j];
                if (Math.abs(denom) <= tol) {
                    continue;
                }
                double ui = (rho - c[j]) / denom;
                double uj = 1.0 - ui;
                if (ui >= -tol && uj >= -tol) {
                    double[] v = new double[c.length];
                    v[i] = clip(ui, 0.0, 1.0);
                    v[j] = 1.0 - v[i];
                    addUniqueVertex(vertices, v, 1.0e-12);
                }
            }
        }
        if (vertices.isEmpty()) {
            throw new IllegalArgumentException("No vertices found for reward slice.");
        }
        return vertices;
    }

    private static void addUniqueVertex(List<double[]> vertices, double[] candidate, double tol) {
        cleanupSmallEntries(candidate, 1.0e-14);
        for (double[] vertex : vertices) {
            if (maxAbsDifference(vertex, candidate) <= tol) {
                return;
            }
        }
        vertices.add(candidate);
    }

    private static List<double[]> nullSpaceBasisOfRewardSlice(double[] c) {
        int n = c.length;
        List<double[]> constraints = new ArrayList<double[]>();
        double[] ones = new double[n];
        Arrays.fill(ones, 1.0);
        normalizeInPlace(ones);
        constraints.add(ones);

        double[] rewardDirection = c.clone();
        orthogonalize(rewardDirection, constraints);
        if (norm(rewardDirection) > 1.0e-10) {
            normalizeInPlace(rewardDirection);
            constraints.add(rewardDirection);
        }

        List<double[]> basis = new ArrayList<double[]>();
        for (int coordinate = 0; coordinate < n; coordinate++) {
            double[] candidate = new double[n];
            candidate[coordinate] = 1.0;
            orthogonalize(candidate, constraints);
            orthogonalize(candidate, basis);
            double norm = norm(candidate);
            if (norm > 1.0e-10) {
                scaleInPlace(candidate, 1.0 / norm);
                basis.add(candidate);
            }
        }
        return basis;
    }

    private static double[] randomNullDirection(List<double[]> nullBasis, Random random) {
        double[] direction = new double[nullBasis.get(0).length];
        while (true) {
            Arrays.fill(direction, 0.0);
            for (double[] basisVector : nullBasis) {
                double z = random.nextGaussian();
                for (int i = 0; i < direction.length; i++) {
                    direction[i] += z * basisVector[i];
                }
            }
            double norm = norm(direction);
            if (norm >= 1.0e-14) {
                for (int i = 0; i < direction.length; i++) {
                    direction[i] /= norm;
                }
                return direction;
            }
        }
    }

    private static double[] sampleDirichlet(int dimension, Random random) {
        if (dimension <= 0) {
            throw new IllegalArgumentException("Dirichlet dimension must be positive.");
        }
        double[] values = new double[dimension];
        double total = 0.0;
        for (int i = 0; i < dimension; i++) {
            double u = Math.max(random.nextDouble(), Double.MIN_VALUE);
            values[i] = -Math.log(u);
            total += values[i];
        }
        if (total <= 0.0) {
            Arrays.fill(values, 1.0 / dimension);
            return values;
        }
        for (int i = 0; i < dimension; i++) {
            values[i] /= total;
        }
        return values;
    }

    private static double logPoissonPmf(int n, double mean) {
        if (mean == 0.0) {
            return n == 0 ? 0.0 : LOG_ZERO;
        }
        return n * Math.log(mean) - mean - Poisson.gammln(n + 1.0);
    }

    private static double poissonTailAfter(int n, double mean) {
        if (mean == 0.0) {
            return 0.0;
        }
        double tail = 1.0 - PoissonDistribution.cdf(n, mean);
        if (tail < 0.0 && tail > -1.0e-14) {
            return 0.0;
        }
        return Math.max(0.0, tail);
    }

    private void validateSamplingInputs(int x0, double rho, Options options) {
        checkStateIndex(x0, "x0");
        if (!Double.isFinite(rho)) {
            throw new IllegalArgumentException("rho must be finite.");
        }
        if (rho < minReward - REWARD_TOL || rho > maxReward + REWARD_TOL) {
            throw new IllegalArgumentException("rho is outside the attainable reward range.");
        }
        if (options.nParticles <= 0) {
            throw new IllegalArgumentException("nParticles must be positive.");
        }
        if (options.nHitRun < 0 || options.nBurnin < 0) {
            throw new IllegalArgumentException("nHitRun and nBurnin must be non-negative.");
        }
        if (options.nDwellSamples <= 0) {
            throw new IllegalArgumentException("nDwellSamples must be positive.");
        }
        if (options.hitRunThin <= 0) {
            throw new IllegalArgumentException("hitRunThin must be positive.");
        }
        if (options.essThreshold <= 0.0 || options.essThreshold > 1.0) {
            throw new IllegalArgumentException("essThreshold must lie in (0, 1].");
        }
        if (options.maxParticleMultiplier < 1) {
            throw new IllegalArgumentException("maxParticleMultiplier must be at least 1.");
        }
        if (options.rewardDensityTemperature <= 0.0) {
            throw new IllegalArgumentException("rewardDensityTemperature must be positive.");
        }
        if (options.component == null || options.dwellSampler == null || options.resamplingMethod == null) {
            throw new IllegalArgumentException("Options enums must be non-null.");
        }
    }

    private void checkStateIndex(int state, String name) {
        if (state < 0 || state >= stateCount) {
            throw new IllegalArgumentException(name + " must be a valid state index.");
        }
    }

    private static void validateStaticInputs(double[][] qMatrix, double[] rewardRates,
                                             int xT, double branchLength, double epsTruncation) {
        if (qMatrix == null || rewardRates == null) {
            throw new IllegalArgumentException("Q and rewardRates must be non-null.");
        }
        if (qMatrix.length == 0) {
            throw new IllegalArgumentException("Q must be non-empty.");
        }
        int states = qMatrix.length;
        for (int i = 0; i < states; i++) {
            if (qMatrix[i] == null || qMatrix[i].length != states) {
                throw new IllegalArgumentException("Q must be a square rate matrix.");
            }
        }
        if (rewardRates.length != states) {
            throw new IllegalArgumentException("rewardRates must contain one value per state.");
        }
        if (xT < 0 || xT >= states) {
            throw new IllegalArgumentException("xT must be a valid state index.");
        }
        if (!Double.isFinite(branchLength) || branchLength <= 0.0) {
            throw new IllegalArgumentException("branchLength must be positive and finite.");
        }
        if (epsTruncation <= 0.0 || epsTruncation >= 1.0) {
            throw new IllegalArgumentException("epsTruncation must lie in (0, 1).");
        }
        for (int i = 0; i < states; i++) {
            if (!Double.isFinite(rewardRates[i])) {
                throw new IllegalArgumentException("rewardRates contains non-finite entries.");
            }
            double rowSum = 0.0;
            for (int j = 0; j < states; j++) {
                double value = qMatrix[i][j];
                if (!Double.isFinite(value)) {
                    throw new IllegalArgumentException("Q contains non-finite entries.");
                }
                if (i != j && value < -REWARD_TOL) {
                    throw new IllegalArgumentException("Q has negative off-diagonal rates.");
                }
                rowSum += value;
            }
            if (qMatrix[i][i] > REWARD_TOL) {
                throw new IllegalArgumentException("Q has positive diagonal entries.");
            }
            if (Math.abs(rowSum) > REWARD_TOL) {
                throw new IllegalArgumentException("Rows of Q must sum to zero.");
            }
        }
    }

    private static int randomChoice(double[] probabilities, Random random) {
        double u = random.nextDouble();
        double cumulative = 0.0;
        for (int i = 0; i < probabilities.length; i++) {
            cumulative += probabilities[i];
            if (u <= cumulative) {
                return i;
            }
        }
        return probabilities.length - 1;
    }

    private static double[][] copyMatrix(double[][] matrix) {
        double[][] copy = new double[matrix.length][];
        for (int i = 0; i < matrix.length; i++) {
            copy[i] = matrix[i].clone();
        }
        return copy;
    }

    private static double[][] sliceRows(double[][] matrix, int rows) {
        double[][] slice = new double[rows][];
        for (int i = 0; i < rows; i++) {
            slice[i] = matrix[i].clone();
        }
        return slice;
    }

    private static double[] matrixVectorProduct(double[][] matrix, double[] vector) {
        double[] result = new double[matrix.length];
        for (int i = 0; i < matrix.length; i++) {
            double sum = 0.0;
            for (int j = 0; j < vector.length; j++) {
                sum += matrix[i][j] * vector[j];
            }
            result[i] = sum;
        }
        return result;
    }

    private static void matrixVectorProduct(double[][] matrix, double[] vector, double[] result) {
        for (int i = 0; i < matrix.length; i++) {
            double sum = 0.0;
            double[] row = matrix[i];
            for (int j = 0; j < vector.length; j++) {
                sum += row[j] * vector[j];
            }
            result[i] = sum;
        }
    }

    private static double[][] subMatrix(double[][] matrix, int[] indices) {
        double[][] sub = new double[indices.length][indices.length];
        for (int i = 0; i < indices.length; i++) {
            for (int j = 0; j < indices.length; j++) {
                sub[i][j] = matrix[indices[i]][indices[j]];
            }
        }
        return sub;
    }

    private static int indexOf(int[] values, int target) {
        for (int i = 0; i < values.length; i++) {
            if (values[i] == target) {
                return i;
            }
        }
        return -1;
    }

    private static int[] trueIndices(boolean[] mask) {
        int count = countTrue(mask);
        int[] indices = new int[count];
        int k = 0;
        for (int i = 0; i < mask.length; i++) {
            if (mask[i]) {
                indices[k++] = i;
            }
        }
        return indices;
    }

    private static int countTrue(boolean[] mask) {
        int count = 0;
        for (boolean value : mask) {
            if (value) {
                count++;
            }
        }
        return count;
    }

    private static double[] rewardsForSkeleton(int[] skeleton, double[] rewardRates) {
        double[] rewards = new double[skeleton.length];
        for (int i = 0; i < skeleton.length; i++) {
            rewards[i] = rewardRates[skeleton[i]];
        }
        return rewards;
    }

    private static double[] scale(double[] values, double scale) {
        double[] scaled = new double[values.length];
        for (int i = 0; i < values.length; i++) {
            scaled[i] = values[i] * scale;
        }
        return scaled;
    }

    private static double min(double[] values) {
        double min = Double.POSITIVE_INFINITY;
        for (double value : values) {
            min = Math.min(min, value);
        }
        return min;
    }

    private static double max(double[] values) {
        double max = Double.NEGATIVE_INFINITY;
        for (double value : values) {
            max = Math.max(max, value);
        }
        return max;
    }

    private static double maxAbsDiagonal(double[][] matrix) {
        double max = 0.0;
        for (int i = 0; i < matrix.length; i++) {
            max = Math.max(max, Math.abs(matrix[i][i]));
        }
        return max;
    }

    private static double sum(double[] values) {
        double sum = 0.0;
        for (double value : values) {
            sum += value;
        }
        return sum;
    }

    private static double dot(double[] x, double[] y) {
        double dot = 0.0;
        for (int i = 0; i < x.length; i++) {
            dot += x[i] * y[i];
        }
        return dot;
    }

    private static double norm(double[] values) {
        return Math.sqrt(dot(values, values));
    }

    private static void normalizeInPlace(double[] values) {
        double norm = norm(values);
        if (norm <= 0.0) {
            return;
        }
        scaleInPlace(values, 1.0 / norm);
    }

    private static void scaleInPlace(double[] values, double scale) {
        for (int i = 0; i < values.length; i++) {
            values[i] *= scale;
        }
    }

    private static void orthogonalize(double[] candidate, List<double[]> basis) {
        for (double[] basisVector : basis) {
            double projection = dot(candidate, basisVector);
            for (int i = 0; i < candidate.length; i++) {
                candidate[i] -= projection * basisVector[i];
            }
        }
    }

    private static double logAdd(double x, double y) {
        if (!Double.isFinite(x)) {
            return y;
        }
        if (!Double.isFinite(y)) {
            return x;
        }
        if (x < y) {
            double tmp = x;
            x = y;
            y = tmp;
        }
        return x + Math.log1p(Math.exp(y - x));
    }

    private static double clip(double value, double lower, double upper) {
        return Math.max(lower, Math.min(upper, value));
    }

    private static double[] cumulative(double[] weights) {
        double[] cumulative = new double[weights.length];
        double sum = 0.0;
        for (int i = 0; i < weights.length; i++) {
            sum += weights[i];
            cumulative[i] = sum;
        }
        return cumulative;
    }

    private static int searchSortedRight(double[] sorted, double value) {
        int low = 0;
        int high = sorted.length;
        while (low < high) {
            int mid = (low + high) >>> 1;
            if (value < sorted[mid]) {
                high = mid;
            } else {
                low = mid + 1;
            }
        }
        return Math.min(low, sorted.length - 1);
    }

    private static int countUnique(int[] values) {
        int[] copy = values.clone();
        Arrays.sort(copy);
        int unique = 0;
        int previous = Integer.MIN_VALUE;
        for (int value : copy) {
            if (unique == 0 || value != previous) {
                unique++;
                previous = value;
            }
        }
        return unique;
    }

    private static double maxAbsDifference(double[] x, double[] y) {
        double max = 0.0;
        for (int i = 0; i < x.length; i++) {
            max = Math.max(max, Math.abs(x[i] - y[i]));
        }
        return max;
    }

    private static double[] meanVector(List<double[]> vectors) {
        double[] mean = new double[vectors.get(0).length];
        for (double[] vector : vectors) {
            for (int i = 0; i < mean.length; i++) {
                mean[i] += vector[i];
            }
        }
        for (int i = 0; i < mean.length; i++) {
            mean[i] /= vectors.size();
        }
        return mean;
    }

    private static int rank(double[][] matrix, double tol) {
        if (matrix.length == 0) {
            return 0;
        }
        double[][] a = copyMatrix(matrix);
        int rows = a.length;
        int cols = a[0].length;
        int rank = 0;
        for (int col = 0; col < cols && rank < rows; col++) {
            int pivot = rank;
            for (int row = rank + 1; row < rows; row++) {
                if (Math.abs(a[row][col]) > Math.abs(a[pivot][col])) {
                    pivot = row;
                }
            }
            if (Math.abs(a[pivot][col]) <= tol) {
                continue;
            }
            double[] tmp = a[rank];
            a[rank] = a[pivot];
            a[pivot] = tmp;
            double pivotValue = a[rank][col];
            for (int j = col; j < cols; j++) {
                a[rank][j] /= pivotValue;
            }
            for (int row = 0; row < rows; row++) {
                if (row == rank) {
                    continue;
                }
                double factor = a[row][col];
                for (int j = col; j < cols; j++) {
                    a[row][j] -= factor * a[rank][j];
                }
            }
            rank++;
        }
        return rank;
    }

    private static int[] farthestPair(List<double[]> vertices) {
        int bestI = 0;
        int bestJ = 1;
        double bestDistance = -1.0;
        for (int i = 0; i < vertices.size(); i++) {
            for (int j = i + 1; j < vertices.size(); j++) {
                double distance = squaredDistance(vertices.get(i), vertices.get(j));
                if (distance > bestDistance) {
                    bestDistance = distance;
                    bestI = i;
                    bestJ = j;
                }
            }
        }
        return new int[]{bestI, bestJ};
    }

    private static double squaredDistance(double[] x, double[] y) {
        double distance = 0.0;
        for (int i = 0; i < x.length; i++) {
            double diff = x[i] - y[i];
            distance += diff * diff;
        }
        return distance;
    }

    private static boolean anyAtMost(double[] values, double threshold) {
        for (double value : values) {
            if (value <= threshold) {
                return true;
            }
        }
        return false;
    }

    private static boolean anyEntryAtMost(double[][] values, double threshold) {
        for (double[] row : values) {
            if (anyAtMost(row, threshold)) {
                return true;
            }
        }
        return false;
    }

    private static void cleanupSmallEntries(double[] values, double threshold) {
        for (int i = 0; i < values.length; i++) {
            if (Math.abs(values[i]) < threshold) {
                values[i] = 0.0;
            }
        }
    }

    private static void cleanupSmallEntries(double[][] values, double threshold) {
        for (double[] row : values) {
            cleanupSmallEntries(row, threshold);
        }
    }

    private static final class ParticleStore {
        private int[][] skeletons;
        private double[] logWeights;
        private int size;

        private ParticleStore(int initialCapacity) {
            int capacity = Math.max(1, initialCapacity);
            this.skeletons = new int[capacity][];
            this.logWeights = new double[capacity];
        }

        private void add(int[] skeleton, double logWeight) {
            ensureCapacity(size + 1);
            skeletons[size] = skeleton;
            logWeights[size] = logWeight;
            size++;
        }

        private int size() {
            return size;
        }

        private int[] getSkeleton(int index) {
            if (index < 0 || index >= size) {
                throw new IndexOutOfBoundsException("particle index " + index);
            }
            return skeletons[index];
        }

        private double getLogWeight(int index) {
            if (index < 0 || index >= size) {
                throw new IndexOutOfBoundsException("particle index " + index);
            }
            return logWeights[index];
        }

        private void ensureCapacity(int targetCapacity) {
            if (targetCapacity <= skeletons.length) {
                return;
            }
            int newCapacity = Math.max(targetCapacity, skeletons.length * 2);
            skeletons = Arrays.copyOf(skeletons, newCapacity);
            logWeights = Arrays.copyOf(logWeights, newCapacity);
        }
    }

    private static final class SkeletonScratch {
        private double[] weights;
        private double[] baseProbs;
        private double[] proposalProbs;
        private double[] feasibleMin;
        private double[] feasibleMax;
        private boolean[] mask;

        private SkeletonScratch(int stateCount) {
            ensureStateCount(stateCount);
        }

        private void ensureStateCount(int stateCount) {
            if (weights != null && weights.length >= stateCount) {
                return;
            }
            weights = new double[stateCount];
            baseProbs = new double[stateCount];
            proposalProbs = new double[stateCount];
            feasibleMin = new double[stateCount];
            feasibleMax = new double[stateCount];
            mask = new boolean[stateCount];
        }
    }

    private static final class RewardDensityScratch {
        private double[] sortedKnots = new double[0];
        private double[] basisLogsA = new double[0];
        private double[] basisLogsB = new double[0];

        private void ensureKnotCount(int knotCount) {
            if (sortedKnots.length < knotCount) {
                sortedKnots = new double[knotCount];
            }
            ensureBasisCount(Math.max(0, knotCount - 1));
        }

        private void ensureBasisCount(int basisCount) {
            if (basisLogsA.length < basisCount) {
                basisLogsA = new double[basisCount];
                basisLogsB = new double[basisCount];
            }
        }
    }

    private static final class Uniformization {
        private final double omega;
        private final double[][] kernel;

        private Uniformization(double omega, double[][] kernel) {
            this.omega = omega;
            this.kernel = kernel;
        }
    }

    private static final class EndpointTruncation {
        private final int[] nMaxByStart;
        private final double[][] endpointMessages;
        private final TruncationInfo[] infoByStart;

        private EndpointTruncation(int[] nMaxByStart, double[][] endpointMessages,
                                   TruncationInfo[] infoByStart) {
            this.nMaxByStart = nMaxByStart;
            this.endpointMessages = endpointMessages;
            this.infoByStart = infoByStart;
        }
    }

    private static final class FutureBounds {
        private final double[][] min;
        private final double[][] max;

        private FutureBounds(double[][] min, double[][] max) {
            this.min = min;
            this.max = max;
        }
    }

    private static final class FutureMoments {
        private final double[][] mean;
        private final double[][] variance;

        private FutureMoments(double[][] mean, double[][] variance) {
            this.mean = mean;
            this.variance = variance;
        }
    }

    private static final class SkeletonResult {
        private final int[] skeleton;
        private final double logBridgeCorrection;

        private SkeletonResult(int[] skeleton, double logBridgeCorrection) {
            this.skeleton = skeleton;
            this.logBridgeCorrection = logBridgeCorrection;
        }
    }

    private static final class DwellResult {
        private final double[] finalSample;
        private final double[][] retainedSamples;
        private final DwellDiagnostics diagnostics;

        private DwellResult(double[] finalSample, double[][] retainedSamples,
                            DwellDiagnostics diagnostics) {
            this.finalSample = finalSample;
            this.retainedSamples = retainedSamples;
            this.diagnostics = diagnostics;
        }
    }

    private static final class ReconstructedPath {
        private final int[] states;
        private final double[] times;

        private ReconstructedPath(int[] states, double[] times) {
            this.states = states;
            this.times = times;
        }
    }

    private static final class WeightNormalization {
        private final double[] weights;
        private final double essFraction;
        private final double essRaw;
        private final int nFiniteParticles;
        private final double maxLogWeight;

        private WeightNormalization(double[] weights, double essFraction,
                                    double essRaw, int nFiniteParticles,
                                    double maxLogWeight) {
            this.weights = weights;
            this.essFraction = essFraction;
            this.essRaw = essRaw;
            this.nFiniteParticles = nFiniteParticles;
            this.maxLogWeight = maxLogWeight;
        }
    }

    private static final class SelectedParticle {
        private final int index;
        private final ResamplingInfo resamplingInfo;

        private SelectedParticle(int index, ResamplingInfo resamplingInfo) {
            this.index = index;
            this.resamplingInfo = resamplingInfo;
        }
    }

    private static final class TwoRewardPartition {
        private final double low;
        private final double high;
        private final boolean[] lowMask;
        private final boolean[] highMask;

        private TwoRewardPartition(double low, double high, boolean[] lowMask, boolean[] highMask) {
            this.low = low;
            this.high = high;
            this.lowMask = lowMask;
            this.highMask = highMask;
        }
    }

    private static final class HitRunState {
        private final double[] u;
        private int skippedMoves;
        private int boundaryHits;

        private HitRunState(double[] u) {
            this.u = u;
        }
    }
}

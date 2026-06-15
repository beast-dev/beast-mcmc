/*
 * RewardsAwareFiftyTipRecoveryExperiment.java
 *
 * Copyright (c) 2002-2026 the BEAST Development Team
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package test.dr.inference.operators;

import dr.app.beagle.tools.BeagleSequenceSimulator;
import dr.app.beagle.tools.Partition;
import dr.evolution.alignment.SimpleAlignment;
import dr.evolution.alignment.SitePatterns;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.branchmodel.HomogeneousBranchModel;
import dr.evomodel.branchmodel.RewardsAwareBranchModel;
import dr.evomodel.branchmodel.RewardsAwareBranchModelRewardRateGradient;
import dr.evomodel.branchmodel.RewardsAwareCtmcSimulation;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.branchratemodel.RewardRates;
import dr.evomodel.branchratemodel.RewardsAwareMixtureBranchRates;
import dr.evomodel.siteratemodel.GammaSiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.LogAdditiveCtmcRateProvider;
import dr.evomodel.substmodel.LogRateSubstitutionModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.substmodel.nucleotide.HKY;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeParameterModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DiscreteDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.RewardsAwareMixtureBranchRatesRewardRateGradient;
import dr.evomodel.treedatalikelihood.discrete.discretetreedataLikelihood.representations.RewardsAwarePartialsRepresentation;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.inference.operators.RewardsMixtureIndicatorAndAtomIndicesOperator;
import dr.math.MathUtils;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Small executable recovery experiment for the reward-aware CTMC machinery.
 *
 * The experiment conditions on the simulated continuous reward proportions
 * and asks whether the indicator/atom operator can recover atomic branches
 * and atom states from the independent reward-process tips, with and without
 * three dependent CTMC likelihoods.
 *
 * Usage:
 *   java ... test.dr.inference.operators.RewardsAwareFiftyTipRecoveryExperiment [outputDir] [sweeps] [burnIn] [sitesPerDependent]
 *
 * Set -DrewardRecovery.dependentUseSpectral=true to diagnose BEAGLE full-likelihood
 * spectral representation for the dependent CTMC likelihoods.
 * Set -DrewardRecovery.printLogTargetComponents=true to print per-likelihood components.
 * Tune reward-rate HMC with -DrewardRecovery.rewardRateHmcStepSize,
 * -DrewardRecovery.rewardRateHmcLeapfrogSteps, and -DrewardRecovery.rewardRateHmcSteps.
 *
 * @author Filippo Monti
 */
public final class RewardsAwareFiftyTipRecoveryExperiment {

    private static final int TIP_COUNT = 50;
    private static final int STATE_COUNT = 4;
    private static final int DEPENDENT_COUNT = 3;
    private static final double BRANCH_HEIGHT_UNIT = 0.18;
    private static final long SEED = 20260611L;
    private static final int DEFAULT_SWEEPS = 600;
    private static final int DEFAULT_BURN_IN = 100;
    private static final int DEFAULT_DEPENDENT_SITE_COUNT = 250;
    private static final boolean DEPENDENT_USE_SPECTRAL =
            Boolean.getBoolean("rewardRecovery.dependentUseSpectral");
    private static final boolean PRINT_LOG_TARGET_COMPONENTS =
            Boolean.getBoolean("rewardRecovery.printLogTargetComponents");
    private static final int REWARD_RATE_HMC_STEPS =
            Integer.getInteger("rewardRecovery.rewardRateHmcSteps", 1);
    private static final int REWARD_RATE_HMC_LEAPFROG_STEPS =
            Integer.getInteger("rewardRecovery.rewardRateHmcLeapfrogSteps", 5);
    private static final double REWARD_RATE_HMC_STEP_SIZE =
            getDoubleProperty("rewardRecovery.rewardRateHmcStepSize", 0.01);

    private static final char[] NUCLEOTIDE_CODES = {'A', 'C', 'G', 'T'};
    private static final double[] INDEPENDENT_ROOT_FREQUENCIES = {0.001, 0.001, 0.499, 0.499};
    private static final double[] TRUE_REWARD_RATE_VALUES = {0.15, 0.85, 0.35, 0.70};
    private static final double[] INITIAL_REWARD_RATE_VARYING_VALUES = {0.25, 0.80};

    private RewardsAwareFiftyTipRecoveryExperiment() {
        // executable utility
    }

    private static double getDoubleProperty(final String name, final double defaultValue) {
        final String raw = System.getProperty(name);
        return raw == null ? defaultValue : Double.parseDouble(raw);
    }

    public static void main(final String[] args) throws Exception {
        final File outputDir = args.length > 0
                ? new File(args[0])
                : new File("../rewardsAwareEmpiricalAnalysis");
        final int sweeps = args.length > 1 ? Integer.parseInt(args[1]) : DEFAULT_SWEEPS;
        final int burnIn = args.length > 2 ? Integer.parseInt(args[2]) : DEFAULT_BURN_IN;
        final int dependentSiteCount = args.length > 3
                ? Integer.parseInt(args[3])
                : DEFAULT_DEPENDENT_SITE_COUNT;

        if (sweeps <= 0) {
            throw new IllegalArgumentException("sweeps must be positive");
        }
        if (burnIn < 0 || burnIn >= sweeps) {
            throw new IllegalArgumentException("burnIn must be in [0, sweeps)");
        }
        if (dependentSiteCount <= 0) {
            throw new IllegalArgumentException("sitesPerDependent must be positive");
        }
        if (!outputDir.exists() && !outputDir.mkdirs()) {
            throw new IOException("Could not create output directory: " + outputDir);
        }

        MathUtils.setSeed(SEED);

        final TreeModel tree = createBalancedTree(TIP_COUNT);
        final RewardRates truthRewardRates = createFixedRewardRates("truth.rewardRates", TRUE_REWARD_RATE_VALUES);
        final LogRateSubstitutionModel independentModel = createIndependentRewardSubstitutionModel("independent.model");
        final double[] qMatrix = new double[STATE_COUNT * STATE_COUNT];
        independentModel.getInfinitesimalMatrix(qMatrix);

        final RewardsAwareCtmcSimulation simulation =
                RewardsAwareCtmcSimulation.withRootFrequencies(
                        "fiftyTip.independentSimulation",
                        STATE_COUNT,
                        truthRewardRates.getValues().getParameterValues(),
                        tree,
                        qMatrix,
                        INDEPENDENT_ROOT_FREQUENCIES
                );

        final double[] trueCts = simulation.getBranchRewardProportions();
        final int[] trueIndicators = simulation.getIndicators();
        final int[] trueAtoms = simulation.getAtomIndices();

        final RewardsAwareMixtureBranchRates trueBranchRates =
                createRewardBranchRates(tree, trueCts, toDouble(trueIndicators), toDouble(trueAtoms),
                        truthRewardRates, "truth");

        final SimpleAlignment independentAlignment =
                createIndependentTipAlignment(tree, simulation.getGeneratedTipStates());
        final SitePatterns independentPatterns = new SitePatterns(independentAlignment, null, 0, -1, 1, true);

        final SimpleAlignment[] dependentAlignments = new SimpleAlignment[DEPENDENT_COUNT];
        final SitePatterns[] dependentPatterns = new SitePatterns[DEPENDENT_COUNT];
        final SubstitutionModel[] dependentModels = new SubstitutionModel[DEPENDENT_COUNT];
        final GammaSiteRateModel[] dependentSiteModels = new GammaSiteRateModel[DEPENDENT_COUNT];
        for (int i = 0; i < DEPENDENT_COUNT; i++) {
            dependentModels[i] = createDependentSubstitutionModel(i);
            dependentSiteModels[i] = new GammaSiteRateModel("dependent.siteModel." + (i + 1));
            dependentAlignments[i] = simulateDependentAlignment(
                    tree,
                    trueBranchRates,
                    dependentModels[i],
                    dependentSiteModels[i],
                    dependentSiteCount,
                    "dependent." + (i + 1)
            );
            dependentPatterns[i] = new SitePatterns(dependentAlignments[i], null, 0, -1, 1, true);
        }

        final RunResult independentOnly = runRecovery(
                "independentOnly",
                tree,
                TRUE_REWARD_RATE_VALUES,
                trueCts,
                independentModel,
                independentPatterns,
                dependentModels,
                dependentSiteModels,
                new SitePatterns[0],
                trueIndicators,
                trueAtoms,
                sweeps,
                burnIn
        );

        MathUtils.setSeed(SEED + 1);
        final RunResult withDependents = runRecovery(
                "independentPlusThreeDependents",
                tree,
                TRUE_REWARD_RATE_VALUES,
                trueCts,
                independentModel,
                independentPatterns,
                dependentModels,
                dependentSiteModels,
                dependentPatterns,
                trueIndicators,
                trueAtoms,
                sweeps,
                burnIn
        );

        writeSummary(
                new File(outputDir, "fifty_tip_reward_recovery_summary.tsv"),
                tree,
                simulation,
                dependentSiteCount,
                sweeps,
                burnIn,
                independentOnly,
                withDependents
        );
        writeBranchSummary(
                new File(outputDir, "fifty_tip_reward_recovery_by_branch.tsv"),
                tree,
                simulation,
                trueCts,
                trueIndicators,
                trueAtoms,
                independentOnly,
                withDependents
        );
        writeNewick(new File(outputDir, "fifty_tip_balanced_tree.nwk"), tree);
        writeAlignment(new File(outputDir, "fifty_tip_independent_tips.fasta"), independentAlignment);
        for (int i = 0; i < DEPENDENT_COUNT; i++) {
            writeAlignment(new File(outputDir, "fifty_tip_dependent_" + (i + 1) + ".fasta"),
                    dependentAlignments[i]);
        }

        System.out.println("Wrote recovery outputs to " + outputDir.getAbsolutePath());
        System.out.println(independentOnly.toConsoleString());
        System.out.println(withDependents.toConsoleString());
    }

    private static RunResult runRecovery(final String label,
                                         final TreeModel tree,
                                         final double[] trueRewardRateValues,
                                         final double[] trueCts,
                                         final LogRateSubstitutionModel independentModel,
                                         final SitePatterns independentPatterns,
                                         final SubstitutionModel[] dependentModels,
                                         final GammaSiteRateModel[] dependentSiteModels,
                                         final SitePatterns[] dependentPatterns,
                                         final int[] trueIndicators,
                                         final int[] trueAtoms,
                                         final int sweeps,
                                         final int burnIn) {

        final double[] initialIndicators = toDouble(trueIndicators);
        final double[] initialAtoms = toDouble(trueAtoms);
        final RewardRates rewardRates = createInferenceRewardRates(label + ".rewardRates",
                trueRewardRateValues[0],
                trueRewardRateValues[1]);

        final Parameter ctsRewards = new Parameter.Default(label + ".ctsRewards", Arrays.copyOf(trueCts, trueCts.length));
        final Parameter indicator = new Parameter.Default(label + ".indicator", initialIndicators);
        final Parameter atomIndices = new Parameter.Default(label + ".atomIndices", initialAtoms);

        final RewardsAwareMixtureBranchRates branchRates = new RewardsAwareMixtureBranchRates(
                tree,
                ctsRewards,
                indicator,
                atomIndices,
                rewardRates,
                new ArbitraryBranchRates.BranchRateTransform.None(),
                false,
                TreeParameterModel.Type.WITHOUT_ROOT
        );
        final RewardsAwareBranchModel rewardsAwareBranchModel = new RewardsAwareBranchModel(
                tree,
                independentModel,
                rewardRates,
                indicator,
                branchRates,
                atomIndices,
                false
        );

        final GammaSiteRateModel independentSiteRateModel = new GammaSiteRateModel(label + ".independent.siteModel");
        final DiscreteDataLikelihoodDelegate independentDelegate = new DiscreteDataLikelihoodDelegate(
                tree,
                independentPatterns,
                rewardsAwareBranchModel,
                independentSiteRateModel,
                false,
                false,
                PartialsRescalingScheme.NONE,
                false,
                new PreOrderSettings(true, false, false, false),
                new RewardsAwarePartialsRepresentation(rewardsAwareBranchModel),
                DiscreteDataLikelihoodDelegate.PartialTransform.IDENTITY,
                DiscreteDataLikelihoodDelegate.PartialTransform.IDENTITY,
                null,
                null
        );
        final TreeDataLikelihood independentLikelihood =
                new TreeDataLikelihood(independentDelegate, tree, new DefaultBranchRateModel());

        final TreeDataLikelihood[] dependentLikelihoods = new TreeDataLikelihood[dependentPatterns.length];
        for (int i = 0; i < dependentLikelihoods.length; i++) {
            final BeagleDataLikelihoodDelegate dependentDelegate = new BeagleDataLikelihoodDelegate(
                    tree,
                    dependentPatterns[i],
                    new HomogeneousBranchModel(dependentModels[i]),
                    dependentSiteModels[i],
                    DEPENDENT_USE_SPECTRAL,
                    false,
                    PartialsRescalingScheme.ALWAYS,
                    false,
                    new PreOrderSettings(true, true, false, false, DEPENDENT_USE_SPECTRAL)
            );
            dependentLikelihoods[i] = new TreeDataLikelihood(dependentDelegate, tree, branchRates);
        }

        final RewardRateHmcKernel rewardRateHmc = createRewardRateHmcKernel(
                independentLikelihood,
                dependentLikelihoods,
                rewardsAwareBranchModel,
                branchRates,
                rewardRates,
                trueRewardRateValues[0],
                trueRewardRateValues[1]);

        final RewardsMixtureIndicatorAndAtomIndicesOperator operator =
                new RewardsMixtureIndicatorAndAtomIndicesOperator(
                        indicator,
                        atomIndices,
                        rewardsAwareBranchModel,
                        independentLikelihood,
                        dependentLikelihoods,
                        1.0 / trueCts.length,
                        false,
                        false,
                        2,
                        0.5,
                        1.0
                );

        final double initialLogTarget = logTarget(independentLikelihood, dependentLikelihoods);
        final int branchCount = trueCts.length;
        final int[] atomicSampleCounts = new int[branchCount];
        final int[][] atomSampleCounts = new int[branchCount][STATE_COUNT];
        final int retainedSampleCount = sweeps - burnIn;
        final double[][] indicatorSamplesByBranch = new double[branchCount][retainedSampleCount];
        final double[] logTargetSamples = new double[retainedSampleCount];
        final int rewardRateDimension = rewardRates.getVaryingValues().getDimension();
        final double[][] rewardRateSamplesByIndex =
                new double[rewardRateDimension][retainedSampleCount];
        int sampleCount = 0;

        for (int sweep = 0; sweep < sweeps; sweep++) {
            operator.doOperation();
            for (int h = 0; h < REWARD_RATE_HMC_STEPS; h++) {
                rewardRateHmc.doOperation();
            }
            if (sweep >= burnIn) {
                final int sampleIndex = sampleCount++;
                for (int branch = 0; branch < branchCount; branch++) {
                    final int z = (int) Math.round(indicator.getParameterValue(branch));
                    indicatorSamplesByBranch[branch][sampleIndex] = z;
                    if (z == 1) {
                        atomicSampleCounts[branch]++;
                        final int atom = (int) Math.round(atomIndices.getParameterValue(branch));
                        if (atom >= 0 && atom < STATE_COUNT) {
                            atomSampleCounts[branch][atom]++;
                        }
                    }
                }
                for (int r = 0; r < rewardRateDimension; r++) {
                    rewardRateSamplesByIndex[r][sampleIndex] =
                            rewardRates.getVaryingValues().getParameterValue(r);
                }
                logTargetSamples[sampleIndex] = logTarget(independentLikelihood, dependentLikelihoods);
            }
        }

        final double finalLogTarget = logTarget(independentLikelihood, dependentLikelihoods);
        final double[] indicatorEss = new double[branchCount];
        for (int branch = 0; branch < branchCount; branch++) {
            indicatorEss[branch] = effectiveSampleSize(indicatorSamplesByBranch[branch]);
        }
        return summarizeRun(
                label,
                trueIndicators,
                trueAtoms,
                atomicSampleCounts,
                atomSampleCounts,
                indicatorEss,
                effectiveSampleSize(logTargetSamples),
                Arrays.copyOfRange(trueRewardRateValues, 2, trueRewardRateValues.length),
                Arrays.copyOf(INITIAL_REWARD_RATE_VARYING_VALUES, INITIAL_REWARD_RATE_VARYING_VALUES.length),
                rewardRateSamplesByIndex,
                rewardRates.getVaryingValues().getParameterValues(),
                rewardRateHmc.getAttemptCount(),
                rewardRateHmc.getAcceptCount(),
                sampleCount,
                initialLogTarget,
                finalLogTarget
        );
    }

    private static RunResult summarizeRun(final String label,
                                          final int[] trueIndicators,
                                          final int[] trueAtoms,
                                          final int[] atomicSampleCounts,
                                          final int[][] atomSampleCounts,
                                          final double[] indicatorEss,
                                          final double logTargetEss,
                                          final double[] trueRewardRateVaryingValues,
                                          final double[] initialRewardRateVaryingValues,
                                          final double[][] rewardRateSamplesByIndex,
                                          final double[] finalRewardRateVaryingValues,
                                          final int rewardRateHmcAttemptCount,
                                          final int rewardRateHmcAcceptCount,
                                          final int sampleCount,
                                          final double initialLogTarget,
                                          final double finalLogTarget) {

        final int branchCount = trueIndicators.length;
        final double[] posteriorAtomicProb = new double[branchCount];
        final int[] atomCalls = new int[branchCount];

        int indicatorCorrect = 0;
        int trueAtomic = 0;
        int calledAtomic = 0;
        int truePositiveAtomic = 0;
        int atomCorrectWhenTrueAtomic = 0;
        int confidentAtomCorrect = 0;
        double brier = 0.0;
        double meanAtomicProbTrueAtomic = 0.0;
        double meanAtomicProbTrueContinuous = 0.0;
        int trueContinuous = 0;
        final double[] rewardRatePosteriorMeans = new double[rewardRateSamplesByIndex.length];
        final double[] rewardRateEss = new double[rewardRateSamplesByIndex.length];

        for (int i = 0; i < rewardRateSamplesByIndex.length; i++) {
            rewardRatePosteriorMeans[i] = mean(rewardRateSamplesByIndex[i]);
            rewardRateEss[i] = effectiveSampleSize(rewardRateSamplesByIndex[i]);
        }

        for (int branch = 0; branch < branchCount; branch++) {
            posteriorAtomicProb[branch] = atomicSampleCounts[branch] / (double) sampleCount;
            atomCalls[branch] = argMax(atomSampleCounts[branch]);

            final int trueZ = trueIndicators[branch];
            final int calledZ = posteriorAtomicProb[branch] >= 0.5 ? 1 : 0;
            if (trueZ == calledZ) {
                indicatorCorrect++;
            }
            if (trueZ == 1) {
                trueAtomic++;
                meanAtomicProbTrueAtomic += posteriorAtomicProb[branch];
                if (atomCalls[branch] == trueAtoms[branch]) {
                    atomCorrectWhenTrueAtomic++;
                }
                if (calledZ == 1 && atomCalls[branch] == trueAtoms[branch]) {
                    confidentAtomCorrect++;
                }
            } else {
                trueContinuous++;
                meanAtomicProbTrueContinuous += posteriorAtomicProb[branch];
            }
            if (calledZ == 1) {
                calledAtomic++;
            }
            if (trueZ == 1 && calledZ == 1) {
                truePositiveAtomic++;
            }
            final double error = posteriorAtomicProb[branch] - trueZ;
            brier += error * error;
        }

        return new RunResult(
                label,
                sampleCount,
                initialLogTarget,
                finalLogTarget,
                posteriorAtomicProb,
                atomCalls,
                indicatorCorrect / (double) branchCount,
                trueAtomic == 0 ? Double.NaN : truePositiveAtomic / (double) trueAtomic,
                calledAtomic == 0 ? Double.NaN : truePositiveAtomic / (double) calledAtomic,
                trueAtomic == 0 ? Double.NaN : atomCorrectWhenTrueAtomic / (double) trueAtomic,
                trueAtomic == 0 ? Double.NaN : confidentAtomCorrect / (double) trueAtomic,
                brier / branchCount,
                trueAtomic == 0 ? Double.NaN : meanAtomicProbTrueAtomic / trueAtomic,
                trueContinuous == 0 ? Double.NaN : meanAtomicProbTrueContinuous / trueContinuous,
                calledAtomic,
                indicatorEss,
                mean(indicatorEss),
                median(indicatorEss),
                min(indicatorEss),
                max(indicatorEss),
                countAtLeast(indicatorEss, 0.95 * sampleCount),
                logTargetEss,
                trueRewardRateVaryingValues,
                initialRewardRateVaryingValues,
                rewardRatePosteriorMeans,
                finalRewardRateVaryingValues,
                rewardRateEss,
                rewardRateHmcAttemptCount,
                rewardRateHmcAcceptCount
        );
    }

    private static double logTarget(final TreeDataLikelihood independentLikelihood,
                                    final TreeDataLikelihood[] dependentLikelihoods) {
        independentLikelihood.makeDirty();
        final double independentLogLikelihood = independentLikelihood.getLogLikelihood();
        double logTarget = independentLogLikelihood;
        if (PRINT_LOG_TARGET_COMPONENTS) {
            System.err.println("logTarget.independent=" + independentLogLikelihood);
        }
        for (int i = 0; i < dependentLikelihoods.length; i++) {
            final TreeDataLikelihood dependentLikelihood = dependentLikelihoods[i];
            dependentLikelihood.makeDirty();
            final double dependentLogLikelihood = dependentLikelihood.getLogLikelihood();
            if (PRINT_LOG_TARGET_COMPONENTS) {
                System.err.println("logTarget.dependent." + (i + 1) + "=" + dependentLogLikelihood);
                if (!Double.isFinite(dependentLogLikelihood)) {
                    printRootPartialDiagnostics("dependent." + (i + 1), dependentLikelihood);
                }
            }
            logTarget += dependentLogLikelihood;
        }
        if (Double.isNaN(logTarget) || logTarget == Double.POSITIVE_INFINITY) {
            throw new IllegalStateException("Invalid log target: " + logTarget);
        }
        return logTarget;
    }

    private static void printRootPartialDiagnostics(final String label,
                                                    final TreeDataLikelihood likelihood) {
        if (!(likelihood.getDataLikelihoodDelegate() instanceof BeagleDataLikelihoodDelegate)) {
            return;
        }

        final BeagleDataLikelihoodDelegate delegate =
                (BeagleDataLikelihoodDelegate) likelihood.getDataLikelihoodDelegate();
        final int stateCount = delegate.getPatternList().getDataType().getStateCount();
        final int patternCount = delegate.getPatternList().getPatternCount();
        final int categoryCount = delegate.getSiteRateModel().getCategoryCount();
        final double[] partials = new double[stateCount * patternCount * categoryCount];
        delegate.getPartials(likelihood.getTree().getRoot().getNumber(), partials);

        final double[] frequencies = delegate.getEvolutionaryProcessDelegate().getRootStateFrequencies();
        final double[] categoryWeights = delegate.getSiteRateModel().getCategoryProportions();

        double minPartial = Double.POSITIVE_INFINITY;
        double maxPartial = Double.NEGATIVE_INFINITY;
        int negativePartialCount = 0;
        int nanPartialCount = 0;
        for (double partial : partials) {
            if (Double.isNaN(partial)) {
                nanPartialCount++;
            } else {
                minPartial = Math.min(minPartial, partial);
                maxPartial = Math.max(maxPartial, partial);
                if (partial < 0.0) {
                    negativePartialCount++;
                }
            }
        }

        double minPatternSum = Double.POSITIVE_INFINITY;
        double maxPatternSum = Double.NEGATIVE_INFINITY;
        int nonPositivePatternCount = 0;
        int nanPatternCount = 0;
        for (int p = 0; p < patternCount; p++) {
            double patternSum = 0.0;
            for (int c = 0; c < categoryCount; c++) {
                final double categoryWeight = categoryWeights[c];
                final int offset = ((c * patternCount) + p) * stateCount;
                double categorySum = 0.0;
                for (int s = 0; s < stateCount; s++) {
                    categorySum += frequencies[s] * partials[offset + s];
                }
                patternSum += categoryWeight * categorySum;
            }
            if (Double.isNaN(patternSum)) {
                nanPatternCount++;
            } else {
                minPatternSum = Math.min(minPatternSum, patternSum);
                maxPatternSum = Math.max(maxPatternSum, patternSum);
                if (!(patternSum > 0.0)) {
                    nonPositivePatternCount++;
                }
            }
        }

        System.err.println("rootDiagnostics." + label +
                ".partialRange=[" + minPartial + ", " + maxPartial + "]" +
                ", negativePartials=" + negativePartialCount +
                ", nanPartials=" + nanPartialCount +
                ", patternSumRange=[" + minPatternSum + ", " + maxPatternSum + "]" +
                ", nonPositivePatterns=" + nonPositivePatternCount +
                ", nanPatterns=" + nanPatternCount);
    }

    private static TreeModel createBalancedTree(final int tipCount) {
        final SimpleNode root = buildBalancedSubtree(0, tipCount, heightForTipCount(tipCount));
        final SimpleTree simpleTree = new SimpleTree(root);
        return new DefaultTreeModel("fiftyTipBalancedTree", simpleTree);
    }

    private static SimpleNode buildBalancedSubtree(final int fromInclusive,
                                                   final int toExclusive,
                                                   final double height) {
        final int count = toExclusive - fromInclusive;
        final SimpleNode node = new SimpleNode();
        node.setHeight(height);

        if (count == 1) {
            node.setTaxon(new Taxon("taxon" + String.format("%02d", fromInclusive + 1)));
            node.setHeight(0.0);
            return node;
        }

        final int leftCount = count / 2;
        final int split = fromInclusive + leftCount;
        node.addChild(buildBalancedSubtree(fromInclusive, split, heightForTipCount(leftCount)));
        node.addChild(buildBalancedSubtree(split, toExclusive, heightForTipCount(count - leftCount)));
        return node;
    }

    private static double heightForTipCount(final int tipCount) {
        int levels = 0;
        int n = 1;
        while (n < tipCount) {
            n <<= 1;
            levels++;
        }
        return levels * BRANCH_HEIGHT_UNIT;
    }

    private static LogRateSubstitutionModel createIndependentRewardSubstitutionModel(final String id) {
        final FrequencyModel frequencyModel =
                new FrequencyModel(Nucleotides.INSTANCE, INDEPENDENT_ROOT_FREQUENCIES);
        final Parameter logRates = new Parameter.Default(
                id + ".logRates",
                new double[]{
                        Math.log(0.10), Math.log(1.00), Math.log(1.00),
                        Math.log(1.00), Math.log(1.00), Math.log(0.80),
                        Math.log(0.10), Math.log(1.0e-5), Math.log(1.0e-5),
                        Math.log(1.0e-5), Math.log(1.0e-5), Math.log(0.80)
                }
        );
        final LogRateSubstitutionModel model = new LogRateSubstitutionModel(
                id,
                Nucleotides.INSTANCE,
                frequencyModel,
                new LogAdditiveCtmcRateProvider.DataAugmented.Basic(id + ".rateProvider", logRates)
        );
        model.setNormalization(true);
        model.setScaleRatesByFrequencies(false);
        return model;
    }

    private static SubstitutionModel createDependentSubstitutionModel(final int index) {
        final double[][] frequencies = new double[][]{
                {0.30, 0.20, 0.25, 0.25},
                {0.20, 0.30, 0.20, 0.30},
                {0.25, 0.25, 0.30, 0.20}
        };
        final double[] kappas = new double[]{2.0, 3.5, 1.25};
        final FrequencyModel frequencyModel =
                new FrequencyModel(Nucleotides.INSTANCE, frequencies[index % frequencies.length]);
        return new HKY(new Parameter.Default("dependent.kappa." + (index + 1), kappas[index % kappas.length]),
                frequencyModel);
    }

    private static RewardRateHmcKernel createRewardRateHmcKernel(
            final TreeDataLikelihood independentLikelihood,
            final TreeDataLikelihood[] dependentLikelihoods,
            final RewardsAwareBranchModel rewardsAwareBranchModel,
            final RewardsAwareMixtureBranchRates branchRates,
            final RewardRates rewardRates,
            final double lower,
            final double upper) {

        final List<GradientWrtParameterProvider> gradients =
                new ArrayList<GradientWrtParameterProvider>();
        gradients.add(new RewardsAwareBranchModelRewardRateGradient(
                independentLikelihood,
                rewardsAwareBranchModel,
                rewardRates.getValues(),
                rewardRates.getVaryingValues(),
                null));
        for (TreeDataLikelihood dependentLikelihood : dependentLikelihoods) {
            gradients.add(new RewardsAwareMixtureBranchRatesRewardRateGradient(
                    dependentLikelihood,
                    branchRates,
                    rewardRates.getValues(),
                    rewardRates.getVaryingValues(),
                    null));
        }

        return new RewardRateHmcKernel(
                rewardRates.getVaryingValues(),
                gradients,
                new LogTargetEvaluator() {
                    @Override
                    public double evaluate() {
                        return logTarget(independentLikelihood, dependentLikelihoods);
                    }
                },
                lower,
                upper,
                REWARD_RATE_HMC_STEP_SIZE,
                REWARD_RATE_HMC_LEAPFROG_STEPS);
    }

    private static RewardRates createFixedRewardRates(final String id, final double[] values) {
        return new RewardRates(
                new Parameter.Default(id + ".values", Arrays.copyOf(values, values.length)),
                null,
                new Parameter.Default(id + ".internal", Arrays.copyOfRange(values, 2, values.length)),
                new Parameter.Default(id + ".mapping", new double[]{0.0, 1.0, 2.0, 3.0})
        );
    }

    private static RewardRates createInferenceRewardRates(final String id,
                                                          final double fixedLower,
                                                          final double fixedUpper) {
        final Parameter fixedValues =
                new Parameter.Default(id + ".fixed", new double[]{fixedLower, fixedUpper});
        final Parameter varyingValues =
                new Parameter.Default(id + ".internal",
                        Arrays.copyOf(INITIAL_REWARD_RATE_VARYING_VALUES,
                                INITIAL_REWARD_RATE_VARYING_VALUES.length));
        varyingValues.addBounds(new Parameter.DefaultBounds(
                fixedUpper,
                fixedLower,
                varyingValues.getDimension()));

        final CompoundParameter values = new CompoundParameter(id + ".values");
        values.addParameter(fixedValues);
        values.addParameter(varyingValues);

        return new RewardRates(
                values,
                fixedValues,
                varyingValues,
                new Parameter.Default(id + ".mapping", new double[]{0.0, 1.0, 2.0, 3.0})
        );
    }

    private static RewardsAwareMixtureBranchRates createRewardBranchRates(final TreeModel tree,
                                                                          final double[] cts,
                                                                          final double[] indicators,
                                                                          final double[] atoms,
                                                                          final RewardRates rewardRates,
                                                                          final String prefix) {
        return new RewardsAwareMixtureBranchRates(
                tree,
                new Parameter.Default(prefix + ".ctsRewards", Arrays.copyOf(cts, cts.length)),
                new Parameter.Default(prefix + ".indicator", Arrays.copyOf(indicators, indicators.length)),
                new Parameter.Default(prefix + ".atomIndices", Arrays.copyOf(atoms, atoms.length)),
                rewardRates,
                new ArbitraryBranchRates.BranchRateTransform.None(),
                false,
                TreeParameterModel.Type.WITHOUT_ROOT
        );
    }

    private static SimpleAlignment createIndependentTipAlignment(final TreeModel tree,
                                                                final int[] generatedTipStates) {
        final SimpleAlignment alignment = new SimpleAlignment();
        alignment.setId("independent.generatedTips");
        alignment.setDataType(Nucleotides.INSTANCE);
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            final NodeRef tip = tree.getExternalNode(i);
            final Taxon taxon = tree.getNodeTaxon(tip);
            final int state = generatedTipStates[i];
            final Sequence sequence = new Sequence(taxon, Character.toString(NUCLEOTIDE_CODES[state]));
            sequence.setDataType(Nucleotides.INSTANCE);
            alignment.addSequence(sequence);
        }
        return alignment;
    }

    private static SimpleAlignment simulateDependentAlignment(final TreeModel tree,
                                                              final RewardsAwareMixtureBranchRates branchRates,
                                                              final SubstitutionModel substitutionModel,
                                                              final GammaSiteRateModel siteRateModel,
                                                              final int siteCount,
                                                              final String id) {
        final ArrayList<Partition> partitions = new ArrayList<Partition>();
        partitions.add(new Partition(
                tree,
                new HomogeneousBranchModel(substitutionModel),
                siteRateModel,
                branchRates,
                substitutionModel.getFrequencyModel(),
                0,
                siteCount - 1,
                1,
                Nucleotides.INSTANCE
        ));

        final SimpleAlignment alignment = new BeagleSequenceSimulator(partitions).simulate(false, false);
        alignment.setId(id + ".alignment");
        alignment.setDataType(Nucleotides.INSTANCE);
        return alignment;
    }

    private static void writeSummary(final File file,
                                     final TreeModel tree,
                                     final RewardsAwareCtmcSimulation simulation,
                                     final int dependentSiteCount,
                                     final int sweeps,
                                     final int burnIn,
                                     final RunResult independentOnly,
                                     final RunResult withDependents) throws IOException {
        try (PrintWriter out = new PrintWriter(file)) {
            out.println("metric\tvalue");
            out.println("seed\t" + SEED);
            out.println("tipCount\t" + tree.getExternalNodeCount());
            out.println("branchCount\t" + (tree.getNodeCount() - 1));
            out.println("dependentCount\t" + DEPENDENT_COUNT);
            out.println("dependentSiteCount\t" + dependentSiteCount);
            out.println("sweeps\t" + sweeps);
            out.println("burnIn\t" + burnIn);
            out.println("simulation.totalJumps\t" + simulation.getSummary().totalJumps);
            out.println("simulation.atomicBranchCount\t" + simulation.getSummary().atomicBranchCount);
            out.println("simulation.continuousBranchCount\t" + simulation.getSummary().continuousBranchCount);
            independentOnly.writeMetrics(out);
            withDependents.writeMetrics(out);
        }
    }

    private static void writeBranchSummary(final File file,
                                           final TreeModel tree,
                                           final RewardsAwareCtmcSimulation simulation,
                                           final double[] trueCts,
                                           final int[] trueIndicators,
                                           final int[] trueAtoms,
                                           final RunResult independentOnly,
                                           final RunResult withDependents) throws IOException {
        try (PrintWriter out = new PrintWriter(file)) {
            out.println("branchIndex\tnodeNumber\ttrueIndicator\ttrueAtom\ttrueCtsReward\tindependentOnly.pAtomic\tindependentOnly.atomCall\tindependentOnly.indicatorESS\twithDependents.pAtomic\twithDependents.atomCall\twithDependents.indicatorESS");
            for (int branch = 0; branch < trueCts.length; branch++) {
                out.println(branch + "\t" +
                        simulation.getNodeNumberForBranchIndex(branch) + "\t" +
                        trueIndicators[branch] + "\t" +
                        trueAtoms[branch] + "\t" +
                        trueCts[branch] + "\t" +
                        independentOnly.posteriorAtomicProb[branch] + "\t" +
                        independentOnly.atomCalls[branch] + "\t" +
                        independentOnly.indicatorEss[branch] + "\t" +
                        withDependents.posteriorAtomicProb[branch] + "\t" +
                        withDependents.atomCalls[branch] + "\t" +
                        withDependents.indicatorEss[branch]);
            }
        }
    }

    private static void writeNewick(final File file, final TreeModel tree) throws IOException {
        try (PrintWriter out = new PrintWriter(file)) {
            out.println(toNewick(tree, tree.getRoot()) + ";");
        }
    }

    private static String toNewick(final Tree tree, final NodeRef node) {
        final StringBuilder sb = new StringBuilder();
        if (tree.isExternal(node)) {
            sb.append(tree.getNodeTaxon(node).getId());
        } else {
            sb.append('(');
            for (int i = 0; i < tree.getChildCount(node); i++) {
                if (i > 0) {
                    sb.append(',');
                }
                sb.append(toNewick(tree, tree.getChild(node, i)));
            }
            sb.append(')');
        }
        if (!tree.isRoot(node)) {
            sb.append(':').append(tree.getBranchLength(node));
        }
        return sb.toString();
    }

    private static void writeAlignment(final File file, final SimpleAlignment alignment) throws IOException {
        try (PrintWriter out = new PrintWriter(file)) {
            for (int i = 0; i < alignment.getSequenceCount(); i++) {
                final Sequence sequence = alignment.getSequence(i);
                out.println('>' + sequence.getTaxon().getId());
                out.println(sequence.getSequenceString());
            }
        }
    }

    private static double[] toDouble(final int[] values) {
        final double[] out = new double[values.length];
        for (int i = 0; i < values.length; i++) {
            out[i] = values[i];
        }
        return out;
    }

    private static int argMax(final int[] counts) {
        int bestIndex = 0;
        int bestValue = counts[0];
        for (int i = 1; i < counts.length; i++) {
            if (counts[i] > bestValue) {
                bestValue = counts[i];
                bestIndex = i;
            }
        }
        return bestIndex;
    }

    private static double effectiveSampleSize(final double[] values) {
        final int n = values.length;
        if (n <= 1) {
            return n;
        }

        double mean = 0.0;
        for (double value : values) {
            mean += value;
        }
        mean /= n;

        double variance = 0.0;
        for (double value : values) {
            final double centered = value - mean;
            variance += centered * centered;
        }
        variance /= n;
        if (!(variance > 0.0) || Double.isNaN(variance)) {
            return n;
        }

        double autocorrelationSum = 0.0;
        for (int lag = 1; lag < n; lag++) {
            double covariance = 0.0;
            for (int i = 0; i < n - lag; i++) {
                covariance += (values[i] - mean) * (values[i + lag] - mean);
            }
            covariance /= (n - lag);
            final double autocorrelation = covariance / variance;
            if (!(autocorrelation > 0.0)) {
                break;
            }
            autocorrelationSum += autocorrelation;
        }

        final double ess = n / (1.0 + 2.0 * autocorrelationSum);
        return Math.max(1.0, Math.min(n, ess));
    }

    private static double mean(final double[] values) {
        if (values.length == 0) {
            return Double.NaN;
        }
        double sum = 0.0;
        for (double value : values) {
            sum += value;
        }
        return sum / values.length;
    }

    private static double median(final double[] values) {
        if (values.length == 0) {
            return Double.NaN;
        }
        final double[] copy = Arrays.copyOf(values, values.length);
        Arrays.sort(copy);
        final int midpoint = copy.length / 2;
        if ((copy.length & 1) == 1) {
            return copy[midpoint];
        }
        return 0.5 * (copy[midpoint - 1] + copy[midpoint]);
    }

    private static double min(final double[] values) {
        double min = Double.POSITIVE_INFINITY;
        for (double value : values) {
            min = Math.min(min, value);
        }
        return min;
    }

    private static double max(final double[] values) {
        double max = Double.NEGATIVE_INFINITY;
        for (double value : values) {
            max = Math.max(max, value);
        }
        return max;
    }

    private static int countAtLeast(final double[] values, final double threshold) {
        int count = 0;
        for (double value : values) {
            if (value >= threshold) {
                count++;
            }
        }
        return count;
    }

    private interface LogTargetEvaluator {
        double evaluate();
    }

    private static final class RewardRateHmcKernel {

        private static final double TRANSFORM_EPSILON = 1.0e-12;

        private final Parameter rewardRateVaryingValues;
        private final List<GradientWrtParameterProvider> gradients;
        private final LogTargetEvaluator target;
        private final double lower;
        private final double upper;
        private final double width;
        private final double stepSize;
        private final int leapfrogSteps;

        private int attemptCount = 0;
        private int acceptCount = 0;

        private RewardRateHmcKernel(final Parameter rewardRateVaryingValues,
                                    final List<GradientWrtParameterProvider> gradients,
                                    final LogTargetEvaluator target,
                                    final double lower,
                                    final double upper,
                                    final double stepSize,
                                    final int leapfrogSteps) {
            if (!(upper > lower)) {
                throw new IllegalArgumentException("Reward-rate HMC bounds must satisfy upper > lower.");
            }
            if (!(stepSize > 0.0)) {
                throw new IllegalArgumentException("Reward-rate HMC step size must be positive.");
            }
            if (leapfrogSteps <= 0) {
                throw new IllegalArgumentException("Reward-rate HMC leapfrog steps must be positive.");
            }
            this.rewardRateVaryingValues = rewardRateVaryingValues;
            this.gradients = gradients;
            this.target = target;
            this.lower = lower;
            this.upper = upper;
            this.width = upper - lower;
            this.stepSize = stepSize;
            this.leapfrogSteps = leapfrogSteps;
        }

        private void doOperation() {
            attemptCount++;

            final double[] oldValues = rewardRateVaryingValues.getParameterValues();
            final double[] theta = valuesToTheta(oldValues);
            final double oldLogDensity = transformedLogDensity(theta);
            if (!Double.isFinite(oldLogDensity)) {
                restoreValues(oldValues);
                return;
            }

            final double[] momentum = new double[theta.length];
            for (int i = 0; i < momentum.length; i++) {
                momentum[i] = MathUtils.nextGaussian();
            }
            final double oldKineticEnergy = kineticEnergy(momentum);

            double[] gradient = transformedGradient(theta);
            if (!allFinite(gradient)) {
                restoreValues(oldValues);
                return;
            }

            for (int i = 0; i < momentum.length; i++) {
                momentum[i] += 0.5 * stepSize * gradient[i];
            }

            for (int step = 0; step < leapfrogSteps; step++) {
                for (int i = 0; i < theta.length; i++) {
                    theta[i] += stepSize * momentum[i];
                }

                gradient = transformedGradient(theta);
                if (!allFinite(gradient)) {
                    restoreValues(oldValues);
                    return;
                }

                final double scale = step == leapfrogSteps - 1 ? 0.5 : 1.0;
                for (int i = 0; i < momentum.length; i++) {
                    momentum[i] += scale * stepSize * gradient[i];
                }
            }

            final double newLogDensity = transformedLogDensity(theta);
            final double newKineticEnergy = kineticEnergy(momentum);
            final double logAcceptanceProbability =
                    newLogDensity - oldLogDensity + oldKineticEnergy - newKineticEnergy;

            if (Double.isFinite(logAcceptanceProbability) &&
                    Math.log(MathUtils.nextDouble()) < logAcceptanceProbability) {
                acceptCount++;
            } else {
                restoreValues(oldValues);
            }
        }

        private double transformedLogDensity(final double[] theta) {
            setValuesFromTheta(theta);
            final double logLikelihood = target.evaluate();
            if (!Double.isFinite(logLikelihood)) {
                return Double.NEGATIVE_INFINITY;
            }

            double logJacobian = 0.0;
            for (double value : theta) {
                final double s = logistic(value);
                logJacobian += Math.log(width) + Math.log(s) + Math.log1p(-s);
            }
            return logLikelihood + logJacobian;
        }

        private double[] transformedGradient(final double[] theta) {
            setValuesFromTheta(theta);
            final double logLikelihood = target.evaluate();
            if (!Double.isFinite(logLikelihood)) {
                final double[] invalid = new double[theta.length];
                Arrays.fill(invalid, Double.NaN);
                return invalid;
            }

            final double[] valueGradient = new double[theta.length];
            for (GradientWrtParameterProvider gradientProvider : gradients) {
                final double[] contribution = gradientProvider.getGradientLogDensity();
                if (contribution.length != valueGradient.length) {
                    throw new IllegalStateException(
                            "Reward-rate gradient dimension mismatch: expected " +
                                    valueGradient.length + " but found " + contribution.length);
                }
                for (int i = 0; i < valueGradient.length; i++) {
                    valueGradient[i] += contribution[i];
                }
            }

            final double[] thetaGradient = new double[theta.length];
            for (int i = 0; i < theta.length; i++) {
                final double s = logistic(theta[i]);
                final double jacobian = width * s * (1.0 - s);
                thetaGradient[i] = valueGradient[i] * jacobian + 1.0 - 2.0 * s;
            }
            return thetaGradient;
        }

        private void setValuesFromTheta(final double[] theta) {
            for (int i = 0; i < theta.length; i++) {
                final double s = logistic(theta[i]);
                rewardRateVaryingValues.setParameterValueQuietly(i, lower + width * s);
            }
            rewardRateVaryingValues.fireParameterChangedEvent();
        }

        private double[] valuesToTheta(final double[] values) {
            final double[] theta = new double[values.length];
            for (int i = 0; i < values.length; i++) {
                double s = (values[i] - lower) / width;
                s = Math.max(TRANSFORM_EPSILON, Math.min(1.0 - TRANSFORM_EPSILON, s));
                theta[i] = Math.log(s) - Math.log1p(-s);
            }
            return theta;
        }

        private void restoreValues(final double[] values) {
            for (int i = 0; i < values.length; i++) {
                rewardRateVaryingValues.setParameterValueQuietly(i, values[i]);
            }
            rewardRateVaryingValues.fireParameterChangedEvent();
        }

        private static double logistic(final double x) {
            if (x >= 0.0) {
                final double z = Math.exp(-x);
                return 1.0 / (1.0 + z);
            } else {
                final double z = Math.exp(x);
                return z / (1.0 + z);
            }
        }

        private static double kineticEnergy(final double[] momentum) {
            double sum = 0.0;
            for (double value : momentum) {
                sum += value * value;
            }
            return 0.5 * sum;
        }

        private static boolean allFinite(final double[] values) {
            for (double value : values) {
                if (!Double.isFinite(value)) {
                    return false;
                }
            }
            return true;
        }

        private int getAttemptCount() {
            return attemptCount;
        }

        private int getAcceptCount() {
            return acceptCount;
        }
    }

    private static final class RunResult {
        final String label;
        final int sampleCount;
        final double initialLogTarget;
        final double finalLogTarget;
        final double[] posteriorAtomicProb;
        final int[] atomCalls;
        final double indicatorAccuracy;
        final double atomicRecall;
        final double atomicPrecision;
        final double atomAccuracyAmongTrueAtomic;
        final double confidentAtomAccuracyAmongTrueAtomic;
        final double indicatorBrier;
        final double meanAtomicProbTrueAtomic;
        final double meanAtomicProbTrueContinuous;
        final int calledAtomicCount;
        final double[] indicatorEss;
        final double indicatorEssMean;
        final double indicatorEssMedian;
        final double indicatorEssMin;
        final double indicatorEssMax;
        final int indicatorEssNearIndependentCount;
        final double logTargetEss;
        final double[] trueRewardRateVaryingValues;
        final double[] initialRewardRateVaryingValues;
        final double[] rewardRatePosteriorMeans;
        final double[] finalRewardRateVaryingValues;
        final double[] rewardRateEss;
        final int rewardRateHmcAttemptCount;
        final int rewardRateHmcAcceptCount;

        private RunResult(final String label,
                          final int sampleCount,
                          final double initialLogTarget,
                          final double finalLogTarget,
                          final double[] posteriorAtomicProb,
                          final int[] atomCalls,
                          final double indicatorAccuracy,
                          final double atomicRecall,
                          final double atomicPrecision,
                          final double atomAccuracyAmongTrueAtomic,
                          final double confidentAtomAccuracyAmongTrueAtomic,
                          final double indicatorBrier,
                          final double meanAtomicProbTrueAtomic,
                          final double meanAtomicProbTrueContinuous,
                          final int calledAtomicCount,
                          final double[] indicatorEss,
                          final double indicatorEssMean,
                          final double indicatorEssMedian,
                          final double indicatorEssMin,
                          final double indicatorEssMax,
                          final int indicatorEssNearIndependentCount,
                          final double logTargetEss,
                          final double[] trueRewardRateVaryingValues,
                          final double[] initialRewardRateVaryingValues,
                          final double[] rewardRatePosteriorMeans,
                          final double[] finalRewardRateVaryingValues,
                          final double[] rewardRateEss,
                          final int rewardRateHmcAttemptCount,
                          final int rewardRateHmcAcceptCount) {
            this.label = label;
            this.sampleCount = sampleCount;
            this.initialLogTarget = initialLogTarget;
            this.finalLogTarget = finalLogTarget;
            this.posteriorAtomicProb = Arrays.copyOf(posteriorAtomicProb, posteriorAtomicProb.length);
            this.atomCalls = Arrays.copyOf(atomCalls, atomCalls.length);
            this.indicatorAccuracy = indicatorAccuracy;
            this.atomicRecall = atomicRecall;
            this.atomicPrecision = atomicPrecision;
            this.atomAccuracyAmongTrueAtomic = atomAccuracyAmongTrueAtomic;
            this.confidentAtomAccuracyAmongTrueAtomic = confidentAtomAccuracyAmongTrueAtomic;
            this.indicatorBrier = indicatorBrier;
            this.meanAtomicProbTrueAtomic = meanAtomicProbTrueAtomic;
            this.meanAtomicProbTrueContinuous = meanAtomicProbTrueContinuous;
            this.calledAtomicCount = calledAtomicCount;
            this.indicatorEss = Arrays.copyOf(indicatorEss, indicatorEss.length);
            this.indicatorEssMean = indicatorEssMean;
            this.indicatorEssMedian = indicatorEssMedian;
            this.indicatorEssMin = indicatorEssMin;
            this.indicatorEssMax = indicatorEssMax;
            this.indicatorEssNearIndependentCount = indicatorEssNearIndependentCount;
            this.logTargetEss = logTargetEss;
            this.trueRewardRateVaryingValues =
                    Arrays.copyOf(trueRewardRateVaryingValues, trueRewardRateVaryingValues.length);
            this.initialRewardRateVaryingValues =
                    Arrays.copyOf(initialRewardRateVaryingValues, initialRewardRateVaryingValues.length);
            this.rewardRatePosteriorMeans =
                    Arrays.copyOf(rewardRatePosteriorMeans, rewardRatePosteriorMeans.length);
            this.finalRewardRateVaryingValues =
                    Arrays.copyOf(finalRewardRateVaryingValues, finalRewardRateVaryingValues.length);
            this.rewardRateEss = Arrays.copyOf(rewardRateEss, rewardRateEss.length);
            this.rewardRateHmcAttemptCount = rewardRateHmcAttemptCount;
            this.rewardRateHmcAcceptCount = rewardRateHmcAcceptCount;
        }

        private void writeMetrics(final PrintWriter out) {
            out.println(label + ".sampleCount\t" + sampleCount);
            out.println(label + ".initialLogTarget\t" + initialLogTarget);
            out.println(label + ".finalLogTarget\t" + finalLogTarget);
            out.println(label + ".logTargetESS\t" + logTargetEss);
            out.println(label + ".indicatorESS.mean\t" + indicatorEssMean);
            out.println(label + ".indicatorESS.median\t" + indicatorEssMedian);
            out.println(label + ".indicatorESS.min\t" + indicatorEssMin);
            out.println(label + ".indicatorESS.max\t" + indicatorEssMax);
            out.println(label + ".indicatorESS.nearIndependentCount\t" +
                    indicatorEssNearIndependentCount);
            out.println(label + ".indicatorAccuracy\t" + indicatorAccuracy);
            out.println(label + ".atomicRecall\t" + atomicRecall);
            out.println(label + ".atomicPrecision\t" + atomicPrecision);
            out.println(label + ".calledAtomicCount\t" + calledAtomicCount);
            out.println(label + ".atomAccuracyAmongTrueAtomic\t" + atomAccuracyAmongTrueAtomic);
            out.println(label + ".confidentAtomAccuracyAmongTrueAtomic\t" +
                    confidentAtomAccuracyAmongTrueAtomic);
            out.println(label + ".indicatorBrier\t" + indicatorBrier);
            out.println(label + ".meanAtomicProbTrueAtomic\t" + meanAtomicProbTrueAtomic);
            out.println(label + ".meanAtomicProbTrueContinuous\t" + meanAtomicProbTrueContinuous);
            out.println(label + ".rewardRateHMC.attemptCount\t" + rewardRateHmcAttemptCount);
            out.println(label + ".rewardRateHMC.acceptCount\t" + rewardRateHmcAcceptCount);
            out.println(label + ".rewardRateHMC.acceptanceRate\t" +
                    (rewardRateHmcAttemptCount == 0
                            ? Double.NaN
                            : rewardRateHmcAcceptCount / (double) rewardRateHmcAttemptCount));
            for (int i = 0; i < rewardRatePosteriorMeans.length; i++) {
                out.println(label + ".rewardRateInternal." + i + ".true\t" +
                        trueRewardRateVaryingValues[i]);
                out.println(label + ".rewardRateInternal." + i + ".initial\t" +
                        initialRewardRateVaryingValues[i]);
                out.println(label + ".rewardRateInternal." + i + ".posteriorMean\t" +
                        rewardRatePosteriorMeans[i]);
                out.println(label + ".rewardRateInternal." + i + ".final\t" +
                        finalRewardRateVaryingValues[i]);
                out.println(label + ".rewardRateInternal." + i + ".ESS\t" + rewardRateEss[i]);
            }
        }

        private String toConsoleString() {
            return label +
                    ": indicatorAccuracy=" + indicatorAccuracy +
                    ", atomicRecall=" + atomicRecall +
                    ", atomicPrecision=" + atomicPrecision +
                    ", atomAccuracyAmongTrueAtomic=" + atomAccuracyAmongTrueAtomic +
                    ", brier=" + indicatorBrier +
                    ", logTargetESS=" + logTargetEss +
                    ", indicatorESS[min/median/mean]=" + indicatorEssMin + "/" +
                    indicatorEssMedian + "/" + indicatorEssMean +
                    ", rewardRateHMC.acceptance=" +
                    (rewardRateHmcAttemptCount == 0
                            ? Double.NaN
                            : rewardRateHmcAcceptCount / (double) rewardRateHmcAttemptCount) +
                    ", rewardRateMeans=" + Arrays.toString(rewardRatePosteriorMeans) +
                    ", finalLogTarget=" + finalLogTarget;
        }
    }
}

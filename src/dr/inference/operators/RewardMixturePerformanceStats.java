package dr.inference.operators;

import java.util.concurrent.atomic.AtomicLong;

/**
 * Optional reward-mixture performance counters.
 *
 * Enable with:
 *   -Ddr.inference.operators.rewardMixturePerf=true
 */
public final class RewardMixturePerformanceStats {

    public static final String ENABLE_PROPERTY = "dr.inference.operators.rewardMixturePerf";
    public static final boolean ENABLED = Boolean.parseBoolean(System.getProperty(ENABLE_PROPERTY, "false"));

    public enum OperationCacheClearReason {
        UNKNOWN,
        OPERATION_START,
        SMOOTH_POSITION_UPDATE,
        ACCEPTED_CATEGORY_CROSSING,
        GIBBS_CATEGORY_CHANGE,
        FINAL_POSITION_REFRESH,
        DYNAMIC_EMBEDDING_REFRESH
    }

    private static final AtomicLong rewardProviderRefreshCalls = new AtomicLong();
    private static final AtomicLong rewardProviderRefreshNanos = new AtomicLong();
    private static final AtomicLong postOrderRefreshCalls = new AtomicLong();
    private static final AtomicLong postOrderRefreshNanos = new AtomicLong();
    private static final AtomicLong preOrderRefreshCalls = new AtomicLong();
    private static final AtomicLong preOrderRefreshNanos = new AtomicLong();
    private static final AtomicLong discretePreOrderEnsureHits = new AtomicLong();
    private static final AtomicLong discretePreOrderEnsureHitNanos = new AtomicLong();
    private static final AtomicLong discretePreOrderEnsureMisses = new AtomicLong();
    private static final AtomicLong discretePreOrderEnsureMissNanos = new AtomicLong();
    private static final AtomicLong dependentPrepareCalls = new AtomicLong();
    private static final AtomicLong dependentPrepareNanos = new AtomicLong();
    private static final AtomicLong discreteLikelihoodTraversalCalls = new AtomicLong();
    private static final AtomicLong discreteLikelihoodTraversalStatsOnlyCalls = new AtomicLong();
    private static final AtomicLong discreteLikelihoodTraversalFullTreeCalls = new AtomicLong();
    private static final AtomicLong discreteLikelihoodTraversalPartialTreeCalls = new AtomicLong();
    private static final AtomicLong discreteLikelihoodTraversalEmptyCalls = new AtomicLong();
    private static final AtomicLong discreteLikelihoodTraversalBranchOps = new AtomicLong();
    private static final AtomicLong discreteLikelihoodTraversalNodeOps = new AtomicLong();
    private static final AtomicLong discreteLikelihoodTraversalSubstitutionDirtyCalls = new AtomicLong();
    private static final AtomicLong discreteLikelihoodTraversalSiteDirtyCalls = new AtomicLong();
    private static final AtomicLong discreteLikelihoodTraversalRootFrequencyDirtyCalls = new AtomicLong();

    private static final AtomicLong operationCacheBegins = new AtomicLong();
    private static final AtomicLong operationCacheClears = new AtomicLong();
    private static final AtomicLong operationCacheClearsUnknown = new AtomicLong();
    private static final AtomicLong operationCacheClearsOperationStart = new AtomicLong();
    private static final AtomicLong operationCacheClearsSmoothPositionUpdate = new AtomicLong();
    private static final AtomicLong operationCacheClearsAcceptedCategoryCrossing = new AtomicLong();
    private static final AtomicLong operationCacheClearsGibbsCategoryChange = new AtomicLong();
    private static final AtomicLong operationCacheClearsFinalPositionRefresh = new AtomicLong();
    private static final AtomicLong operationCacheClearsDynamicEmbeddingRefresh = new AtomicLong();
    private static final AtomicLong branchWeightCacheHits = new AtomicLong();
    private static final AtomicLong branchWeightCacheMisses = new AtomicLong();
    private static final AtomicLong branchWeightCacheBypasses = new AtomicLong();
    private static final AtomicLong likelihoodMessageMisses = new AtomicLong();
    private static final AtomicLong branchWeightComputations = new AtomicLong();
    private static final AtomicLong branchWeightComputationNanos = new AtomicLong();
    private static final AtomicLong branchMessageLoadCalls = new AtomicLong();
    private static final AtomicLong branchMessageLoadNanos = new AtomicLong();
    private static final AtomicLong branchMessageCacheFillCalls = new AtomicLong();
    private static final AtomicLong branchMessageCacheFillNanos = new AtomicLong();
    private static final AtomicLong branchPrePartialLoadCalls = new AtomicLong();
    private static final AtomicLong branchPrePartialLoadNanos = new AtomicLong();
    private static final AtomicLong branchPostPartialLoadCalls = new AtomicLong();
    private static final AtomicLong branchPostPartialLoadNanos = new AtomicLong();
    private static final AtomicLong branchPreScaleLoadCalls = new AtomicLong();
    private static final AtomicLong branchPreScaleLoadNanos = new AtomicLong();
    private static final AtomicLong branchPostScaleLoadCalls = new AtomicLong();
    private static final AtomicLong branchPostScaleLoadNanos = new AtomicLong();
    private static final AtomicLong branchAtomicWeightLoopCalls = new AtomicLong();
    private static final AtomicLong branchAtomicWeightLoopNanos = new AtomicLong();
    private static final AtomicLong branchContinuousMatrixAccessCalls = new AtomicLong();
    private static final AtomicLong branchContinuousMatrixAccessNanos = new AtomicLong();
    private static final AtomicLong branchContinuousWeightCalls = new AtomicLong();
    private static final AtomicLong branchContinuousWeightNanos = new AtomicLong();
    private static final AtomicLong dependentEvidenceCalls = new AtomicLong();
    private static final AtomicLong dependentEvidenceNanos = new AtomicLong();

    private static final AtomicLong beaglePrepareCalls = new AtomicLong();
    private static final AtomicLong beaglePrepareNanos = new AtomicLong();
    private static final AtomicLong beagleRefreshSpectralNanos = new AtomicLong();
    private static final AtomicLong beagleFillPostNanos = new AtomicLong();
    private static final AtomicLong beagleFillTopNanos = new AtomicLong();
    private static final AtomicLong beagleFillPostTopNanos = new AtomicLong();
    private static final AtomicLong beagleFillPreBottomNanos = new AtomicLong();
    private static final AtomicLong beagleFillPreOrderNanos = new AtomicLong();
    private static final AtomicLong beagleLogEvidenceCalls = new AtomicLong();
    private static final AtomicLong beagleLogEvidenceNanos = new AtomicLong();
    private static final AtomicLong beagleLogEvidenceAtomicCalls = new AtomicLong();
    private static final AtomicLong beagleLogEvidenceContinuousCalls = new AtomicLong();
    private static final AtomicLong beagleCachedMessageEvidenceCalls = new AtomicLong();
    private static final AtomicLong beagleCachedMessageEvidenceNanos = new AtomicLong();
    private static final AtomicLong beagleEdgeInnerProductCalls = new AtomicLong();
    private static final AtomicLong beagleEdgeInnerProductNanos = new AtomicLong();
    private static final AtomicLong beagleMessageCopyCalls = new AtomicLong();
    private static final AtomicLong beagleMessageCopyNanos = new AtomicLong();

    private static final AtomicLong discontinuousSteps = new AtomicLong();
    private static final AtomicLong discontinuousAcceptedCrossings = new AtomicLong();
    private static final AtomicLong discontinuousReflections = new AtomicLong();
    private static final AtomicLong cacheClearsAfterAcceptedCrossing = new AtomicLong();
    private static final AtomicLong skippedCacheClearsAfterNoCrossing = new AtomicLong();
    private static final AtomicLong categoricalGibbsBranchUpdates = new AtomicLong();
    private static final AtomicLong categoricalGibbsCategoryChanges = new AtomicLong();
    private static final AtomicLong categoricalGibbsSameCategoryUpdates = new AtomicLong();
    private static final AtomicLong categoricalGibbsSkippedCacheClearsAfterSameCategory = new AtomicLong();

    static {
        if (ENABLED) {
            Runtime.getRuntime().addShutdownHook(new Thread(new Runnable() {
                @Override
                public void run() {
                    printReport();
                }
            }, "reward-mixture-performance-report"));
        }
    }

    private RewardMixturePerformanceStats() {
    }

    public static long startTimer() {
        return ENABLED ? System.nanoTime() : 0L;
    }

    public static long elapsed(final long startNanos) {
        return ENABLED ? System.nanoTime() - startNanos : 0L;
    }

    public static void recordRewardProviderRefresh(final long nanos) {
        if (ENABLED) {
            rewardProviderRefreshCalls.incrementAndGet();
            rewardProviderRefreshNanos.addAndGet(nanos);
        }
    }

    public static void recordPostOrderRefresh(final long nanos) {
        if (ENABLED) {
            postOrderRefreshCalls.incrementAndGet();
            postOrderRefreshNanos.addAndGet(nanos);
        }
    }

    public static void recordPreOrderRefresh(final long nanos) {
        if (ENABLED) {
            preOrderRefreshCalls.incrementAndGet();
            preOrderRefreshNanos.addAndGet(nanos);
        }
    }

    public static void recordDiscretePreOrderEnsure(final boolean hit, final long nanos) {
        if (ENABLED) {
            if (hit) {
                discretePreOrderEnsureHits.incrementAndGet();
                discretePreOrderEnsureHitNanos.addAndGet(nanos);
            } else {
                discretePreOrderEnsureMisses.incrementAndGet();
                discretePreOrderEnsureMissNanos.addAndGet(nanos);
            }
        }
    }

    public static void recordDependentPrepare(final long nanos) {
        if (ENABLED) {
            dependentPrepareCalls.incrementAndGet();
            dependentPrepareNanos.addAndGet(nanos);
        }
    }

    public static void recordDiscreteLikelihoodTraversal(final int branchOperations,
                                                         final int nodeOperations,
                                                         final int nodeCount,
                                                         final int tipCount,
                                                         final boolean statsOnly,
                                                         final boolean substitutionDirty,
                                                         final boolean siteDirty,
                                                         final boolean rootFrequencyDirty) {
        if (ENABLED) {
            discreteLikelihoodTraversalCalls.incrementAndGet();
            if (statsOnly) {
                discreteLikelihoodTraversalStatsOnlyCalls.incrementAndGet();
            }
            discreteLikelihoodTraversalBranchOps.addAndGet(branchOperations);
            discreteLikelihoodTraversalNodeOps.addAndGet(nodeOperations);

            if (branchOperations == 0 && nodeOperations == 0) {
                discreteLikelihoodTraversalEmptyCalls.incrementAndGet();
            } else if (branchOperations >= Math.max(0, nodeCount - 1) &&
                    nodeOperations >= Math.max(0, nodeCount - tipCount)) {
                discreteLikelihoodTraversalFullTreeCalls.incrementAndGet();
            } else {
                discreteLikelihoodTraversalPartialTreeCalls.incrementAndGet();
            }

            if (substitutionDirty) {
                discreteLikelihoodTraversalSubstitutionDirtyCalls.incrementAndGet();
            }
            if (siteDirty) {
                discreteLikelihoodTraversalSiteDirtyCalls.incrementAndGet();
            }
            if (rootFrequencyDirty) {
                discreteLikelihoodTraversalRootFrequencyDirtyCalls.incrementAndGet();
            }
        }
    }

    public static void recordOperationCacheBegin() {
        if (ENABLED) {
            operationCacheBegins.incrementAndGet();
        }
    }

    public static void recordOperationCacheClear() {
        recordOperationCacheClear(OperationCacheClearReason.UNKNOWN);
    }

    public static void recordOperationCacheClear(final OperationCacheClearReason reason) {
        if (ENABLED) {
            operationCacheClears.incrementAndGet();
            switch (reason == null ? OperationCacheClearReason.UNKNOWN : reason) {
                case OPERATION_START:
                    operationCacheClearsOperationStart.incrementAndGet();
                    break;
                case SMOOTH_POSITION_UPDATE:
                    operationCacheClearsSmoothPositionUpdate.incrementAndGet();
                    break;
                case ACCEPTED_CATEGORY_CROSSING:
                    operationCacheClearsAcceptedCategoryCrossing.incrementAndGet();
                    break;
                case GIBBS_CATEGORY_CHANGE:
                    operationCacheClearsGibbsCategoryChange.incrementAndGet();
                    break;
                case FINAL_POSITION_REFRESH:
                    operationCacheClearsFinalPositionRefresh.incrementAndGet();
                    break;
                case DYNAMIC_EMBEDDING_REFRESH:
                    operationCacheClearsDynamicEmbeddingRefresh.incrementAndGet();
                    break;
                case UNKNOWN:
                default:
                    operationCacheClearsUnknown.incrementAndGet();
                    break;
            }
        }
    }

    public static void recordBranchWeightCacheHit() {
        if (ENABLED) {
            branchWeightCacheHits.incrementAndGet();
        }
    }

    public static void recordBranchWeightCacheMiss() {
        if (ENABLED) {
            branchWeightCacheMisses.incrementAndGet();
        }
    }

    public static void recordBranchWeightCacheBypass() {
        if (ENABLED) {
            branchWeightCacheBypasses.incrementAndGet();
        }
    }

    public static void recordLikelihoodMessageMiss() {
        if (ENABLED) {
            likelihoodMessageMisses.incrementAndGet();
        }
    }

    public static void recordBranchWeightComputation(final long nanos) {
        if (ENABLED) {
            branchWeightComputations.incrementAndGet();
            branchWeightComputationNanos.addAndGet(nanos);
        }
    }

    public static void recordBranchMessageLoad(final long nanos) {
        if (ENABLED) {
            branchMessageLoadCalls.incrementAndGet();
            branchMessageLoadNanos.addAndGet(nanos);
        }
    }

    public static void recordBranchMessageCacheFill(final long nanos) {
        if (ENABLED) {
            branchMessageCacheFillCalls.incrementAndGet();
            branchMessageCacheFillNanos.addAndGet(nanos);
        }
    }

    public static void recordBranchPrePartialLoad(final long nanos) {
        if (ENABLED) {
            branchPrePartialLoadCalls.incrementAndGet();
            branchPrePartialLoadNanos.addAndGet(nanos);
        }
    }

    public static void recordBranchPostPartialLoad(final long nanos) {
        if (ENABLED) {
            branchPostPartialLoadCalls.incrementAndGet();
            branchPostPartialLoadNanos.addAndGet(nanos);
        }
    }

    public static void recordBranchPreScaleLoad(final long nanos) {
        if (ENABLED) {
            branchPreScaleLoadCalls.incrementAndGet();
            branchPreScaleLoadNanos.addAndGet(nanos);
        }
    }

    public static void recordBranchPostScaleLoad(final long nanos) {
        if (ENABLED) {
            branchPostScaleLoadCalls.incrementAndGet();
            branchPostScaleLoadNanos.addAndGet(nanos);
        }
    }

    public static void recordBranchAtomicWeightLoop(final long nanos) {
        if (ENABLED) {
            branchAtomicWeightLoopCalls.incrementAndGet();
            branchAtomicWeightLoopNanos.addAndGet(nanos);
        }
    }

    public static void recordBranchContinuousMatrixAccess(final long nanos) {
        if (ENABLED) {
            branchContinuousMatrixAccessCalls.incrementAndGet();
            branchContinuousMatrixAccessNanos.addAndGet(nanos);
        }
    }

    public static void recordBranchContinuousWeight(final long nanos) {
        if (ENABLED) {
            branchContinuousWeightCalls.incrementAndGet();
            branchContinuousWeightNanos.addAndGet(nanos);
        }
    }

    public static void recordDependentEvidence(final long nanos) {
        if (ENABLED) {
            dependentEvidenceCalls.incrementAndGet();
            dependentEvidenceNanos.addAndGet(nanos);
        }
    }

    public static void recordBeaglePrepare(final long nanos) {
        if (ENABLED) {
            beaglePrepareCalls.incrementAndGet();
            beaglePrepareNanos.addAndGet(nanos);
        }
    }

    public static void addBeagleRefreshSpectralNanos(final long nanos) {
        if (ENABLED) {
            beagleRefreshSpectralNanos.addAndGet(nanos);
        }
    }

    public static void addBeagleFillPostNanos(final long nanos) {
        if (ENABLED) {
            beagleFillPostNanos.addAndGet(nanos);
        }
    }

    public static void addBeagleFillTopNanos(final long nanos) {
        if (ENABLED) {
            beagleFillTopNanos.addAndGet(nanos);
        }
    }

    public static void addBeagleFillPostTopNanos(final long nanos) {
        if (ENABLED) {
            beagleFillPostTopNanos.addAndGet(nanos);
        }
    }

    public static void addBeagleFillPreBottomNanos(final long nanos) {
        if (ENABLED) {
            beagleFillPreBottomNanos.addAndGet(nanos);
        }
    }

    public static void addBeagleFillPreOrderNanos(final long nanos) {
        if (ENABLED) {
            beagleFillPreOrderNanos.addAndGet(nanos);
        }
    }

    public static void recordBeagleLogEvidence(final long nanos, final boolean atomicCandidate) {
        if (ENABLED) {
            beagleLogEvidenceCalls.incrementAndGet();
            beagleLogEvidenceNanos.addAndGet(nanos);
            if (atomicCandidate) {
                beagleLogEvidenceAtomicCalls.incrementAndGet();
            } else {
                beagleLogEvidenceContinuousCalls.incrementAndGet();
            }
        }
    }

    public static void recordBeagleCachedMessageEvidence(final long nanos) {
        if (ENABLED) {
            beagleCachedMessageEvidenceCalls.incrementAndGet();
            beagleCachedMessageEvidenceNanos.addAndGet(nanos);
        }
    }

    public static void recordBeagleEdgeInnerProduct(final long nanos) {
        if (ENABLED) {
            beagleEdgeInnerProductCalls.incrementAndGet();
            beagleEdgeInnerProductNanos.addAndGet(nanos);
        }
    }

    public static void recordBeagleMessageCopy(final long nanos) {
        if (ENABLED) {
            beagleMessageCopyCalls.incrementAndGet();
            beagleMessageCopyNanos.addAndGet(nanos);
        }
    }

    public static void recordDiscontinuousStep(final boolean crossed, final boolean reflected) {
        if (ENABLED) {
            discontinuousSteps.incrementAndGet();
            if (crossed) {
                discontinuousAcceptedCrossings.incrementAndGet();
            }
            if (reflected) {
                discontinuousReflections.incrementAndGet();
            }
        }
    }

    public static void recordCacheClearAfterAcceptedCrossing() {
        if (ENABLED) {
            cacheClearsAfterAcceptedCrossing.incrementAndGet();
        }
    }

    public static void recordSkippedCacheClearAfterNoCrossing() {
        if (ENABLED) {
            skippedCacheClearsAfterNoCrossing.incrementAndGet();
        }
    }

    public static void recordCategoricalGibbsBranchUpdate(final boolean categoryChanged) {
        if (ENABLED) {
            categoricalGibbsBranchUpdates.incrementAndGet();
            if (categoryChanged) {
                categoricalGibbsCategoryChanges.incrementAndGet();
            } else {
                categoricalGibbsSameCategoryUpdates.incrementAndGet();
            }
        }
    }

    public static void recordCategoricalGibbsSkippedCacheClearAfterSameCategory() {
        if (ENABLED) {
            categoricalGibbsSkippedCacheClearsAfterSameCategory.incrementAndGet();
        }
    }

    private static void printReport() {
        final long hits = branchWeightCacheHits.get();
        final long misses = branchWeightCacheMisses.get();
        final long cachedLookups = hits + misses;

        System.err.println("Reward-mixture performance counters:");
        printCountTime("  refreshLikelihoodMessages", rewardProviderRefreshCalls, rewardProviderRefreshNanos);
        printCountTime("  primary postorder refresh", postOrderRefreshCalls, postOrderRefreshNanos);
        printCountTime("  primary preorder refresh", preOrderRefreshCalls, preOrderRefreshNanos);
        printCountTime("  discrete preorder ensure hits", discretePreOrderEnsureHits, discretePreOrderEnsureHitNanos);
        printCountTime("  discrete preorder ensure recomputes", discretePreOrderEnsureMisses, discretePreOrderEnsureMissNanos);
        printCountTime("  dependent provider prepare", dependentPrepareCalls, dependentPrepareNanos);
        System.err.println("  discrete likelihood traversals calls=" +
                discreteLikelihoodTraversalCalls.get() +
                " stats_only=" + discreteLikelihoodTraversalStatsOnlyCalls.get() +
                " full_tree=" + discreteLikelihoodTraversalFullTreeCalls.get() +
                " partial_tree=" + discreteLikelihoodTraversalPartialTreeCalls.get() +
                " empty=" + discreteLikelihoodTraversalEmptyCalls.get() +
                " branch_ops_total=" + discreteLikelihoodTraversalBranchOps.get() +
                " branch_ops_mean=" + mean(discreteLikelihoodTraversalBranchOps.get(),
                discreteLikelihoodTraversalCalls.get()) +
                " node_ops_total=" + discreteLikelihoodTraversalNodeOps.get() +
                " node_ops_mean=" + mean(discreteLikelihoodTraversalNodeOps.get(),
                discreteLikelihoodTraversalCalls.get()));
        System.err.println("  discrete likelihood dirty flags before traversal: substitution=" +
                discreteLikelihoodTraversalSubstitutionDirtyCalls.get() +
                " site=" + discreteLikelihoodTraversalSiteDirtyCalls.get() +
                " root_frequency=" + discreteLikelihoodTraversalRootFrequencyDirtyCalls.get());
        System.err.println("  operation cache begins=" + operationCacheBegins.get() +
                " clears=" + operationCacheClears.get() +
                " likelihood_message_misses=" + likelihoodMessageMisses.get());
        System.err.println("  operation cache clear reasons: operation_start=" +
                operationCacheClearsOperationStart.get() +
                " smooth_position_update=" + operationCacheClearsSmoothPositionUpdate.get() +
                " accepted_category_crossing=" + operationCacheClearsAcceptedCategoryCrossing.get() +
                " gibbs_category_change=" + operationCacheClearsGibbsCategoryChange.get() +
                " final_position_refresh=" + operationCacheClearsFinalPositionRefresh.get() +
                " dynamic_embedding_refresh=" + operationCacheClearsDynamicEmbeddingRefresh.get() +
                " unknown=" + operationCacheClearsUnknown.get());
        System.err.println("  branch weight cache lookups=" + cachedLookups +
                " hits=" + hits +
                " misses=" + misses +
                " bypasses=" + branchWeightCacheBypasses.get() +
                " hit_rate=" + rate(hits, cachedLookups));
        printCountTime("  branch weight computations", branchWeightComputations, branchWeightComputationNanos);
        printCountTime("  branch message loads", branchMessageLoadCalls, branchMessageLoadNanos);
        printCountTime("  branch message cache fills", branchMessageCacheFillCalls, branchMessageCacheFillNanos);
        printCountTime("  branch pre-partial loads", branchPrePartialLoadCalls, branchPrePartialLoadNanos);
        printCountTime("  branch post-partial loads", branchPostPartialLoadCalls, branchPostPartialLoadNanos);
        printCountTime("  branch pre-scale loads", branchPreScaleLoadCalls, branchPreScaleLoadNanos);
        printCountTime("  branch post-scale loads", branchPostScaleLoadCalls, branchPostScaleLoadNanos);
        printCountTime("  branch atomic weight loops", branchAtomicWeightLoopCalls, branchAtomicWeightLoopNanos);
        printCountTime("  branch continuous matrix access", branchContinuousMatrixAccessCalls, branchContinuousMatrixAccessNanos);
        printCountTime("  branch continuous weight algebra", branchContinuousWeightCalls, branchContinuousWeightNanos);
        printCountTime("  dependent evidence calls", dependentEvidenceCalls, dependentEvidenceNanos);
        printCountTime("  BEAGLE dependent prepare", beaglePrepareCalls, beaglePrepareNanos);
        printTimeOnly("  BEAGLE refresh spectral", beagleRefreshSpectralNanos);
        printTimeOnly("  BEAGLE fill postorder", beagleFillPostNanos);
        printTimeOnly("  BEAGLE fill preorder-top", beagleFillTopNanos);
        printTimeOnly("  BEAGLE fill post-top", beagleFillPostTopNanos);
        printTimeOnly("  BEAGLE fill pre-bottom", beagleFillPreBottomNanos);
        printTimeOnly("  BEAGLE native preorder", beagleFillPreOrderNanos);
        System.err.println("  BEAGLE logEvidence calls=" + beagleLogEvidenceCalls.get() +
                " atomic=" + beagleLogEvidenceAtomicCalls.get() +
                " continuous=" + beagleLogEvidenceContinuousCalls.get() +
                " total_ms=" + millis(beagleLogEvidenceNanos.get()));
        printCountTime("  BEAGLE cached-message evidence", beagleCachedMessageEvidenceCalls, beagleCachedMessageEvidenceNanos);
        printCountTime("  BEAGLE edge inner product", beagleEdgeInnerProductCalls, beagleEdgeInnerProductNanos);
        printCountTime("  BEAGLE message copies", beagleMessageCopyCalls, beagleMessageCopyNanos);
        System.err.println("  discontinuous steps=" + discontinuousSteps.get() +
                " crossed=" + discontinuousAcceptedCrossings.get() +
                " reflected=" + discontinuousReflections.get() +
                " cache_clears_after_crossing=" + cacheClearsAfterAcceptedCrossing.get() +
                " skipped_cache_clears_after_no_crossing=" + skippedCacheClearsAfterNoCrossing.get());
        System.err.println("  categorical Gibbs branch updates=" + categoricalGibbsBranchUpdates.get() +
                " category_changes=" + categoricalGibbsCategoryChanges.get() +
                " same_category=" + categoricalGibbsSameCategoryUpdates.get() +
                " skipped_cache_clears_after_same_category=" +
                categoricalGibbsSkippedCacheClearsAfterSameCategory.get());
    }

    private static void printCountTime(final String label,
                                       final AtomicLong count,
                                       final AtomicLong nanos) {
        final long c = count.get();
        final long n = nanos.get();
        System.err.println(label + " calls=" + c +
                " total_ms=" + millis(n) +
                " mean_us=" + meanMicros(n, c));
    }

    private static void printTimeOnly(final String label, final AtomicLong nanos) {
        System.err.println(label + " total_ms=" + millis(nanos.get()));
    }

    private static String millis(final long nanos) {
        return String.format("%.3f", nanos / 1.0e6);
    }

    private static String meanMicros(final long nanos, final long count) {
        if (count == 0L) {
            return "NA";
        }
        return String.format("%.3f", nanos / (count * 1.0e3));
    }

    private static String rate(final long numerator, final long denominator) {
        if (denominator == 0L) {
            return "NA";
        }
        return String.format("%.4f", ((double) numerator) / denominator);
    }

    private static String mean(final long total, final long count) {
        if (count == 0L) {
            return "NA";
        }
        return String.format("%.3f", ((double) total) / count);
    }
}

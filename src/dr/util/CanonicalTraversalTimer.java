package dr.util;

import java.util.Locale;

/**
 * Debug-only aggregate timing for canonical OU tree traversal internals.
 */
public final class CanonicalTraversalTimer {

    public static final String PROPERTY = "beast.debug.canonicalTraversalTiming";
    public static final String REPORT_EVERY_PROPERTY = "beast.debug.canonicalTraversalTimingReportEvery";

    private static final boolean ENABLED = Boolean.getBoolean(PROPERTY);
    private static final int REPORT_EVERY = Math.max(1, Integer.getInteger(REPORT_EVERY_PROPERTY, 50));

    private static long postorderTotalNanos;
    private static long postorderPreloadNanos;
    private static long postorderTransitionNanos;
    private static long postorderTipNanos;
    private static long postorderInternalPushNanos;
    private static long postorderCombineNanos;
    private static long postorderRootNanos;
    private static long preorderTotalNanos;
    private static long preorderRootInitNanos;
    private static long preorderOutsideTransitionNanos;
    private static long preorderOutsidePropagateNanos;
    private static long preorderSiblingProductNanos;
    private static long preorderStoreNanos;
    private static long postorderCount;
    private static long preorderCount;
    private static long postorderTips;
    private static long postorderInternalChildren;
    private static long postorderTransitionCalls;
    private static long preorderEdges;
    private static long preorderSiblingEdges;
    private static long preorderOutsideTransitionCalls;
    private static long postorderCacheRequests;
    private static long postorderCacheMisses;
    private static long preorderCacheRequests;
    private static long preorderCacheMisses;
    private static int lastNodeCount;
    private static int lastLeafCount;
    private static int lastInternalCount;
    private static int lastLevelCount;
    private static int lastMaxLevelWidth;
    private static double lastMeanLevelWidth;

    private CanonicalTraversalTimer() { }

    public static boolean isEnabled() {
        return ENABLED;
    }

    public static long start() {
        return ENABLED ? System.nanoTime() : 0L;
    }

    public static void recordTreeShape(final int nodeCount,
                                       final int leafCount,
                                       final int internalCount,
                                       final int levelCount,
                                       final int maxLevelWidth,
                                       final double meanLevelWidth) {
        if (ENABLED) {
            lastNodeCount = nodeCount;
            lastLeafCount = leafCount;
            lastInternalCount = internalCount;
            lastLevelCount = levelCount;
            lastMaxLevelWidth = maxLevelWidth;
            lastMeanLevelWidth = meanLevelWidth;
        }
    }

    public static void finishPostorder(final long startNanos,
                                       final long cacheRequests,
                                       final long cacheMisses) {
        if (ENABLED) {
            postorderTotalNanos += elapsedSince(startNanos);
            postorderCacheRequests += cacheRequests;
            postorderCacheMisses += cacheMisses;
            postorderCount++;
        }
    }

    public static void finishPreorder(final long startNanos,
                                      final long cacheRequests,
                                      final long cacheMisses) {
        if (!ENABLED) {
            return;
        }
        preorderTotalNanos += elapsedSince(startNanos);
        preorderCacheRequests += cacheRequests;
        preorderCacheMisses += cacheMisses;
        preorderCount++;
        if (preorderCount % REPORT_EVERY == 0) {
            reportAndReset();
        }
    }

    public static void recordPostorderTransition(final long startNanos) {
        if (ENABLED) {
            postorderTransitionNanos += elapsedSince(startNanos);
            postorderTransitionCalls++;
        }
    }

    public static void recordPostorderPreload(final long startNanos) {
        if (ENABLED) {
            postorderPreloadNanos += elapsedSince(startNanos);
        }
    }

    public static void recordPostorderTipMessage(final long startNanos) {
        if (ENABLED) {
            postorderTipNanos += elapsedSince(startNanos);
            postorderTips++;
        }
    }

    public static void recordPostorderInternalPush(final long startNanos) {
        if (ENABLED) {
            postorderInternalPushNanos += elapsedSince(startNanos);
            postorderInternalChildren++;
        }
    }

    public static void recordPostorderCombine(final long startNanos) {
        if (ENABLED) {
            postorderCombineNanos += elapsedSince(startNanos);
        }
    }

    public static void recordPostorderRoot(final long startNanos) {
        if (ENABLED) {
            postorderRootNanos += elapsedSince(startNanos);
        }
    }

    public static void recordPreorderRootInit(final long startNanos) {
        if (ENABLED) {
            preorderRootInitNanos += elapsedSince(startNanos);
        }
    }

    public static void recordPreorderOutsideTransition(final long startNanos) {
        if (ENABLED) {
            preorderOutsideTransitionNanos += elapsedSince(startNanos);
            preorderOutsideTransitionCalls++;
        }
    }

    public static void recordPreorderOutsidePropagate(final long startNanos) {
        if (ENABLED) {
            preorderOutsidePropagateNanos += elapsedSince(startNanos);
        }
    }

    public static void recordPreorderSiblingProduct(final long startNanos,
                                                    final int siblingEdges) {
        if (ENABLED) {
            preorderSiblingProductNanos += elapsedSince(startNanos);
            preorderSiblingEdges += siblingEdges;
        }
    }

    public static void recordPreorderStore(final long startNanos) {
        if (ENABLED) {
            preorderStoreNanos += elapsedSince(startNanos);
            preorderEdges++;
        }
    }

    private static long elapsedSince(final long startNanos) {
        return System.nanoTime() - startNanos;
    }

    private static void reportAndReset() {
        final String message = String.format(Locale.US,
                "[canonical-traversal-timing] post.count=%d pre.count=%d nodes=%d leaves=%d internal=%d "
                        + "levels=%d maxWidth=%d meanWidth=%.1f "
                        + "post.total=%.3fms post.preload=%.3fms post.transition=%.3fms post.tip=%.3fms "
                        + "post.internalPush=%.3fms post.combine=%.3fms post.root=%.3fms "
                        + "post.edges=transition/tip/internal=%d/%d/%d "
                        + "post.cache=requests/misses=%d/%d "
                        + "pre.total=%.3fms pre.rootInit=%.3fms pre.outsideTransition=%.3fms "
                        + "pre.outsidePropagate=%.3fms pre.siblingProduct=%.3fms pre.store=%.3fms "
                        + "pre.edges=outside/transition/sibling=%d/%d/%d "
                        + "pre.cache=requests/misses=%d/%d",
                postorderCount,
                preorderCount,
                lastNodeCount,
                lastLeafCount,
                lastInternalCount,
                lastLevelCount,
                lastMaxLevelWidth,
                lastMeanLevelWidth,
                toMillis(postorderTotalNanos),
                toMillis(postorderPreloadNanos),
                toMillis(postorderTransitionNanos),
                toMillis(postorderTipNanos),
                toMillis(postorderInternalPushNanos),
                toMillis(postorderCombineNanos),
                toMillis(postorderRootNanos),
                postorderTransitionCalls,
                postorderTips,
                postorderInternalChildren,
                postorderCacheRequests,
                postorderCacheMisses,
                toMillis(preorderTotalNanos),
                toMillis(preorderRootInitNanos),
                toMillis(preorderOutsideTransitionNanos),
                toMillis(preorderOutsidePropagateNanos),
                toMillis(preorderSiblingProductNanos),
                toMillis(preorderStoreNanos),
                preorderEdges,
                preorderOutsideTransitionCalls,
                preorderSiblingEdges,
                preorderCacheRequests,
                preorderCacheMisses);
        System.err.println(message);
        reset();
    }

    private static double toMillis(final long nanos) {
        return nanos / 1.0e6;
    }

    private static void reset() {
        postorderTotalNanos = 0L;
        postorderPreloadNanos = 0L;
        postorderTransitionNanos = 0L;
        postorderTipNanos = 0L;
        postorderInternalPushNanos = 0L;
        postorderCombineNanos = 0L;
        postorderRootNanos = 0L;
        preorderTotalNanos = 0L;
        preorderRootInitNanos = 0L;
        preorderOutsideTransitionNanos = 0L;
        preorderOutsidePropagateNanos = 0L;
        preorderSiblingProductNanos = 0L;
        preorderStoreNanos = 0L;
        postorderCount = 0L;
        preorderCount = 0L;
        postorderTips = 0L;
        postorderInternalChildren = 0L;
        postorderTransitionCalls = 0L;
        preorderEdges = 0L;
        preorderSiblingEdges = 0L;
        preorderOutsideTransitionCalls = 0L;
        postorderCacheRequests = 0L;
        postorderCacheMisses = 0L;
        preorderCacheRequests = 0L;
        preorderCacheMisses = 0L;
    }
}

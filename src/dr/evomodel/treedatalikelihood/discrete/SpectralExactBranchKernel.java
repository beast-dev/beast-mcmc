package dr.evomodel.treedatalikelihood.discrete;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.substmodel.EigenDecomposition;
import dr.evomodel.substmodel.LogAdditiveCtmcRateProvider;
import dr.evomodel.substmodel.LogRateSubstitutionModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.treedatalikelihood.DiscreteDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.discrete.discretetreedataLikelihood.ComplexBlockKernelUtils;
import dr.evomodel.treedatalikelihood.discrete.discretetreedataLikelihood.RealKernelUtils;
import dr.evomodel.treedatalikelihood.discrete.discretetreedataLikelihood.representations.RealBranchExponentialsProvider;
import dr.inference.model.Parameter;
import dr.util.Transform;

import java.util.Arrays;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.*;

final class SpectralExactBranchKernel {

    /** System property to control parallelism; defaults to 1 (single-threaded). */
    private static final String THREADS_PROPERTY = "beast.gradient.threads";

    private final Tree tree;
    private final DiscreteDataLikelihoodDelegate likelihoodDelegate;
    private final BranchModel branchModel;
    private final int stateCount;
    private final int substitutionModelCount;

    // Per-model accumulation in the rotated (eigen) basis (used only in single-thread path)
    private final double[][] eigenBasisAccumByModel;

    // Scratch for the final rotation step
    private final double[] midBuffer;
    private double[] branchDifferentials;
    private double[] branchLogRateScores;

    // Parallelism
    private final int nThreads;
    private final ExecutorService pool;
    private final ThreadWorkspace[] workspaces;

    // When true, skip repeated RealKernelUtils.isAllReal() detection per gradient call.
    private final boolean forceAllReal;
    // Pre-allocated isAllReal array — reused across gradient calls to avoid hot-path allocation.
    private final boolean[] isAllRealCache;
    private final EigenDecomposition[] eigenDecomps;
    private final int[][] directedRateIndex;
    private final double[] generatorBuffer;

    // -------------------------------------------------------------------------
    // Per-thread workspace
    // -------------------------------------------------------------------------

    private final class ThreadWorkspace {
        final double[][] eigenBasisAccum;
        final double[] tmpPreTop;
        final double[] tmpPreBottom;
        final double[] tmpPostTop;
        final double[] tmpPostBottom;
        final double[] rotatedPre;
        final double[] rotatedPost;
        final double[] branchEigenBasisAccum;
        final double[] branchStandardAccum;
        final double[] eigenRowProjection;
        final double[] standardRowProjection;
        final double[] midBuffer;
        final ComplexBlockKernelUtils.ComplexKernelPlan[] planByModel;
        final RealKernelUtils.RealKernelPlan[] realPlanByModel;
        final RealBranchExponentialsProvider.BorrowedSlice realBranchExponentials;

        ThreadWorkspace() {
            final int K2 = stateCount * stateCount;
            eigenBasisAccum  = new double[substitutionModelCount][K2];
            tmpPreTop        = new double[stateCount];
            tmpPreBottom     = new double[stateCount];
            tmpPostTop       = new double[stateCount];
            tmpPostBottom    = new double[stateCount];
            rotatedPre       = new double[stateCount];
            rotatedPost      = new double[stateCount];
            branchEigenBasisAccum = new double[K2];
            branchStandardAccum = new double[K2];
            eigenRowProjection = new double[stateCount];
            standardRowProjection = new double[stateCount];
            midBuffer       = new double[K2];
            planByModel      = new ComplexBlockKernelUtils.ComplexKernelPlan[substitutionModelCount];
            realPlanByModel  = new RealKernelUtils.RealKernelPlan[substitutionModelCount];
            realBranchExponentials = new RealBranchExponentialsProvider.BorrowedSlice();
            for (int i = 0; i < substitutionModelCount; i++) {
                planByModel[i]     = new ComplexBlockKernelUtils.ComplexKernelPlan(stateCount);
                realPlanByModel[i] = new RealKernelUtils.RealKernelPlan(stateCount);
            }
        }
    }

    // -------------------------------------------------------------------------

    SpectralExactBranchKernel(Tree tree,
                              DiscreteDataLikelihoodDelegate likelihoodDelegate,
                              int stateCount,
                              boolean forceAllReal) {
        this.tree = tree;
        this.likelihoodDelegate = likelihoodDelegate;
        this.branchModel = likelihoodDelegate.getBranchModel();
        this.stateCount = stateCount;
        this.substitutionModelCount = branchModel.getSubstitutionModels().size();

        final int K2 = stateCount * stateCount;
        this.eigenBasisAccumByModel = new double[substitutionModelCount][K2];
        this.midBuffer = new double[K2];

        this.forceAllReal   = forceAllReal;
        this.isAllRealCache = new boolean[substitutionModelCount];
        this.eigenDecomps = new EigenDecomposition[substitutionModelCount];
        this.directedRateIndex = makeDirectedRateIndex(stateCount);
        this.generatorBuffer = new double[stateCount * stateCount];
        if (forceAllReal) {
            Arrays.fill(isAllRealCache, true);
        }

        this.nThreads = Integer.getInteger(THREADS_PROPERTY, 1);
        if (nThreads > 1) {
            this.pool = Executors.newFixedThreadPool(nThreads, r -> {
                Thread t = new Thread(r, "gradient-worker");
                t.setDaemon(true);
                return t;
            });
        } else {
            this.pool = null;
        }
        this.workspaces = new ThreadWorkspace[Math.max(nThreads, 1)];
        for (int t = 0; t < workspaces.length; t++) {
            workspaces[t] = new ThreadWorkspace();
        }
    }

    int getGradientLength() {
        return stateCount * stateCount * substitutionModelCount;
    }

    int getBranchDifferentialLength() {
        return tree.getNodeCount() * substitutionModelCount * stateCount * stateCount;
    }

    int getBranchLogRateScoreLength() {
        return tree.getNodeCount() * substitutionModelCount * getDirectedRateDimension();
    }

    void computeGradient(double[] first) {
        final int K2 = stateCount * stateCount;
        Arrays.fill(first, 0.0);

        likelihoodDelegate.ensurePreOrderComputed();

        final List<SubstitutionModel> models = branchModel.getSubstitutionModels();
        final double[] patternWeights  = likelihoodDelegate.getPatternWeights();
        final double[] categoryWeights = likelihoodDelegate.getCategoryWeights();
        final double[] categoryRates   = likelihoodDelegate.getSiteRateModel().getCategoryRates();
        final boolean useInternalRotatedMessages = likelihoodDelegate.isSpectralRepresentation();

        // Cache all eigen decompositions once — avoids 598× synchronized calls per gradient eval.
        for (int m = 0; m < substitutionModelCount; m++) {
            eigenDecomps[m] = models.get(m).getEigenDecomposition();
        }

        // Detect all-real eigensystems per model (reversible CTMCs).
        // When forceAllReal=true the cache is already all-true; skip the O(K) detection calls.
        if (!forceAllReal) {
            for (int m = 0; m < substitutionModelCount; m++) {
                isAllRealCache[m] = eigenDecomps[m] != null &&
                        RealKernelUtils.isAllReal(eigenDecomps[m], stateCount);
            }
        }

        // Fill structural plan (time-independent) into every workspace.
        // Real models use RealKernelUtils (packed K² table); complex models use ComplexBlockKernelUtils.
        for (int m = 0; m < substitutionModelCount; m++) {
            if (eigenDecomps[m] != null) {
                for (ThreadWorkspace ws : workspaces) {
                    if (isAllRealCache[m]) {
                        RealKernelUtils.fillStructure(ws.realPlanByModel[m], eigenDecomps[m], stateCount);
                    } else {
                        ComplexBlockKernelUtils.fillStructure(ws.planByModel[m], eigenDecomps[m]);
                    }
                }
            }
        }

        // Reset per-workspace accumulation buffers
        for (ThreadWorkspace ws : workspaces) {
            for (double[] buf : ws.eigenBasisAccum) {
                Arrays.fill(buf, 0.0);
            }
        }

        if (nThreads <= 1) {
            processBranchRange(0, tree.getNodeCount(), workspaces[0],
                    eigenDecomps, isAllRealCache, patternWeights, categoryWeights, categoryRates, useInternalRotatedMessages);
        } else {
            final int nodeCount = tree.getNodeCount();
            final int chunkSize = (nodeCount + nThreads - 1) / nThreads;
            final List<Callable<Void>> tasks = new ArrayList<>(nThreads);
            for (int t = 0; t < nThreads; t++) {
                final int start = t * chunkSize;
                final int end   = Math.min(start + chunkSize, nodeCount);
                if (start >= nodeCount) break;
                final ThreadWorkspace ws = workspaces[t];
                tasks.add(() -> {
                    processBranchRange(start, end, ws,
                            eigenDecomps, isAllRealCache, patternWeights, categoryWeights, categoryRates, useInternalRotatedMessages);
                    return null;
                });
            }
            try {
                final List<Future<Void>> futures = pool.invokeAll(tasks);
                for (Future<Void> f : futures) f.get();
            } catch (InterruptedException | ExecutionException e) {
                throw new RuntimeException("Parallel gradient computation failed", e);
            }
        }

        // Reduce: sum per-workspace accumulations into eigenBasisAccumByModel
        for (int m = 0; m < substitutionModelCount; m++) {
            final double[] dest = eigenBasisAccumByModel[m];
            Arrays.fill(dest, 0.0);
            for (ThreadWorkspace ws : workspaces) {
                final double[] src = ws.eigenBasisAccum[m];
                for (int k = 0; k < K2; k++) {
                    dest[k] += src[k];
                }
            }
        }

        // Rotate each model's accumulation matrix and write into first
        // gradient = R^{-T} * eigenBasisAccum * R^T
        for (int m = 0; m < substitutionModelCount; m++) {
            if (eigenDecomps[m] == null) continue;
            final double[] evec = eigenDecomps[m].getEigenVectors();
            final double[] ievc = eigenDecomps[m].getInverseEigenVectors();
            rotateIntoOutput(eigenBasisAccumByModel[m], evec, ievc, first, m * K2);
        }
    }

    double[] getBranchDifferentials() {
        final int requiredLength = getBranchDifferentialLength();
        if (branchDifferentials == null || branchDifferentials.length != requiredLength) {
            branchDifferentials = new double[requiredLength];
        }

        computeBranchDifferentials(branchDifferentials);
        return branchDifferentials;
    }

    double[] getBranchLogRateScores() {
        final int requiredLength = getBranchLogRateScoreLength();
        if (branchLogRateScores == null || branchLogRateScores.length != requiredLength) {
            branchLogRateScores = new double[requiredLength];
        }

        computeBranchLogRateScores(branchLogRateScores);
        return branchLogRateScores;
    }

    private void computeBranchLogRateScores(double[] destination) {
        Arrays.fill(destination, 0.0);

        likelihoodDelegate.ensurePreOrderComputed();

        final List<SubstitutionModel> models = branchModel.getSubstitutionModels();
        final double[] patternWeights  = likelihoodDelegate.getPatternWeights();
        final double[] categoryWeights = likelihoodDelegate.getCategoryWeights();
        final double[] categoryRates   = likelihoodDelegate.getSiteRateModel().getCategoryRates();
        final boolean useInternalRotatedMessages = likelihoodDelegate.isSpectralRepresentation();

        for (int m = 0; m < substitutionModelCount; m++) {
            eigenDecomps[m] = models.get(m).getEigenDecomposition();
        }

        if (!forceAllReal) {
            for (int m = 0; m < substitutionModelCount; m++) {
                isAllRealCache[m] = eigenDecomps[m] != null &&
                        RealKernelUtils.isAllReal(eigenDecomps[m], stateCount);
            }
        }

        for (int m = 0; m < substitutionModelCount; m++) {
            if (eigenDecomps[m] != null) {
                for (ThreadWorkspace ws : workspaces) {
                    if (isAllRealCache[m]) {
                        RealKernelUtils.fillStructure(ws.realPlanByModel[m], eigenDecomps[m], stateCount);
                    } else {
                        ComplexBlockKernelUtils.fillStructure(ws.planByModel[m], eigenDecomps[m]);
                    }
                }
            }
        }

        if (nThreads <= 1) {
            processBranchLogRateScoreRange(0, tree.getNodeCount(), workspaces[0],
                    models, eigenDecomps, isAllRealCache, patternWeights, categoryWeights, categoryRates,
                    useInternalRotatedMessages, destination);
        } else {
            final int nodeCount = tree.getNodeCount();
            final int chunkSize = (nodeCount + nThreads - 1) / nThreads;
            final List<Callable<Void>> tasks = new ArrayList<>(nThreads);
            for (int t = 0; t < nThreads; t++) {
                final int start = t * chunkSize;
                final int end = Math.min(start + chunkSize, nodeCount);
                if (start >= nodeCount) break;
                final ThreadWorkspace ws = workspaces[t];
                tasks.add(() -> {
                    processBranchLogRateScoreRange(start, end, ws,
                            models, eigenDecomps, isAllRealCache, patternWeights, categoryWeights, categoryRates,
                            useInternalRotatedMessages, destination);
                    return null;
                });
            }
            try {
                final List<Future<Void>> futures = pool.invokeAll(tasks);
                for (Future<Void> f : futures) f.get();
            } catch (InterruptedException | ExecutionException e) {
                throw new RuntimeException("Parallel branch log-rate score computation failed", e);
            }
        }

        finishRawLogRateScores(destination, models);
    }

    private void computeBranchDifferentials(double[] destination) {
        Arrays.fill(destination, 0.0);

        likelihoodDelegate.ensurePreOrderComputed();

        final List<SubstitutionModel> models = branchModel.getSubstitutionModels();
        final double[] patternWeights  = likelihoodDelegate.getPatternWeights();
        final double[] categoryWeights = likelihoodDelegate.getCategoryWeights();
        final double[] categoryRates   = likelihoodDelegate.getSiteRateModel().getCategoryRates();
        final boolean useInternalRotatedMessages = likelihoodDelegate.isSpectralRepresentation();

        for (int m = 0; m < substitutionModelCount; m++) {
            eigenDecomps[m] = models.get(m).getEigenDecomposition();
        }

        if (!forceAllReal) {
            for (int m = 0; m < substitutionModelCount; m++) {
                isAllRealCache[m] = eigenDecomps[m] != null &&
                        RealKernelUtils.isAllReal(eigenDecomps[m], stateCount);
            }
        }

        for (int m = 0; m < substitutionModelCount; m++) {
            if (eigenDecomps[m] != null) {
                for (ThreadWorkspace ws : workspaces) {
                    if (isAllRealCache[m]) {
                        RealKernelUtils.fillStructure(ws.realPlanByModel[m], eigenDecomps[m], stateCount);
                    } else {
                        ComplexBlockKernelUtils.fillStructure(ws.planByModel[m], eigenDecomps[m]);
                    }
                }
            }
        }

        if (nThreads <= 1) {
            processBranchDifferentialRange(0, tree.getNodeCount(), workspaces[0],
                    eigenDecomps, isAllRealCache, patternWeights, categoryWeights, categoryRates,
                    useInternalRotatedMessages, destination);
        } else {
            final int nodeCount = tree.getNodeCount();
            final int chunkSize = (nodeCount + nThreads - 1) / nThreads;
            final List<Callable<Void>> tasks = new ArrayList<>(nThreads);
            for (int t = 0; t < nThreads; t++) {
                final int start = t * chunkSize;
                final int end   = Math.min(start + chunkSize, nodeCount);
                if (start >= nodeCount) break;
                final ThreadWorkspace ws = workspaces[t];
                tasks.add(() -> {
                    processBranchDifferentialRange(start, end, ws,
                            eigenDecomps, isAllRealCache, patternWeights, categoryWeights, categoryRates,
                            useInternalRotatedMessages, destination);
                    return null;
                });
            }
            try {
                final List<Future<Void>> futures = pool.invokeAll(tasks);
                for (Future<Void> f : futures) f.get();
            } catch (InterruptedException | ExecutionException e) {
                throw new RuntimeException("Parallel branch-differential computation failed", e);
            }
        }
    }

    private void processBranchRange(int start, int end,
                                    ThreadWorkspace ws,
                                    EigenDecomposition[] eigenDecomps,
                                    boolean[] isAllReal,
                                    double[] patternWeights,
                                    double[] categoryWeights,
                                    double[] categoryRates,
                                    boolean useInternalRotatedMessages) {
        for (int childNumber = start; childNumber < end; childNumber++) {
            final NodeRef child = tree.getNode(childNumber);
            if (tree.isRoot(child)) {
                continue;
            }

            final double baseBranchLength = likelihoodDelegate.getEffectiveBranchLength(childNumber);

            final BranchModel.Mapping mapping = branchModel.getBranchModelMapping(child);
            if (mapping == BranchModel.DEFAULT) {
                processBranchComponent(childNumber, 0, baseBranchLength,
                        ws, eigenDecomps, isAllReal, patternWeights, categoryWeights, categoryRates,
                        useInternalRotatedMessages);
                continue;
            }

            final int[]    order   = mapping.getOrder();
            final double[] weights = mapping.getWeights();

            for (int mixtureIndex = 0; mixtureIndex < order.length; mixtureIndex++) {
                final int modelNumber    = order[mixtureIndex];
                final double modelWeight = relativeWeight(mixtureIndex, weights);
                final double branchLengthForModel = baseBranchLength * modelWeight;

                processBranchComponent(childNumber, modelNumber, branchLengthForModel,
                        ws, eigenDecomps, isAllReal, patternWeights, categoryWeights, categoryRates,
                        useInternalRotatedMessages);
            }
        }
    }

    private void processBranchDifferentialRange(int start, int end,
                                                ThreadWorkspace ws,
                                                EigenDecomposition[] eigenDecomps,
                                                boolean[] isAllReal,
                                                double[] patternWeights,
                                                double[] categoryWeights,
                                                double[] categoryRates,
                                                boolean useInternalRotatedMessages,
                                                double[] destination) {
        for (int childNumber = start; childNumber < end; childNumber++) {
            final NodeRef child = tree.getNode(childNumber);
            if (tree.isRoot(child)) {
                continue;
            }

            final double baseBranchLength = likelihoodDelegate.getEffectiveBranchLength(childNumber);

            final BranchModel.Mapping mapping = branchModel.getBranchModelMapping(child);
            if (mapping == BranchModel.DEFAULT) {
                processBranchDifferentialComponent(childNumber, 0, baseBranchLength,
                        ws, eigenDecomps, isAllReal, patternWeights, categoryWeights, categoryRates,
                        useInternalRotatedMessages, destination);
                continue;
            }

            final int[]    order   = mapping.getOrder();
            final double[] weights = mapping.getWeights();

            for (int mixtureIndex = 0; mixtureIndex < order.length; mixtureIndex++) {
                final int modelNumber    = order[mixtureIndex];
                final double modelWeight = relativeWeight(mixtureIndex, weights);
                final double branchLengthForModel = baseBranchLength * modelWeight;

                processBranchDifferentialComponent(childNumber, modelNumber, branchLengthForModel,
                        ws, eigenDecomps, isAllReal, patternWeights, categoryWeights, categoryRates,
                        useInternalRotatedMessages, destination);
            }
        }
    }

    private void processBranchLogRateScoreRange(int start, int end,
                                                ThreadWorkspace ws,
                                                List<SubstitutionModel> models,
                                                EigenDecomposition[] eigenDecomps,
                                                boolean[] isAllReal,
                                                double[] patternWeights,
                                                double[] categoryWeights,
                                                double[] categoryRates,
                                                boolean useInternalRotatedMessages,
                                                double[] destination) {
        for (int childNumber = start; childNumber < end; childNumber++) {
            final NodeRef child = tree.getNode(childNumber);
            if (tree.isRoot(child)) {
                continue;
            }

            final double baseBranchLength = likelihoodDelegate.getEffectiveBranchLength(childNumber);

            final BranchModel.Mapping mapping = branchModel.getBranchModelMapping(child);
            if (mapping == BranchModel.DEFAULT) {
                processBranchLogRateScoreComponent(childNumber, 0, baseBranchLength,
                        ws, models, eigenDecomps, isAllReal, patternWeights, categoryWeights, categoryRates,
                        useInternalRotatedMessages, destination);
                continue;
            }

            final int[] order = mapping.getOrder();
            final double[] weights = mapping.getWeights();

            for (int mixtureIndex = 0; mixtureIndex < order.length; mixtureIndex++) {
                final int modelNumber = order[mixtureIndex];
                final double modelWeight = relativeWeight(mixtureIndex, weights);
                final double branchLengthForModel = baseBranchLength * modelWeight;

                processBranchLogRateScoreComponent(childNumber, modelNumber, branchLengthForModel,
                        ws, models, eigenDecomps, isAllReal, patternWeights, categoryWeights, categoryRates,
                        useInternalRotatedMessages, destination);
            }
        }
    }

    private void processBranchComponent(int childNumber,
                                        int modelNumber,
                                        double branchLengthForModel,
                                        ThreadWorkspace ws,
                                        EigenDecomposition[] eigenDecomps,
                                        boolean[] isAllReal,
                                        double[] patternWeights,
                                        double[] categoryWeights,
                                        double[] categoryRates,
                                        boolean useInternalRotatedMessages) {
        final EigenDecomposition eigenDecomp = eigenDecomps[modelNumber];
        if (eigenDecomp == null) {
            return;
        }

        final double[] evec = eigenDecomp.getEigenVectors();
        final double[] ievc = eigenDecomp.getInverseEigenVectors();

        if (isAllReal[modelNumber]) {
            accumulateBranchContributionReal(
                    childNumber,
                    branchLengthForModel,
                    categoryRates,
                    categoryWeights,
                    patternWeights,
                    eigenDecomp,
                    ws.realPlanByModel[modelNumber],
                    evec,
                    ievc,
                    ws.eigenBasisAccum[modelNumber],
                    ws,
                    useInternalRotatedMessages
            );
        } else {
            accumulateBranchContribution(
                    childNumber,
                    branchLengthForModel,
                    categoryRates,
                    categoryWeights,
                    patternWeights,
                    eigenDecomp,
                    ws.planByModel[modelNumber],
                    evec,
                    ievc,
                    ws.eigenBasisAccum[modelNumber],
                    ws,
                    useInternalRotatedMessages
            );
        }
    }

    private void processBranchDifferentialComponent(int childNumber,
                                                    int modelNumber,
                                                    double branchLengthForModel,
                                                    ThreadWorkspace ws,
                                                    EigenDecomposition[] eigenDecomps,
                                                    boolean[] isAllReal,
                                                    double[] patternWeights,
                                                    double[] categoryWeights,
                                                    double[] categoryRates,
                                                    boolean useInternalRotatedMessages,
                                                    double[] destination) {
        final EigenDecomposition eigenDecomp = eigenDecomps[modelNumber];
        if (eigenDecomp == null) {
            return;
        }

        final double[] evec = eigenDecomp.getEigenVectors();
        final double[] ievc = eigenDecomp.getInverseEigenVectors();
        Arrays.fill(ws.branchEigenBasisAccum, 0.0);

        if (isAllReal[modelNumber]) {
            accumulateBranchContributionReal(
                    childNumber,
                    branchLengthForModel,
                    categoryRates,
                    categoryWeights,
                    patternWeights,
                    eigenDecomp,
                    ws.realPlanByModel[modelNumber],
                    evec,
                    ievc,
                    ws.branchEigenBasisAccum,
                    ws,
                    useInternalRotatedMessages
            );
        } else {
            accumulateBranchContribution(
                    childNumber,
                    branchLengthForModel,
                    categoryRates,
                    categoryWeights,
                    patternWeights,
                    eigenDecomp,
                    ws.planByModel[modelNumber],
                    evec,
                    ievc,
                    ws.branchEigenBasisAccum,
                    ws,
                    useInternalRotatedMessages
            );
        }

        final int offset = (childNumber * substitutionModelCount + modelNumber) * stateCount * stateCount;
        Arrays.fill(ws.branchStandardAccum, 0.0);
        rotateIntoOutput(ws.branchEigenBasisAccum, evec, ievc, ws.branchStandardAccum, 0, ws.midBuffer);
        for (int i = 0; i < ws.branchStandardAccum.length; ++i) {
            destination[offset + i] += ws.branchStandardAccum[i];
        }
    }

    private void processBranchLogRateScoreComponent(int childNumber,
                                                    int modelNumber,
                                                    double branchLengthForModel,
                                                    ThreadWorkspace ws,
                                                    List<SubstitutionModel> models,
                                                    EigenDecomposition[] eigenDecomps,
                                                    boolean[] isAllReal,
                                                    double[] patternWeights,
                                                    double[] categoryWeights,
                                                    double[] categoryRates,
                                                    boolean useInternalRotatedMessages,
                                                    double[] destination) {
        if (!isDirectLogRateModel(models.get(modelNumber))) {
            return;
        }

        final EigenDecomposition eigenDecomp = eigenDecomps[modelNumber];
        if (eigenDecomp == null) {
            return;
        }

        final double[] evec = eigenDecomp.getEigenVectors();
        final double[] ievc = eigenDecomp.getInverseEigenVectors();
        Arrays.fill(ws.branchEigenBasisAccum, 0.0);

        if (isAllReal[modelNumber]) {
            accumulateBranchContributionReal(
                    childNumber,
                    branchLengthForModel,
                    categoryRates,
                    categoryWeights,
                    patternWeights,
                    eigenDecomp,
                    ws.realPlanByModel[modelNumber],
                    evec,
                    ievc,
                    ws.branchEigenBasisAccum,
                    ws,
                    useInternalRotatedMessages
            );
        } else {
            accumulateBranchContribution(
                    childNumber,
                    branchLengthForModel,
                    categoryRates,
                    categoryWeights,
                    patternWeights,
                    eigenDecomp,
                    ws.planByModel[modelNumber],
                    evec,
                    ievc,
                    ws.branchEigenBasisAccum,
                    ws,
                    useInternalRotatedMessages
            );
        }

        final int offset = (childNumber * substitutionModelCount + modelNumber) * getDirectedRateDimension();
        addRawLogRateScoresFromEigenBasis(ws.branchEigenBasisAccum, evec, ievc,
                destination, offset, ws.eigenRowProjection, ws.standardRowProjection);
    }

    /**
     * For a single branch (one mixture component), accumulate the Fréchet integral
     * contributions across all categories and patterns into eigenBasisAccum.
     */
    private void accumulateBranchContribution(int childNumber,
                                               double branchLength,
                                               double[] categoryRates,
                                               double[] categoryWeights,
                                               double[] patternWeights,
                                               EigenDecomposition eigenDecomp,
                                               ComplexBlockKernelUtils.ComplexKernelPlan plan,
                                               double[] evec,
                                               double[] ievc,
                                               double[] eigenBasisAccum,
                                               ThreadWorkspace ws,
                                               boolean useInternalRotatedMessages) {
        final int categoryCount = likelihoodDelegate.getCategoryCount();
        final int patternCount  = likelihoodDelegate.getPatternCount();

        for (int c = 0; c < categoryCount; c++) {
            final double wc = categoryWeights[c];
            final double tc = branchLength * categoryRates[c];

            if (patternCount == 1) {
                // Fast path for single-pattern data (e.g. phylogeography): fuse coefficient
                // computation with the outer-product accumulation in one pass, eliminating
                // the K²×16 intermediate coefficient store/load (~100 MB/gradient at K=26).
                final double wp = patternWeights[0];
                if (wp == 0.0) continue;

                double denom = 0.0;
                if (useInternalRotatedMessages) {
                    likelihoodDelegate.getInternalPreOrderBranchTopInto(childNumber, c, 0, ws.rotatedPre);
                    likelihoodDelegate.getPostOrderBranchTopInto(childNumber, c, 0, ws.tmpPostTop);
                    likelihoodDelegate.getInternalPostOrderBranchBottomInto(childNumber, c, 0, ws.rotatedPost);
                    for (int s = 0; s < stateCount; s++) {
                        denom += ws.rotatedPre[s] * ws.tmpPostTop[s];
                    }
                } else {
                    likelihoodDelegate.getPreOrderBranchTopInto(childNumber, c, 0, ws.tmpPreTop);
                    likelihoodDelegate.getPreOrderBranchBottomInto(childNumber, c, 0, ws.tmpPreBottom);
                    likelihoodDelegate.getPostOrderBranchBottomInto(childNumber, c, 0, ws.tmpPostBottom);
                    for (int s = 0; s < stateCount; s++) {
                        denom += ws.tmpPreBottom[s] * ws.tmpPostBottom[s];
                    }
                }
                if (denom <= 0.0 || Double.isNaN(denom) || Double.isInfinite(denom)) continue;

                if (!useInternalRotatedMessages) {
                    multiplyTransposeMatrixVector(evec, ws.tmpPreTop, ws.rotatedPre);
                    multiplyMatrixVector(ievc, ws.tmpPostBottom, ws.rotatedPost);
                }

                ComplexBlockKernelUtils.fillAndApplyToOuterProduct(
                        plan, eigenDecomp, tc, ws.rotatedPre, ws.rotatedPost,
                        (wp * wc) / denom, eigenBasisAccum);
            } else {
                // Multi-pattern path: fill coefficients once, apply across all patterns.
                ComplexBlockKernelUtils.fillTimeDependentCoefficients(plan, eigenDecomp, tc);

                for (int p = 0; p < patternCount; p++) {
                    final double wp = patternWeights[p];
                    if (wp == 0.0) continue;

                    double denom = 0.0;
                    if (useInternalRotatedMessages) {
                        likelihoodDelegate.getInternalPreOrderBranchTopInto(childNumber, c, p, ws.rotatedPre);
                        likelihoodDelegate.getPostOrderBranchTopInto(childNumber, c, p, ws.tmpPostTop);
                        likelihoodDelegate.getInternalPostOrderBranchBottomInto(childNumber, c, p, ws.rotatedPost);
                        for (int s = 0; s < stateCount; s++) {
                            denom += ws.rotatedPre[s] * ws.tmpPostTop[s];
                        }
                    } else {
                        likelihoodDelegate.getPreOrderBranchTopInto(childNumber, c, p, ws.tmpPreTop);
                        likelihoodDelegate.getPreOrderBranchBottomInto(childNumber, c, p, ws.tmpPreBottom);
                        likelihoodDelegate.getPostOrderBranchBottomInto(childNumber, c, p, ws.tmpPostBottom);
                        for (int s = 0; s < stateCount; s++) {
                            denom += ws.tmpPreBottom[s] * ws.tmpPostBottom[s];
                        }
                    }
                    if (denom <= 0.0 || Double.isNaN(denom) || Double.isInfinite(denom)) continue;

                    if (!useInternalRotatedMessages) {
                        multiplyTransposeMatrixVector(evec, ws.tmpPreTop, ws.rotatedPre);
                        multiplyMatrixVector(ievc, ws.tmpPostBottom, ws.rotatedPost);
                    }

                    ComplexBlockKernelUtils.applyPlanToOuterProduct(
                            plan, ws.rotatedPre, ws.rotatedPost, (wp * wc) / denom, eigenBasisAccum);
                }
            }
        }
    }

    /**
     * Same as {@link #accumulateBranchContribution} but uses the packed real-eigenvalue
     * kernel (RealKernelUtils) for models whose eigensystem is all-real.
     */
    private void accumulateBranchContributionReal(int childNumber,
                                                   double branchLength,
                                                   double[] categoryRates,
                                                   double[] categoryWeights,
                                                   double[] patternWeights,
                                                   EigenDecomposition eigenDecomp,
                                                   RealKernelUtils.RealKernelPlan plan,
                                                   double[] evec,
                                                   double[] ievc,
                                                   double[] eigenBasisAccum,
                                                   ThreadWorkspace ws,
                                                   boolean useInternalRotatedMessages) {
        final int categoryCount = likelihoodDelegate.getCategoryCount();
        final int patternCount  = likelihoodDelegate.getPatternCount();

        for (int c = 0; c < categoryCount; c++) {
            final double wc = categoryWeights[c];
            final double tc = branchLength * categoryRates[c];

            if (patternCount == 1) {
                final double wp = patternWeights[0];
                if (wp == 0.0) continue;

                double denom = 0.0;
                if (useInternalRotatedMessages) {
                    likelihoodDelegate.getInternalPreOrderBranchTopInto(childNumber, c, 0, ws.rotatedPre);
                    likelihoodDelegate.getPostOrderBranchTopInto(childNumber, c, 0, ws.tmpPostTop);
                    likelihoodDelegate.getInternalPostOrderBranchBottomInto(childNumber, c, 0, ws.rotatedPost);
                    for (int s = 0; s < stateCount; s++) {
                        denom += ws.rotatedPre[s] * ws.tmpPostTop[s];
                    }
                } else {
                    likelihoodDelegate.getPreOrderBranchTopInto(childNumber, c, 0, ws.tmpPreTop);
                    likelihoodDelegate.getPreOrderBranchBottomInto(childNumber, c, 0, ws.tmpPreBottom);
                    likelihoodDelegate.getPostOrderBranchBottomInto(childNumber, c, 0, ws.tmpPostBottom);
                    for (int s = 0; s < stateCount; s++) {
                        denom += ws.tmpPreBottom[s] * ws.tmpPostBottom[s];
                    }
                }
                if (denom <= 0.0 || Double.isNaN(denom) || Double.isInfinite(denom)) continue;

                if (!useInternalRotatedMessages) {
                    multiplyTransposeMatrixVector(evec, ws.tmpPreTop, ws.rotatedPre);
                    multiplyMatrixVector(ievc, ws.tmpPostBottom, ws.rotatedPost);
                }

                if (useInternalRotatedMessages &&
                        likelihoodDelegate.borrowRealBranchExponentials(childNumber, tc, ws.realBranchExponentials)) {
                    RealKernelUtils.applyToOuterProductFromExponentials(
                            plan,
                            ws.realBranchExponentials.values(),
                            ws.realBranchExponentials.offset(),
                            tc,
                            ws.rotatedPre,
                            ws.rotatedPost,
                            (wp * wc) / denom,
                            eigenBasisAccum,
                            stateCount);
                } else {
                    RealKernelUtils.fillAndApplyToOuterProduct(
                            plan, eigenDecomp, tc, ws.rotatedPre, ws.rotatedPost,
                            (wp * wc) / denom, eigenBasisAccum, stateCount);
                }
            } else {
                if (useInternalRotatedMessages &&
                        likelihoodDelegate.borrowRealBranchExponentials(childNumber, tc, ws.realBranchExponentials)) {
                    RealKernelUtils.fillCoefficientsFromExponentials(
                            plan,
                            ws.realBranchExponentials.values(),
                            ws.realBranchExponentials.offset(),
                            tc,
                            stateCount);
                } else {
                    RealKernelUtils.fillCoefficients(plan, eigenDecomp, tc, stateCount);
                }

                for (int p = 0; p < patternCount; p++) {
                    final double wp = patternWeights[p];
                    if (wp == 0.0) continue;

                    double denom = 0.0;
                    if (useInternalRotatedMessages) {
                        likelihoodDelegate.getInternalPreOrderBranchTopInto(childNumber, c, p, ws.rotatedPre);
                        likelihoodDelegate.getPostOrderBranchTopInto(childNumber, c, p, ws.tmpPostTop);
                        likelihoodDelegate.getInternalPostOrderBranchBottomInto(childNumber, c, p, ws.rotatedPost);
                        for (int s = 0; s < stateCount; s++) {
                            denom += ws.rotatedPre[s] * ws.tmpPostTop[s];
                        }
                    } else {
                        likelihoodDelegate.getPreOrderBranchTopInto(childNumber, c, p, ws.tmpPreTop);
                        likelihoodDelegate.getPreOrderBranchBottomInto(childNumber, c, p, ws.tmpPreBottom);
                        likelihoodDelegate.getPostOrderBranchBottomInto(childNumber, c, p, ws.tmpPostBottom);
                        for (int s = 0; s < stateCount; s++) {
                            denom += ws.tmpPreBottom[s] * ws.tmpPostBottom[s];
                        }
                    }
                    if (denom <= 0.0 || Double.isNaN(denom) || Double.isInfinite(denom)) continue;

                    if (!useInternalRotatedMessages) {
                        multiplyTransposeMatrixVector(evec, ws.tmpPreTop, ws.rotatedPre);
                        multiplyMatrixVector(ievc, ws.tmpPostBottom, ws.rotatedPost);
                    }

                    RealKernelUtils.applyToOuterProduct(
                            plan, ws.rotatedPre, ws.rotatedPost, (wp * wc) / denom, eigenBasisAccum, stateCount);
                }
            }
        }
    }

    private void addRawLogRateScoresFromEigenBasis(double[] eigenBasisAccum,
                                                   double[] evec,
                                                   double[] ievc,
                                                   double[] destination,
                                                   int offset,
                                                   double[] eigenRowProjection,
                                                   double[] standardRowProjection) {
        final int K = stateCount;

        for (int source = 0; source < K; source++) {
            Arrays.fill(eigenRowProjection, 0.0);

            for (int eigenRow = 0; eigenRow < K; eigenRow++) {
                final double left = ievc[eigenRow * K + source];
                final int accumOffset = eigenRow * K;
                for (int eigenCol = 0; eigenCol < K; eigenCol++) {
                    eigenRowProjection[eigenCol] += left * eigenBasisAccum[accumOffset + eigenCol];
                }
            }

            for (int destinationState = 0; destinationState < K; destinationState++) {
                final int columnOffset = destinationState * K;
                double total = 0.0;
                for (int eigenCol = 0; eigenCol < K; eigenCol++) {
                    total += eigenRowProjection[eigenCol] * evec[columnOffset + eigenCol];
                }
                standardRowProjection[destinationState] = total;
            }

            for (int destinationState = 0; destinationState < K; destinationState++) {
                if (source == destinationState) {
                    continue;
                }
                destination[offset + directedRateIndex[source][destinationState]] +=
                        standardRowProjection[destinationState] - standardRowProjection[source];
            }
        }
    }

    private void finishRawLogRateScores(double[] destination, List<SubstitutionModel> models) {
        final int rateDimension = getDirectedRateDimension();
        final double[] generator = generatorBuffer;

        for (int modelNumber = 0; modelNumber < substitutionModelCount; modelNumber++) {
            final SubstitutionModel substitutionModel = models.get(modelNumber);
            if (!isDirectLogRateModel(substitutionModel)) {
                continue;
            }

            final LogRateSubstitutionModel model = (LogRateSubstitutionModel) substitutionModel;
            final LogAdditiveCtmcRateProvider.DataAugmented rateProvider =
                    (LogAdditiveCtmcRateProvider.DataAugmented) model.getRateProvider();
            final Parameter transformedParameter = rateProvider.getLogRateParameter();
            final boolean normalize = model.getNormalization();
            final boolean scaleByFrequencies = model.getScaleRatesByFrequencies();
            final Transform transform = model.getTransform();
            final double normalizationScalar = normalize && transform != null ? 1 / model.setupMatrix() : 0.0;

            model.getInfinitesimalMatrix(generator);

            for (int nodeIndex = 0; nodeIndex < tree.getNodeCount(); nodeIndex++) {
                if (tree.isRoot(tree.getNode(nodeIndex))) {
                    continue;
                }

                final int offset = (nodeIndex * substitutionModelCount + modelNumber) * rateDimension;
                final double normalizationConstant = computeNormalizationConstant(destination, offset, generator);

                for (int source = 0; source < stateCount; source++) {
                    for (int destinationState = 0; destinationState < stateCount; destinationState++) {
                        if (source == destinationState) {
                            continue;
                        }

                        final int index = directedRateIndex[source][destinationState];
                        final int ij = source * stateCount + destinationState;
                        final double sourceFrequency = model.getFrequencyModel().getFrequency(source);
                        double element;
                        if (transform == null) {
                            element = generator[ij];
                        } else {
                            element = transform.gradient(transformedParameter.getParameterValue(index));
                            if (normalize) {
                                element *= normalizationScalar;
                            }
                            if (scaleByFrequencies) {
                                element *= sourceFrequency;
                            }
                        }

                        double score = destination[offset + index] * element;
                        if (normalize) {
                            score -= element * sourceFrequency * normalizationConstant;
                        }
                        destination[offset + index] = score;
                    }
                }
            }
        }
    }

    private double computeNormalizationConstant(double[] rawScores, int offset, double[] generator) {
        double total = 0.0;
        for (int source = 0; source < stateCount; source++) {
            for (int destinationState = 0; destinationState < stateCount; destinationState++) {
                if (source == destinationState) {
                    continue;
                }
                final int index = directedRateIndex[source][destinationState];
                total += rawScores[offset + index] * generator[source * stateCount + destinationState];
            }
        }
        return total;
    }

    private boolean isDirectLogRateModel(SubstitutionModel model) {
        if (!(model instanceof LogRateSubstitutionModel)) {
            return false;
        }

        final LogAdditiveCtmcRateProvider rateProvider =
                ((LogRateSubstitutionModel) model).getRateProvider();
        if (!(rateProvider instanceof LogAdditiveCtmcRateProvider.DataAugmented)) {
            return false;
        }

        final Parameter parameter =
                ((LogAdditiveCtmcRateProvider.DataAugmented) rateProvider).getLogRateParameter();
        return parameter != null && parameter.getDimension() == getDirectedRateDimension();
    }

    private int getDirectedRateDimension() {
        return stateCount * (stateCount - 1);
    }

    private static int[][] makeDirectedRateIndex(int stateCount) {
        final int[][] index = new int[stateCount][stateCount];
        for (int i = 0; i < stateCount; i++) {
            Arrays.fill(index[i], -1);
        }

        int k = 0;
        for (int i = 0; i < stateCount; ++i) {
            for (int j = i + 1; j < stateCount; ++j) {
                index[i][j] = k++;
            }
        }

        for (int j = 0; j < stateCount; ++j) {
            for (int i = j + 1; i < stateCount; ++i) {
                index[i][j] = k++;
            }
        }

        return index;
    }

    /**
     * Computes: out[modelOffset .. modelOffset+K²) = R^{-T} * eigenBasisAccum * R^T
     *
     * Step 1: mid = R^{-T} * accum
     *   mid[i,col] = sum_k R^{-T}[i,k] * accum[k,col]
     *              = sum_k ievc[k*K+i] * accum[k*K+col]
     *
     * Step 2: out = mid * R^T
     *   out[row,col] = sum_j mid[row,j] * R^T[j,col]
     *                = sum_j mid[row*K+j] * evec[col*K+j]
     */
    private void rotateIntoOutput(double[] eigenBasisAccum,
                                  double[] evec,
                                  double[] ievc,
                                  double[] out,
                                  int modelOffset) {
        rotateIntoOutput(eigenBasisAccum, evec, ievc, out, modelOffset, midBuffer);
    }

    private void rotateIntoOutput(double[] eigenBasisAccum,
                                  double[] evec,
                                  double[] ievc,
                                  double[] out,
                                  int modelOffset,
                                  double[] rotationBuffer) {
        final int K = stateCount;

        // Step 1: mid = R^{-T} * accum
        for (int row = 0; row < K; row++) {
            final int rowOff = row * K;
            Arrays.fill(rotationBuffer, rowOff, rowOff + K, 0.0);
            for (int k = 0; k < K; k++) {
                final double a = ievc[k * K + row];
                final int accumOff = k * K;
                for (int col = 0; col < K; col++) {
                    rotationBuffer[rowOff + col] += a * eigenBasisAccum[accumOff + col];
                }
            }
        }

        // Step 2: out = mid * R^T
        for (int row = 0; row < K; row++) {
            final int rowOff = row * K;
            final int outRowOff = modelOffset + rowOff;
            for (int col = 0; col < K; col++) {
                double sum = 0.0;
                final int colOff = col * K;
                for (int j = 0; j < K; j++) {
                    sum += rotationBuffer[rowOff + j] * evec[colOff + j];
                }
                out[outRowOff + col] = sum;
            }
        }
    }

    private static double relativeWeight(int k, double[] weights) {
        double sum = 0.0;
        for (double w : weights) {
            sum += w;
        }
        return weights[k] / sum;
    }

    /** out[j] += matrix[i*K+j] * vector[i]  for all i  (M^T * v) */
    private void multiplyTransposeMatrixVector(double[] matrix, double[] vector, double[] out) {
        Arrays.fill(out, 0.0);
        for (int i = 0; i < stateCount; i++) {
            final double vi = vector[i];
            if (vi == 0.0) continue;
            final int rowBase = i * stateCount;
            for (int j = 0; j < stateCount; j++) {
                out[j] += matrix[rowBase + j] * vi;
            }
        }
    }

    /** out[i] = sum_j matrix[i*K+j] * vector[j]  (M * v) */
    private void multiplyMatrixVector(double[] matrix, double[] vector, double[] out) {
        for (int i = 0; i < stateCount; i++) {
            double sum = 0.0;
            final int rowBase = i * stateCount;
            for (int j = 0; j < stateCount; j++) {
                sum += matrix[rowBase + j] * vector[j];
            }
            out[i] = sum;
        }
    }
}

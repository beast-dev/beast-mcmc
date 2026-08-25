package dr.inference.operators;

import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.UncertainSiteList;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.RewardMixtureBranchRateModel;
import dr.evomodel.branchratemodel.RewardMixtureCategoricalBranchRateModel;
import dr.evomodel.branchratemodel.RewardsAwareCategoricalMixtureBranchRatesDynamic;
import dr.evomodel.branchratemodel.RewardsAwareMixtureBranchRates;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evomodel.substmodel.EigenDecomposition;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.ProcessSimulation;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.discrete.discretetreedataLikelihood.ComplexBlockKernelUtils;
import dr.evomodel.treedatalikelihood.discrete.discretetreedataLikelihood.RealKernelUtils;
import dr.evomodel.treedatalikelihood.preorder.AbstractBeagleGradientDelegate;
import dr.evomodel.treedatalikelihood.preorder.DiscretePartialsType;
import dr.inference.hmc.EmbeddedOrdinalParameter;
import dr.inference.model.Parameter;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.List;
import java.util.Locale;

/**
 * BEAGLE-backed CTMC edge-evidence provider for dependent reward processes.
 *
 * Postorder messages are read from BEAGLE, branch-start top messages are
 * reconstructed locally by default, and an optional BEAGLE-preorder mode can
 * use BEAGLE TOP preorder buffers directly when the linked runtime exposes the
 * v5 preorder entry point. Candidate kernels are evaluated locally so all
 * continuous and atomic candidates for a branch share the same cached messages.
 */
public final class BeagleRewardDependentCtmcEdgeEvidenceProvider
        implements RewardDependentCtmcEdgeEvidenceProvider {

    private static final String PROPERTY_PREFIX =
            "dr.inference.operators.BeagleRewardDependentCtmcEdgeEvidenceProvider.";
    private static final double ATOM_MATCH_TOLERANCE = 1.0e-12;

    private final TreeDataLikelihood treeDataLikelihood;
    private final BeagleDataLikelihoodDelegate likelihoodDelegate;
    private final Tree tree;
    private final BranchModel branchModel;
    private final RewardMixtureBranchRateModel rewardBranchRates;
    private final RewardStateAdapter rewardStateAdapter;
    private final SiteRateModel siteRateModel;
    private final BeaglePreOrderDelegate preOrderTopDelegate;
    private final ProcessSimulation preOrderTopSimulation;
    private final BeaglePreOrderDelegate preOrderBottomDelegate;
    private final ProcessSimulation preOrderBottomSimulation;

    private final int stateCount;
    private final int patternCount;
    private final int categoryCount;
    private final int flattenedLength;

    private final double[] prePartials;
    private final double[] postPartials;
    private final double[] rotatedPre;
    private final double[] rotatedPost;
    private final double[] transitionMatrix;
    private final double[][] topPartialsByNode;
    private final double[][] preBottomPartialsByNode;
    private final double[][] postPartialsByNode;
    private final double[][] postTopPartialsByNode;
    private final double[][] beagleTopPartialsByNode;
    private final double[][] beagleBottomPartialsByNode;
    private final boolean[] topPartialsKnown;
    private final boolean[] postPartialsKnown;
    private final double[] diagnosticMatrix;
    private RewardStateSnapshot baselineRewardState;

    private EigenDecomposition[] eigenDecompositions;
    private boolean[] allRealEigen;
    private ComplexBlockKernelUtils.ComplexKernelPlan[] complexPlans;

    private final Diagnostics diagnostics;
    private boolean beaglePreOrderAvailable = false;
    private String beaglePreOrderFailure = "";
    private long prepareCount = 0L;
    private long diagnosticRowsWritten = 0L;
    private double baselineLogLikelihood = Double.NaN;

    public BeagleRewardDependentCtmcEdgeEvidenceProvider(final TreeDataLikelihood treeDataLikelihood) {
        this(treeDataLikelihood, Diagnostics.fromSystemProperties());
    }

    public BeagleRewardDependentCtmcEdgeEvidenceProvider(final TreeDataLikelihood treeDataLikelihood,
                                                        final Diagnostics diagnostics) {
        if (treeDataLikelihood == null) {
            throw new IllegalArgumentException("treeDataLikelihood must be non-null");
        }
        this.diagnostics = diagnostics == null ? Diagnostics.disabled() : diagnostics;

        final DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();
        if (!(delegate instanceof BeagleDataLikelihoodDelegate)) {
            throw new IllegalArgumentException(
                    "Dependent CTMC reward evidence requires BeagleDataLikelihoodDelegate, found " +
                            (delegate == null ? "null" : delegate.getClass().getName())
            );
        }

        this.treeDataLikelihood = treeDataLikelihood;
        this.likelihoodDelegate = (BeagleDataLikelihoodDelegate) delegate;
        this.tree = treeDataLikelihood.getTree();
        this.branchModel = likelihoodDelegate.getBranchModel();
        this.siteRateModel = likelihoodDelegate.getSiteRateModel();

        final BranchRateModel branchRateModel = treeDataLikelihood.getBranchRateModel();
        if (!(branchRateModel instanceof RewardMixtureBranchRateModel)) {
            throw new IllegalArgumentException(
                    "Dependent CTMC reward evidence requires a RewardMixtureBranchRateModel, found " +
                            branchRateModel.getClass().getName()
            );
        }
        this.rewardBranchRates = (RewardMixtureBranchRateModel) branchRateModel;
        this.rewardStateAdapter = createRewardStateAdapter(rewardBranchRates);

        if (!likelihoodDelegate.isUsePreOrder()) {
            throw new IllegalArgumentException(
                    "Dependent CTMC TreeDataLikelihood must be configured with usePreOrder=\"true\""
            );
        }
        if (likelihoodDelegate.getPreOrderSettings().isUseSpectralRepresentation()) {
            throw new IllegalArgumentException(
                    "Dependent CTMC reward evidence requires a non-spectral TreeDataLikelihood; " +
                            "set useSpectralRepresentation=\"false\" for dependent CTMC likelihoods."
            );
        }
        // BEAGLE's preorder partial buffers use the same state, pattern, and category layout
        // exposed by the delegate's PatternList and SiteRateModel.
        this.stateCount = likelihoodDelegate.getPatternList().getDataType().getStateCount();
        this.patternCount = likelihoodDelegate.getPatternList().getPatternCount();
        this.categoryCount = likelihoodDelegate.getSiteRateModel().getCategoryCount();
        this.flattenedLength = stateCount * patternCount * categoryCount;

        this.prePartials = new double[flattenedLength];
        this.postPartials = new double[flattenedLength];
        this.rotatedPre = new double[stateCount];
        this.rotatedPost = new double[stateCount];
        this.transitionMatrix = new double[stateCount * stateCount];
        this.topPartialsByNode = new double[tree.getNodeCount()][flattenedLength];
        this.postPartialsByNode = new double[tree.getNodeCount()][flattenedLength];
        this.topPartialsKnown = new boolean[tree.getNodeCount()];
        this.postPartialsKnown = new boolean[tree.getNodeCount()];
        // preBottomPartialsByNode, postTopPartialsByNode, and beagleBottomPartialsByNode are read
        // only from writeDiagnosticRow, and beagleTopPartialsByNode only from the BEAGLE-preorder
        // evidence/comparison paths; at real alignment scale (hundreds of nodes x thousands of
        // patterns) each full table is tens to hundreds of MB, so skip allocating the ones this
        // provider's configuration will never populate rather than paying for them unconditionally.
        this.preBottomPartialsByNode = this.diagnostics.enabled
                ? new double[tree.getNodeCount()][flattenedLength] : null;
        this.postTopPartialsByNode = this.diagnostics.enabled
                ? new double[tree.getNodeCount()][flattenedLength] : null;
        this.beagleTopPartialsByNode = (this.diagnostics.compareBeaglePreOrder || this.diagnostics.useBeaglePreOrderEvidence)
                ? new double[tree.getNodeCount()][flattenedLength] : null;
        this.beagleBottomPartialsByNode = this.diagnostics.enabled
                ? new double[tree.getNodeCount()][flattenedLength] : null;
        this.diagnosticMatrix = new double[stateCount * stateCount * categoryCount];

        final String providerName =
                treeDataLikelihood.getId() == null ? "dependentRewardEvidence" : treeDataLikelihood.getId();
        this.preOrderTopDelegate = new BeagleTopPreOrderDelegate(
                providerName,
                tree,
                likelihoodDelegate
        );
        this.preOrderTopSimulation = new ProcessSimulation(treeDataLikelihood, preOrderTopDelegate);
        this.preOrderBottomDelegate = new BeagleBottomPreOrderDelegate(
                providerName,
                tree,
                likelihoodDelegate
        );
        this.preOrderBottomSimulation = new ProcessSimulation(treeDataLikelihood, preOrderBottomDelegate);

        refreshSpectralStructures();
    }

    @Override
    public void prepare() {
        final long prepareStart = RewardMixturePerformanceStats.startTimer();
        try {
            baselineLogLikelihood = treeDataLikelihood.getLogLikelihood();
            baselineRewardState = rewardStateAdapter.snapshot();

            long start = RewardMixturePerformanceStats.startTimer();
            refreshSpectralStructures();
            RewardMixturePerformanceStats.addBeagleRefreshSpectralNanos(
                    RewardMixturePerformanceStats.elapsed(start));

            resetLocalMessageCacheFlags();

            if (useLazyLocalMessages()) {
                start = RewardMixturePerformanceStats.startTimer();
                fillRootTopPartials();
                RewardMixturePerformanceStats.addBeagleFillTopNanos(
                        RewardMixturePerformanceStats.elapsed(start));
            } else {
                start = RewardMixturePerformanceStats.startTimer();
                fillPostPartialsForAllNodes();
                RewardMixturePerformanceStats.addBeagleFillPostNanos(
                        RewardMixturePerformanceStats.elapsed(start));

                if (postTopPartialsByNode != null) {
                    start = RewardMixturePerformanceStats.startTimer();
                    fillPostTopPartialsForAllNodes();
                    RewardMixturePerformanceStats.addBeagleFillPostTopNanos(
                            RewardMixturePerformanceStats.elapsed(start));
                }

                start = RewardMixturePerformanceStats.startTimer();
                fillTopPartialsFromRoot();
                RewardMixturePerformanceStats.addBeagleFillTopNanos(
                        RewardMixturePerformanceStats.elapsed(start));

                if (preBottomPartialsByNode != null) {
                    start = RewardMixturePerformanceStats.startTimer();
                    fillPreBottomPartialsForAllNodes();
                    RewardMixturePerformanceStats.addBeagleFillPreBottomNanos(
                            RewardMixturePerformanceStats.elapsed(start));
                }
            }

            beaglePreOrderAvailable = false;
            beaglePreOrderFailure = "";
            if (diagnostics.compareBeaglePreOrder || diagnostics.useBeaglePreOrderEvidence) {
                start = RewardMixturePerformanceStats.startTimer();
                fillBeaglePreOrderPartialsIfAvailable();
                RewardMixturePerformanceStats.addBeagleFillPreOrderNanos(
                        RewardMixturePerformanceStats.elapsed(start));
            }
            prepareCount++;
        } finally {
            RewardMixturePerformanceStats.recordBeaglePrepare(
                    RewardMixturePerformanceStats.elapsed(prepareStart));
        }
    }

    @Override
    public double logEvidence(final int branchNodeNumber, final double rawReward) {
        final long start = RewardMixturePerformanceStats.startTimer();
        final boolean atomicCandidate = RewardMixturePerformanceStats.ENABLED && matchingAtomState(rawReward) >= 0;
        try {
            validateBranchNodeNumber(branchNodeNumber);
            if (diagnostics.enabled || !diagnostics.useBeaglePreOrderEvidence) {
                ensureLocalMessagesForBranch(branchNodeNumber);
            }

            final double manualLogEvidence = diagnostics.enabled || !diagnostics.useBeaglePreOrderEvidence
                    ? logEvidenceFromCachedMessages(
                            branchNodeNumber, rawReward,
                            topPartialsByNode[branchNodeNumber],
                            postPartialsByNode[branchNodeNumber])
                    : Double.NaN;

            final double logEvidence;
            if (diagnostics.useBeaglePreOrderEvidence) {
                if (!beaglePreOrderAvailable) {
                    throw new IllegalStateException("BEAGLE preorder evidence requested but unavailable: " +
                            beaglePreOrderFailure);
                }
                logEvidence = logEvidenceFromCachedMessages(
                        branchNodeNumber, rawReward,
                        beagleTopPartialsByNode[branchNodeNumber],
                        postPartialsByNode[branchNodeNumber]);
            } else {
                logEvidence = manualLogEvidence;
            }

            if (diagnostics.enabled) {
                writeDiagnosticRow(branchNodeNumber, rawReward, manualLogEvidence);
            }

            return logEvidence;
        } finally {
            RewardMixturePerformanceStats.recordBeagleLogEvidence(
                    RewardMixturePerformanceStats.elapsed(start), atomicCandidate);
        }
    }

    private double logEvidenceFromCachedMessages(final int branchNodeNumber,
                                                 final double rawReward,
                                                 final double[] topPartials,
                                                 final double[] postPartials) {
        final long evidenceStart = RewardMixturePerformanceStats.startTimer();
        try {
            final NodeRef node = validateBranchNodeNumber(branchNodeNumber);

            final BranchModel.Mapping mapping = branchModel.getBranchModelMapping(node);
            final int[] order = mapping.getOrder();
            if (order.length != 1) {
                throw new UnsupportedOperationException(
                        "Dependent CTMC reward evidence currently supports exactly one substitution model per branch; " +
                                "branch " + branchNodeNumber + " has " + order.length
                );
            }

            final int modelNumber = order[0];
            final double modelWeight = relativeWeight(0, mapping.getWeights());
            final double branchRate = rewardBranchRates.getBranchRateForRawReward(tree, node, rawReward);
            final double candidateBranchLength = tree.getBranchLength(node) * branchRate * modelWeight;

            final long copyStart = RewardMixturePerformanceStats.startTimer();
            System.arraycopy(topPartials, 0, prePartials, 0, flattenedLength);
            System.arraycopy(postPartials, 0, this.postPartials, 0, flattenedLength);
            RewardMixturePerformanceStats.recordBeagleMessageCopy(
                    RewardMixturePerformanceStats.elapsed(copyStart));

            final double[] patternWeights = likelihoodDelegate.getPatternList().getPatternWeights();
            final double[] categoryWeights = siteRateModel.getCategoryProportions();
            final double[] categoryRates = siteRateModel.getCategoryRates();

            double logEvidence = 0.0;

            for (int p = 0; p < patternCount; p++) {
                final double wp = patternWeights[p];
                if (wp == 0.0) {
                    continue;
                }

                double patternEvidence = 0.0;
                for (int c = 0; c < categoryCount; c++) {
                    final double wc = categoryWeights[c];
                    if (wc == 0.0) {
                        continue;
                    }
                    final double rate = categoryRates == null ? 1.0 : categoryRates[c];
                    final int offset = ((c * patternCount) + p) * stateCount;
                    final double time = candidateBranchLength * rate;
                    final double inner = edgeInnerProduct(modelNumber, time, offset);
                    patternEvidence += wc * inner;
                }

                if (!(patternEvidence > 0.0) || Double.isNaN(patternEvidence)) {
                    return Double.NEGATIVE_INFINITY;
                }
                logEvidence += wp * Math.log(patternEvidence);
            }

            return logEvidence;
        } finally {
            RewardMixturePerformanceStats.recordBeagleCachedMessageEvidence(
                    RewardMixturePerformanceStats.elapsed(evidenceStart));
        }
    }

    private boolean useLazyLocalMessages() {
        return !diagnostics.enabled &&
                !diagnostics.compareBeaglePreOrder &&
                !diagnostics.useBeaglePreOrderEvidence;
    }

    private void resetLocalMessageCacheFlags() {
        Arrays.fill(postPartialsKnown, false);
        Arrays.fill(topPartialsKnown, false);
    }

    private void ensureLocalMessagesForBranch(final int branchNodeNumber) {
        ensurePostPartialsForNode(branchNodeNumber);
        ensureTopPartialsForNode(branchNodeNumber);
    }

    private void ensurePostPartialsForNode(final int nodeNumber) {
        if (postPartialsKnown[nodeNumber]) {
            return;
        }
        final long start = RewardMixturePerformanceStats.startTimer();
        try {
            fillPostPartials(tree.getNode(nodeNumber), postPartialsByNode[nodeNumber]);
            postPartialsKnown[nodeNumber] = true;
        } finally {
            RewardMixturePerformanceStats.addBeagleFillPostNanos(
                    RewardMixturePerformanceStats.elapsed(start));
        }
    }

    private void ensureTopPartialsForNode(final int nodeNumber) {
        if (topPartialsKnown[nodeNumber]) {
            return;
        }

        final NodeRef node = tree.getNode(nodeNumber);
        if (tree.isRoot(node)) {
            final long start = RewardMixturePerformanceStats.startTimer();
            try {
                fillRootTopPartials();
            } finally {
                RewardMixturePerformanceStats.addBeagleFillTopNanos(
                        RewardMixturePerformanceStats.elapsed(start));
            }
            return;
        }

        final NodeRef parent = tree.getParent(node);
        ensureTopPartialsForNode(parent.getNumber());

        final long start = RewardMixturePerformanceStats.startTimer();
        try {
            fillTopPartialsForChild(parent, node);
        } finally {
            RewardMixturePerformanceStats.addBeagleFillTopNanos(
                    RewardMixturePerformanceStats.elapsed(start));
        }
    }

    private NodeRef validateBranchNodeNumber(final int branchNodeNumber) {
        if (branchNodeNumber < 0 || branchNodeNumber >= tree.getNodeCount()) {
            throw new IllegalArgumentException("branchNodeNumber out of range: " + branchNodeNumber);
        }
        final NodeRef node = tree.getNode(branchNodeNumber);
        if (tree.isRoot(node)) {
            throw new IllegalArgumentException("Root node has no branch: " + branchNodeNumber);
        }
        return node;
    }

    private void fillPostPartialsForAllNodes() {
        for (int i = 0; i < tree.getNodeCount(); i++) {
            fillPostPartials(tree.getNode(i), postPartialsByNode[i]);
            postPartialsKnown[i] = true;
        }
    }

    private void fillPostTopPartialsForAllNodes() {
        for (int i = 0; i < tree.getNodeCount(); i++) {
            Arrays.fill(postTopPartialsByNode[i], 0.0);
            final NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node)) {
                propagateBranchBottomToTop(node, postPartialsByNode[i], postTopPartialsByNode[i]);
            }
        }
    }

    private void fillPostPartials(final NodeRef node, final double[] out) {
        if (!tree.isExternal(node)) {
            likelihoodDelegate.getPartials(node.getNumber(), out);
            for (int i = 0; i < flattenedLength; i++) {
                out[i] = messageMagnitude(out[i]);
            }
            return;
        }

        final PatternList patternList = likelihoodDelegate.getPatternList();
        final String taxonId = tree.getNodeTaxon(node).getId();
        final int taxonIndex = patternList.getTaxonIndex(taxonId);
        if (taxonIndex < 0) {
            throw new IllegalArgumentException("Taxon " + taxonId + " is not found in pattern list " +
                    patternList.getId());
        }

        int offset = 0;
        for (int p = 0; p < patternCount; p++) {
            if (patternList instanceof UncertainSiteList) {
                ((UncertainSiteList) patternList).fillPartials(taxonIndex, p, out, offset);
                offset += stateCount;
            } else if (patternList.areUncertain()) {
                final double[] probabilities = patternList.getUncertainPatternState(taxonIndex, p);
                System.arraycopy(probabilities, 0, out, offset, stateCount);
                offset += stateCount;
            } else {
                final int state = patternList.getPatternState(taxonIndex, p);
                final boolean[] stateSet = patternList.getDataType().getStateSet(state);
                for (int s = 0; s < stateCount; s++) {
                    out[offset++] = stateSet[s] ? 1.0 : 0.0;
                }
            }
        }

        final int categoryBlockLength = patternCount * stateCount;
        for (int c = 1; c < categoryCount; c++) {
            System.arraycopy(out, 0, out, c * categoryBlockLength, categoryBlockLength);
        }
    }

    private void fillTopPartialsFromRoot() {
        for (int i = 0; i < tree.getNodeCount(); i++) {
            Arrays.fill(topPartialsByNode[i], 0.0);
        }

        fillRootTopPartials();
        fillTopPartialsBelow(tree.getRoot());
    }

    private void fillRootTopPartials() {
        final NodeRef root = tree.getRoot();
        final double[] rootFrequencies = likelihoodDelegate.getEvolutionaryProcessDelegate().getRootStateFrequencies();
        final double[] rootTop = topPartialsByNode[root.getNumber()];

        int offset = 0;
        for (int c = 0; c < categoryCount; c++) {
            for (int p = 0; p < patternCount; p++) {
                System.arraycopy(rootFrequencies, 0, rootTop, offset, stateCount);
                offset += stateCount;
            }
        }
        topPartialsKnown[root.getNumber()] = true;
    }

    private void fillPreBottomPartialsForAllNodes() {
        for (int i = 0; i < tree.getNodeCount(); i++) {
            Arrays.fill(preBottomPartialsByNode[i], 0.0);
            final NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node)) {
                propagateBranchTopToBottom(node, topPartialsByNode[i], preBottomPartialsByNode[i]);
            }
        }
    }

    private void fillBeaglePreOrderPartialsIfAvailable() {
        try {
            preOrderTopSimulation.cacheSimulatedTraits(null);
            for (int i = 0; i < tree.getNodeCount(); i++) {
                final double[] out = beagleTopPartialsByNode[i];
                preOrderTopDelegate.getPreorderPartials(i, DiscretePartialsType.TOP, out);
                for (int j = 0; j < flattenedLength; j++) {
                    out[j] = messageMagnitude(out[j]);
                }
            }

            if (beagleBottomPartialsByNode != null) {
                preOrderBottomSimulation.cacheSimulatedTraits(null);
                for (int i = 0; i < tree.getNodeCount(); i++) {
                    final double[] out = beagleBottomPartialsByNode[i];
                    preOrderBottomDelegate.getPreorderPartials(i, DiscretePartialsType.BOTTOM, out);
                    for (int j = 0; j < flattenedLength; j++) {
                        out[j] = messageMagnitude(out[j]);
                    }
                }
            }
            beaglePreOrderAvailable = true;
        } catch (UnsatisfiedLinkError e) {
            beaglePreOrderFailure = e.getClass().getName() + ": " + e.getMessage();
            beaglePreOrderAvailable = false;
        } catch (RuntimeException e) {
            beaglePreOrderFailure = e.getClass().getName() + ": " + e.getMessage();
            beaglePreOrderAvailable = false;
        }
    }

    private void propagateBranchBottomToTop(final NodeRef node,
                                            final double[] bottom,
                                            final double[] top) {
        final BranchModel.Mapping mapping = branchModel.getBranchModelMapping(node);
        final int[] order = mapping.getOrder();
        if (order.length != 1) {
            throw new UnsupportedOperationException(
                    "Dependent CTMC reward evidence currently supports exactly one substitution model per branch; " +
                            "branch " + node.getNumber() + " has " + order.length
            );
        }

        final int modelNumber = order[0];
        final double modelWeight = relativeWeight(0, mapping.getWeights());
        final double branchLength = tree.getBranchLength(node) *
                rewardBranchRates.getBranchRate(tree, node) * modelWeight;
        final SubstitutionModel model = branchModel.getSubstitutionModels().get(modelNumber);
        final double[] categoryRates = siteRateModel.getCategoryRates();

        for (int c = 0; c < categoryCount; c++) {
            final double rate = categoryRates == null ? 1.0 : categoryRates[c];
            model.getTransitionProbabilities(branchLength * rate, transitionMatrix);

            for (int p = 0; p < patternCount; p++) {
                final int offset = ((c * patternCount) + p) * stateCount;
                for (int parentState = 0; parentState < stateCount; parentState++) {
                    final int row = parentState * stateCount;
                    double sum = 0.0;
                    for (int childState = 0; childState < stateCount; childState++) {
                        sum += transitionMatrix[row + childState] * bottom[offset + childState];
                    }
                    top[offset + parentState] = messageMagnitude(sum);
                }
            }
        }
    }

    private void propagateBranchTopToBottom(final NodeRef node,
                                            final double[] top,
                                            final double[] bottom) {
        final BranchModel.Mapping mapping = branchModel.getBranchModelMapping(node);
        final int[] order = mapping.getOrder();
        if (order.length != 1) {
            throw new UnsupportedOperationException(
                    "Dependent CTMC reward evidence currently supports exactly one substitution model per branch; " +
                            "branch " + node.getNumber() + " has " + order.length
            );
        }

        final int modelNumber = order[0];
        final double modelWeight = relativeWeight(0, mapping.getWeights());
        final double branchLength = tree.getBranchLength(node) *
                rewardBranchRates.getBranchRate(tree, node) * modelWeight;
        final SubstitutionModel model = branchModel.getSubstitutionModels().get(modelNumber);
        final double[] categoryRates = siteRateModel.getCategoryRates();

        for (int c = 0; c < categoryCount; c++) {
            final double rate = categoryRates == null ? 1.0 : categoryRates[c];
            model.getTransitionProbabilities(branchLength * rate, transitionMatrix);

            for (int p = 0; p < patternCount; p++) {
                final int offset = ((c * patternCount) + p) * stateCount;
                for (int childState = 0; childState < stateCount; childState++) {
                    double sum = 0.0;
                    for (int parentState = 0; parentState < stateCount; parentState++) {
                        sum += top[offset + parentState] *
                                transitionMatrix[parentState * stateCount + childState];
                    }
                    bottom[offset + childState] = messageMagnitude(sum);
                }
            }
        }
    }

    private void fillTopPartialsBelow(final NodeRef parent) {
        final int childCount = tree.getChildCount(parent);
        for (int i = 0; i < childCount; i++) {
            final NodeRef child = tree.getChild(parent, i);
            fillTopPartialsForChild(parent, child);
            if (!tree.isExternal(child)) {
                fillTopPartialsBelow(child);
            }
        }
    }

    private void fillTopPartialsForChild(final NodeRef parent, final NodeRef child) {
        final double[] childTop = topPartialsByNode[child.getNumber()];
        fillParentNodeContext(parent, childTop);

        final int childCount = tree.getChildCount(parent);
        for (int i = 0; i < childCount; i++) {
            final NodeRef sibling = tree.getChild(parent, i);
            if (sibling.getNumber() != child.getNumber()) {
                multiplyBySiblingContribution(sibling, childTop);
            }
        }
        topPartialsKnown[child.getNumber()] = true;
    }

    private void fillParentNodeContext(final NodeRef parent, final double[] out) {
        final double[] parentTop = topPartialsByNode[parent.getNumber()];
        if (tree.isRoot(parent)) {
            System.arraycopy(parentTop, 0, out, 0, flattenedLength);
            return;
        }

        final BranchModel.Mapping mapping = branchModel.getBranchModelMapping(parent);
        final int[] order = mapping.getOrder();
        if (order.length != 1) {
            throw new UnsupportedOperationException(
                    "Dependent CTMC reward evidence currently supports exactly one substitution model per branch; " +
                            "parent branch " + parent.getNumber() + " has " + order.length
            );
        }

        final int modelNumber = order[0];
        final double modelWeight = relativeWeight(0, mapping.getWeights());
        final double branchLength = tree.getBranchLength(parent) *
                rewardBranchRates.getBranchRate(tree, parent) * modelWeight;
        final SubstitutionModel model = branchModel.getSubstitutionModels().get(modelNumber);
        final double[] categoryRates = siteRateModel.getCategoryRates();

        for (int c = 0; c < categoryCount; c++) {
            final double rate = categoryRates == null ? 1.0 : categoryRates[c];
            model.getTransitionProbabilities(branchLength * rate, transitionMatrix);

            for (int p = 0; p < patternCount; p++) {
                final int offset = ((c * patternCount) + p) * stateCount;
                for (int endState = 0; endState < stateCount; endState++) {
                    double sum = 0.0;
                    for (int startState = 0; startState < stateCount; startState++) {
                        sum += parentTop[offset + startState] *
                                transitionMatrix[startState * stateCount + endState];
                    }
                    out[offset + endState] = messageMagnitude(sum);
                }
            }
        }
    }

    private void multiplyBySiblingContribution(final NodeRef sibling, final double[] top) {
        ensurePostPartialsForNode(sibling.getNumber());

        final BranchModel.Mapping mapping = branchModel.getBranchModelMapping(sibling);
        final int[] order = mapping.getOrder();
        if (order.length != 1) {
            throw new UnsupportedOperationException(
                    "Dependent CTMC reward evidence currently supports exactly one substitution model per branch; " +
                            "sibling branch " + sibling.getNumber() + " has " + order.length
            );
        }

        final int modelNumber = order[0];
        final double modelWeight = relativeWeight(0, mapping.getWeights());
        final double branchLength = tree.getBranchLength(sibling) *
                rewardBranchRates.getBranchRate(tree, sibling) * modelWeight;
        final SubstitutionModel model = branchModel.getSubstitutionModels().get(modelNumber);
        final double[] categoryRates = siteRateModel.getCategoryRates();
        final double[] siblingPost = postPartialsByNode[sibling.getNumber()];

        for (int c = 0; c < categoryCount; c++) {
            final double rate = categoryRates == null ? 1.0 : categoryRates[c];
            model.getTransitionProbabilities(branchLength * rate, transitionMatrix);

            for (int p = 0; p < patternCount; p++) {
                final int offset = ((c * patternCount) + p) * stateCount;
                for (int parentState = 0; parentState < stateCount; parentState++) {
                    final int row = parentState * stateCount;
                    double contribution = 0.0;
                    for (int childState = 0; childState < stateCount; childState++) {
                        contribution += transitionMatrix[row + childState] *
                                siblingPost[offset + childState];
                    }
                    top[offset + parentState] = messageMagnitude(top[offset + parentState]) *
                            messageMagnitude(contribution);
                }
            }
        }
    }

    private void writeDiagnosticRow(final int branchNodeNumber,
                                    final double rawReward,
                                    final double manualLogEvidence) {
        if (diagnosticRowsWritten >= diagnostics.maxRows ||
                diagnosticRowsWritten >= diagnostics.maxRowsPerProvider) {
            return;
        }
        if (diagnostics.maxBranches > 0 && branchNodeNumber >= diagnostics.maxBranches) {
            return;
        }

        final NodeRef node = tree.getNode(branchNodeNumber);
        if (tree.isRoot(node)) {
            return;
        }

        final int parameterIndex = rewardBranchRates.getParameterIndexFromNode(node);
        final int atomState = matchingAtomState(rawReward);
        final String candidateKind = atomState >= 0 ? "atomic" : "continuous";
        final int matrixBufferIndex =
                likelihoodDelegate.getEvolutionaryProcessDelegate().getMatrixIndex(branchNodeNumber);

        final double beagleLogEvidence = beaglePreOrderAvailable
                ? logEvidenceFromCachedMessages(
                        branchNodeNumber, rawReward,
                        beagleTopPartialsByNode[branchNodeNumber],
                        postPartialsByNode[branchNodeNumber])
                : Double.NaN;

        final ExactComparison exact = diagnostics.compareExact
                ? exactComparison(branchNodeNumber, rawReward, manualLogEvidence, beagleLogEvidence)
                : ExactComparison.notComputed();

        final StringBuilder sb = new StringBuilder(4096);
        sb.append(nullSafeId(treeDataLikelihood.getId())).append('\t');
        sb.append(prepareCount).append('\t');
        sb.append(diagnosticRowsWritten).append('\t');
        sb.append(branchNodeNumber).append('\t');
        sb.append(parameterIndex).append('\t');
        sb.append(candidateKind).append('\t');
        sb.append(atomState).append('\t');
        appendDouble(sb, rawReward).append('\t');
        sb.append(matrixBufferIndex).append('\t');
        appendDouble(sb, manualLogEvidence).append('\t');
        appendDouble(sb, beagleLogEvidence).append('\t');
        sb.append(beaglePreOrderAvailable ? "available" : beaglePreOrderFailure).append('\t');
        appendDouble(sb, exact.exactLogLikelihood).append('\t');
        appendDouble(sb, exact.exactDelta).append('\t');
        appendDouble(sb, exact.manualDelta).append('\t');
        appendDouble(sb, exact.beagleDelta).append('\t');
        appendDouble(sb, exact.manualMinusExactDelta).append('\t');
        appendDouble(sb, exact.beagleMinusExactDelta).append('\t');
        sb.append(formatArray(postPartialsByNode[branchNodeNumber])).append('\t');
        sb.append(formatArray(topPartialsByNode[branchNodeNumber])).append('\t');
        sb.append(beaglePreOrderAvailable ? formatArray(beagleTopPartialsByNode[branchNodeNumber]) : "NA").append('\t');
        sb.append(formatArray(preBottomPartialsByNode[branchNodeNumber])).append('\t');
        sb.append(formatArray(postTopPartialsByNode[branchNodeNumber])).append('\t');
        sb.append(beaglePreOrderAvailable ? formatArray(beagleBottomPartialsByNode[branchNodeNumber]) : "NA").append('\t');
        sb.append(formatCandidateTransitionMatrices(branchNodeNumber, rawReward)).append('\t');
        sb.append(formatCurrentBeagleTransitionMatrix(matrixBufferIndex));

        diagnostics.writeLine(diagnosticRowsWritten == 0L, sb.toString());
        diagnosticRowsWritten++;
    }

    private ExactComparison exactComparison(final int branchNodeNumber,
                                           final double rawReward,
                                           final double manualLogEvidence,
                                           final double beagleLogEvidence) {
        final double currentRawReward = currentRawReward(branchNodeNumber);
        final double manualCurrentEvidence = logEvidenceFromCachedMessages(
                branchNodeNumber, currentRawReward,
                topPartialsByNode[branchNodeNumber],
                postPartialsByNode[branchNodeNumber]);
        final double manualDelta = finiteDifference(manualLogEvidence, manualCurrentEvidence);

        final double beagleDelta;
        if (beaglePreOrderAvailable) {
            final double beagleCurrentEvidence = logEvidenceFromCachedMessages(
                    branchNodeNumber, currentRawReward,
                    beagleTopPartialsByNode[branchNodeNumber],
                    postPartialsByNode[branchNodeNumber]);
            beagleDelta = finiteDifference(beagleLogEvidence, beagleCurrentEvidence);
        } else {
            beagleDelta = Double.NaN;
        }

        final double exactLogLikelihood = exactLogLikelihoodForCandidate(branchNodeNumber, rawReward);
        final double exactDelta = finiteDifference(exactLogLikelihood, baselineLogLikelihood);

        return new ExactComparison(
                exactLogLikelihood,
                exactDelta,
                manualDelta,
                beagleDelta,
                finiteDifference(manualDelta, exactDelta),
                finiteDifference(beagleDelta, exactDelta)
        );
    }

    private double exactLogLikelihoodForCandidate(final int branchNodeNumber, final double rawReward) {
        final NodeRef node = tree.getNode(branchNodeNumber);
        final int parameterIndex = rewardBranchRates.getParameterIndexFromNode(node);
        final int atomState = matchingAtomState(rawReward);

        final RewardStateSnapshot current = rewardStateAdapter.snapshot();
        try {
            baselineRewardState.restore();
            rewardStateAdapter.setCandidate(parameterIndex, rawReward, atomState);
            treeDataLikelihood.makeDirty();
            return treeDataLikelihood.getLogLikelihood();
        } finally {
            current.restore();
            treeDataLikelihood.makeDirty();
            treeDataLikelihood.getLogLikelihood();
        }
    }

    private double currentRawReward(final int branchNodeNumber) {
        final NodeRef node = tree.getNode(branchNodeNumber);
        final int parameterIndex = rewardBranchRates.getParameterIndexFromNode(node);
        return baselineRewardState.rawRewardForParameterIndex(parameterIndex);
    }

    private int matchingAtomState(final double rawReward) {
        final int stateCount = rewardBranchRates.getRewardRates().getStateIndices().getDimension();
        for (int state = 0; state < stateCount; state++) {
            final double atomReward = rewardBranchRates.getRawRewardForAtomState(state);
            if (Math.abs(atomReward - rawReward) <= ATOM_MATCH_TOLERANCE) {
                return state;
            }
        }
        return -1;
    }

    private String formatCandidateTransitionMatrices(final int branchNodeNumber, final double rawReward) {
        final NodeRef node = tree.getNode(branchNodeNumber);
        final BranchModel.Mapping mapping = branchModel.getBranchModelMapping(node);
        final int[] order = mapping.getOrder();
        if (order.length != 1) {
            return "unsupportedMapping";
        }
        final int modelNumber = order[0];
        final double modelWeight = relativeWeight(0, mapping.getWeights());
        final double branchRate = rewardBranchRates.getBranchRateForRawReward(tree, node, rawReward);
        final double candidateBranchLength = tree.getBranchLength(node) * branchRate * modelWeight;
        final SubstitutionModel model = branchModel.getSubstitutionModels().get(modelNumber);
        final double[] categoryRates = siteRateModel.getCategoryRates();

        int offset = 0;
        for (int c = 0; c < categoryCount; c++) {
            final double rate = categoryRates == null ? 1.0 : categoryRates[c];
            model.getTransitionProbabilities(candidateBranchLength * rate, transitionMatrix);
            System.arraycopy(transitionMatrix, 0, diagnosticMatrix, offset, transitionMatrix.length);
            offset += transitionMatrix.length;
        }
        return formatArray(diagnosticMatrix, 0, offset);
    }

    private String formatCurrentBeagleTransitionMatrix(final int matrixBufferIndex) {
        try {
            likelihoodDelegate.getBeagleInstance().getTransitionMatrix(matrixBufferIndex, diagnosticMatrix);
            return formatArray(diagnosticMatrix, 0, categoryCount * stateCount * stateCount);
        } catch (RuntimeException e) {
            return e.getClass().getSimpleName() + ":" + e.getMessage();
        } catch (UnsatisfiedLinkError e) {
            return e.getClass().getSimpleName() + ":" + e.getMessage();
        }
    }

    private static double finiteDifference(final double x, final double y) {
        if (Double.isNaN(x) || Double.isNaN(y)) {
            return Double.NaN;
        }
        if (Double.isInfinite(x) || Double.isInfinite(y)) {
            return x - y;
        }
        return x - y;
    }

    private static String formatArray(final double[] values) {
        return formatArray(values, 0, values.length);
    }

    private static String formatArray(final double[] values, final int offset, final int length) {
        final StringBuilder sb = new StringBuilder(length * 12 + 2);
        sb.append('[');
        for (int i = 0; i < length; i++) {
            if (i > 0) {
                sb.append(',');
            }
            appendDouble(sb, values[offset + i]);
        }
        sb.append(']');
        return sb.toString();
    }

    private static StringBuilder appendDouble(final StringBuilder sb, final double value) {
        if (Double.isNaN(value)) {
            sb.append("NaN");
        } else if (value == Double.POSITIVE_INFINITY) {
            sb.append("Infinity");
        } else if (value == Double.NEGATIVE_INFINITY) {
            sb.append("-Infinity");
        } else {
            sb.append(String.format(Locale.US, "%.17g", value));
        }
        return sb;
    }

    private static String nullSafeId(final String id) {
        return id == null ? "TreeDataLikelihood" : id;
    }

    private static double messageMagnitude(final double value) {
        // Scaled BEAGLE buffers can expose signed entries; edge evidence needs
        // likelihood-message magnitudes when building local branch contexts.
        return Double.isNaN(value) ? 0.0 : Math.abs(value);
    }

    private double edgeInnerProduct(final int modelNumber, final double time, final int offset) {
        final long start = RewardMixturePerformanceStats.startTimer();
        try {
            final EigenDecomposition eigen = eigenDecompositions[modelNumber];
            if (eigen == null) {
                return directEdgeInnerProduct(modelNumber, time, offset);
            }

            rotatePre(eigen.getEigenVectors(), prePartials, offset, rotatedPre);
            rotatePost(eigen.getInverseEigenVectors(), postPartials, offset, rotatedPost);

            final double spectralInner;
            if (allRealEigen[modelNumber]) {
                spectralInner = realSpectralInnerProduct(eigen, time, rotatedPre, rotatedPost);
            } else {
                final ComplexBlockKernelUtils.ComplexKernelPlan plan = complexPlans[modelNumber];
                ComplexBlockKernelUtils.fillTransitionCoefficients(plan, eigen, time);
                spectralInner = ComplexBlockKernelUtils.blockDiagonalTransitionInnerProduct(
                        plan, rotatedPre, 0, rotatedPost, 0
                );
            }

            if (spectralInner > 0.0 && !Double.isNaN(spectralInner)) {
                return spectralInner;
            }

            return directEdgeInnerProduct(modelNumber, time, offset);
        } finally {
            RewardMixturePerformanceStats.recordBeagleEdgeInnerProduct(
                    RewardMixturePerformanceStats.elapsed(start));
        }
    }

    private double directEdgeInnerProduct(final int modelNumber, final double time, final int offset) {
        final SubstitutionModel model = branchModel.getSubstitutionModels().get(modelNumber);
        model.getTransitionProbabilities(time, transitionMatrix);
        return RewardsMixtureBranchResamplingHelper.bilinearFormStable(
                prePartials, offset, transitionMatrix, postPartials, offset, stateCount
        );
    }

    private void refreshSpectralStructures() {
        final List<SubstitutionModel> models = branchModel.getSubstitutionModels();
        if (eigenDecompositions == null || eigenDecompositions.length != models.size()) {
            eigenDecompositions = new EigenDecomposition[models.size()];
            allRealEigen = new boolean[models.size()];
            complexPlans = new ComplexBlockKernelUtils.ComplexKernelPlan[models.size()];
        }

        for (int i = 0; i < models.size(); i++) {
            final SubstitutionModel model = models.get(i);
            final EigenDecomposition eigen = model.getEigenDecomposition();
            eigenDecompositions[i] = eigen;
            if (eigen == null) {
                allRealEigen[i] = false;
                complexPlans[i] = null;
                continue;
            }

            allRealEigen[i] = RealKernelUtils.isAllReal(eigen, stateCount);
            if (!allRealEigen[i]) {
                if (complexPlans[i] == null) {
                    complexPlans[i] = new ComplexBlockKernelUtils.ComplexKernelPlan(stateCount);
                }
                ComplexBlockKernelUtils.fillStructure(complexPlans[i], eigen, stateCount);
            }
        }
    }

    private static double realSpectralInnerProduct(final EigenDecomposition eigen,
                                                   final double time,
                                                   final double[] pre,
                                                   final double[] post) {
        final double[] eigenValues = eigen.getEigenValues();
        double sum = 0.0;
        for (int i = 0; i < pre.length; i++) {
            sum += pre[i] * Math.exp(time * eigenValues[i]) * post[i];
        }
        return sum;
    }

    private void rotatePre(final double[] eigenVectors,
                           final double[] source,
                           final int offset,
                           final double[] out) {
        for (int i = 0; i < stateCount; i++) {
            double sum = 0.0;
            for (int s = 0; s < stateCount; s++) {
                sum += eigenVectors[s * stateCount + i] * source[offset + s];
            }
            out[i] = sum;
        }
    }

    private void rotatePost(final double[] inverseEigenVectors,
                            final double[] source,
                            final int offset,
                            final double[] out) {
        for (int i = 0; i < stateCount; i++) {
            double sum = 0.0;
            final int row = i * stateCount;
            for (int s = 0; s < stateCount; s++) {
                sum += inverseEigenVectors[row + s] * source[offset + s];
            }
            out[i] = sum;
        }
    }

    private static double relativeWeight(final int index, final double[] weights) {
        double sum = 0.0;
        for (double weight : weights) {
            sum += weight;
        }
        return weights[index] / sum;
    }

    private static final class ExactComparison {
        final double exactLogLikelihood;
        final double exactDelta;
        final double manualDelta;
        final double beagleDelta;
        final double manualMinusExactDelta;
        final double beagleMinusExactDelta;

        private ExactComparison(final double exactLogLikelihood,
                                final double exactDelta,
                                final double manualDelta,
                                final double beagleDelta,
                                final double manualMinusExactDelta,
                                final double beagleMinusExactDelta) {
            this.exactLogLikelihood = exactLogLikelihood;
            this.exactDelta = exactDelta;
            this.manualDelta = manualDelta;
            this.beagleDelta = beagleDelta;
            this.manualMinusExactDelta = manualMinusExactDelta;
            this.beagleMinusExactDelta = beagleMinusExactDelta;
        }

        private static ExactComparison notComputed() {
            return new ExactComparison(Double.NaN, Double.NaN, Double.NaN,
                    Double.NaN, Double.NaN, Double.NaN);
        }
    }

    private static RewardStateAdapter createRewardStateAdapter(final RewardMixtureBranchRateModel branchRates) {
        if (branchRates instanceof RewardsAwareMixtureBranchRates) {
            return new LegacyRewardStateAdapter((RewardsAwareMixtureBranchRates) branchRates);
        }
        if (branchRates instanceof RewardMixtureCategoricalBranchRateModel) {
            return new CategoricalRewardStateAdapter((RewardMixtureCategoricalBranchRateModel) branchRates);
        }
        throw new IllegalArgumentException(
                "Unsupported reward-mixture branch-rate state class: " + branchRates.getClass().getName());
    }

    private interface RewardStateAdapter {
        RewardStateSnapshot snapshot();

        void setCandidate(int parameterIndex, double rawReward, int atomState);
    }

    private interface RewardStateSnapshot {
        void restore();

        double rawRewardForParameterIndex(int parameterIndex);
    }

    private static final class LegacyRewardStateAdapter implements RewardStateAdapter {
        private final RewardsAwareMixtureBranchRates branchRates;

        private LegacyRewardStateAdapter(final RewardsAwareMixtureBranchRates branchRates) {
            this.branchRates = branchRates;
        }

        @Override
        public RewardStateSnapshot snapshot() {
            return new LegacyRewardStateSnapshot(branchRates);
        }

        @Override
        public void setCandidate(final int parameterIndex, final double rawReward, final int atomState) {
            if (atomState >= 0) {
                branchRates.getAtomIndices().setParameterValue(parameterIndex, atomState);
                branchRates.getIndicator().setParameterValue(parameterIndex, 1.0);
            } else {
                branchRates.getRateParameter().setParameterValue(parameterIndex, rawReward);
                branchRates.getIndicator().setParameterValue(parameterIndex, 0.0);
            }
        }
    }

    private static final class LegacyRewardStateSnapshot implements RewardStateSnapshot {
        private final RewardsAwareMixtureBranchRates branchRates;
        final double[] indicators;
        final double[] atomIndices;
        final double[] ctsRewards;

        private LegacyRewardStateSnapshot(final RewardsAwareMixtureBranchRates branchRates) {
            this.branchRates = branchRates;
            final int dim = branchRates.getRateParameter().getDimension();
            this.indicators = new double[dim];
            this.atomIndices = new double[dim];
            this.ctsRewards = new double[dim];
            final Parameter indicator = branchRates.getIndicator();
            final Parameter atoms = branchRates.getAtomIndices();
            final Parameter cts = branchRates.getRateParameter();
            for (int i = 0; i < dim; i++) {
                indicators[i] = indicator.getParameterValue(i);
                atomIndices[i] = atoms.getParameterValue(i);
                ctsRewards[i] = cts.getParameterValue(i);
            }
        }

        @Override
        public void restore() {
            final Parameter indicator = branchRates.getIndicator();
            final Parameter atoms = branchRates.getAtomIndices();
            final Parameter cts = branchRates.getRateParameter();
            for (int i = 0; i < indicators.length; i++) {
                cts.setParameterValue(i, ctsRewards[i]);
                atoms.setParameterValue(i, atomIndices[i]);
                indicator.setParameterValue(i, indicators[i]);
            }
        }

        @Override
        public double rawRewardForParameterIndex(final int parameterIndex) {
            final double indicator = indicators[parameterIndex];
            if (Math.abs(indicator - 1.0) <= 1.0e-9) {
                return branchRates.getRawRewardForAtomState((int) Math.round(atomIndices[parameterIndex]));
            }
            return ctsRewards[parameterIndex];
        }
    }

    private static final class CategoricalRewardStateAdapter implements RewardStateAdapter {
        private final RewardMixtureCategoricalBranchRateModel branchRates;

        private CategoricalRewardStateAdapter(final RewardMixtureCategoricalBranchRateModel branchRates) {
            this.branchRates = branchRates;
        }

        @Override
        public RewardStateSnapshot snapshot() {
            return new CategoricalRewardStateSnapshot(branchRates);
        }

        @Override
        public void setCandidate(final int parameterIndex, final double rawReward, final int atomState) {
            if (atomState >= 0) {
                branchRates.getCategoryParameter().setParameterValue(
                        parameterIndex,
                        representativeValueForCategory(atomState + 1));
            } else {
                branchRates.getRateParameter().setParameterValue(parameterIndex, rawReward);
                branchRates.getCategoryParameter().setParameterValue(
                        parameterIndex,
                        representativeValueForCategory(0));
            }
        }

        private double representativeValueForCategory(final int category) {
            if (branchRates instanceof RewardsAwareCategoricalMixtureBranchRatesDynamic) {
                throw new UnsupportedOperationException(
                        "Exact-perturbation diagnostics (dependentCtmcCompareExact / dependentCtmcDiagnostics) " +
                                "are not yet supported with dynamic reward-category ordering: which axis bucket " +
                                "a category occupies is branch-specific there, and this candidate-setting path " +
                                "does not thread a branch index through. Disable exact-comparison diagnostics " +
                                "for rewardsAwareCategoricalMixtureBranchRatesDynamic runs.");
            }
            final Parameter cuts = branchRates.getCategoryCutParameter();
            if (category < 0 || category + 1 >= cuts.getDimension()) {
                throw new IllegalArgumentException("Reward category out of range: " + category);
            }
            final double lower = cuts.getParameterValue(category);
            final double upper = cuts.getParameterValue(category + 1);
            if (Double.isFinite(lower) && Double.isFinite(upper)) {
                return lower + 0.5 * (upper - lower);
            }
            if (Double.isFinite(lower)) {
                return lower + 1.0;
            }
            if (Double.isFinite(upper)) {
                return upper - 1.0;
            }
            return category;
        }
    }

    private static final class CategoricalRewardStateSnapshot implements RewardStateSnapshot {
        private final RewardMixtureCategoricalBranchRateModel branchRates;
        private final double[] categoryState;
        private final double[] ctsRewards;
        private final double[] categoryCuts;

        private CategoricalRewardStateSnapshot(final RewardMixtureCategoricalBranchRateModel branchRates) {
            this.branchRates = branchRates;
            this.categoryState = branchRates.getCategoryParameter().getParameterValues();
            this.ctsRewards = branchRates.getRateParameter().getParameterValues();
            this.categoryCuts = branchRates.getCategoryCutParameter().getParameterValues();
        }

        @Override
        public void restore() {
            final Parameter category = branchRates.getCategoryParameter();
            final Parameter cts = branchRates.getRateParameter();
            for (int i = 0; i < categoryState.length; i++) {
                category.setParameterValue(i, categoryState[i]);
                cts.setParameterValue(i, ctsRewards[i]);
            }
        }

        @Override
        public double rawRewardForParameterIndex(final int parameterIndex) {
            if (branchRates instanceof RewardsAwareCategoricalMixtureBranchRatesDynamic) {
                throw new UnsupportedOperationException(
                        "Exact-perturbation diagnostics (dependentCtmcCompareExact / dependentCtmcDiagnostics) " +
                                "are not yet supported with dynamic reward-category ordering: decoding a raw " +
                                "axis value to a category is branch-specific there. Disable exact-comparison " +
                                "diagnostics for rewardsAwareCategoricalMixtureBranchRatesDynamic runs.");
            }
            final int category = new EmbeddedOrdinalParameter(categoryCuts).getStateIndex(
                    categoryState[parameterIndex]);
            if (category == 0) {
                return ctsRewards[parameterIndex];
            }
            return branchRates.getRawRewardForAtomState(category - 1);
        }
    }

    public static final class Diagnostics {
        public static final String ENABLED = "enabled";
        public static final String FILE_NAME = "fileName";
        public static final String COMPARE_BEAGLE_PREORDER = "compareBeaglePreorder";
        public static final String COMPARE_EXACT = "compareExact";
        public static final String USE_BEAGLE_PREORDER_EVIDENCE = "useBeaglePreorderEvidence";
        public static final String MAX_BRANCHES = "maxBranches";
        public static final String MAX_ROWS = "maxRows";

        private static final String DEFAULT_FILE_NAME =
                "dependent_ctmc_edge_evidence_diagnostics.tsv";
        private static final String HEADER =
                "treeDataLikelihood\tprepare\trow\tbranchNode\tbranchParameter\tcandidateKind\tatomState" +
                        "\trawReward\tmatrixBufferIndex\tmanualLogEvidence\tbeagleLogEvidence" +
                        "\tbeaglePreorderStatus\texactLogLikelihood\texactDelta\tmanualDelta\tbeagleDelta" +
                        "\tmanualMinusExactDelta\tbeagleMinusExactDelta\tpostPartials\tmanualTopPartials" +
                        "\tbeagleTopPartials\tmanualPreBottomPartials\tmanualPostTopPartials" +
                        "\tbeagleBottomPartials\tcandidateTransitionMatrices\tbeagleCurrentTransitionMatrices";

        final boolean enabled;
        final String fileName;
        final boolean compareBeaglePreOrder;
        final boolean compareExact;
        final boolean useBeaglePreOrderEvidence;
        final int maxBranches;
        final long maxRows;
        final long maxRowsPerProvider;

        private Diagnostics(final boolean enabled,
                            final String fileName,
                            final boolean compareBeaglePreOrder,
                            final boolean compareExact,
                            final boolean useBeaglePreOrderEvidence,
                            final int maxBranches,
                            final long maxRows,
                            final long maxRowsPerProvider) {
            this.enabled = enabled;
            this.fileName = fileName;
            this.compareBeaglePreOrder = enabled && compareBeaglePreOrder;
            this.compareExact = enabled && compareExact;
            this.useBeaglePreOrderEvidence = useBeaglePreOrderEvidence;
            this.maxBranches = maxBranches;
            this.maxRows = maxRows;
            this.maxRowsPerProvider = maxRowsPerProvider;
        }

        public static Diagnostics disabled() {
            return new Diagnostics(false, null, false, false, false, 0, 0L, 0L);
        }

        public static Diagnostics fromSystemProperties() {
            final boolean enabled = Boolean.parseBoolean(
                    System.getProperty(PROPERTY_PREFIX + ENABLED, "false"));
            final boolean useBeagle = Boolean.parseBoolean(
                    System.getProperty(PROPERTY_PREFIX + USE_BEAGLE_PREORDER_EVIDENCE, "false"));
            if (!enabled && !useBeagle) {
                return disabled();
            }
            final String fileName = System.getProperty(PROPERTY_PREFIX + FILE_NAME, DEFAULT_FILE_NAME);
            final boolean compareBeagle = Boolean.parseBoolean(
                    System.getProperty(PROPERTY_PREFIX + COMPARE_BEAGLE_PREORDER, "true"));
            final boolean compareExact = Boolean.parseBoolean(
                    System.getProperty(PROPERTY_PREFIX + COMPARE_EXACT, "true"));
            final int maxBranches = Integer.parseInt(
                    System.getProperty(PROPERTY_PREFIX + MAX_BRANCHES, "2147483647"));
            final long maxRows = Long.parseLong(
                    System.getProperty(PROPERTY_PREFIX + MAX_ROWS, "9223372036854775807"));
            return new Diagnostics(enabled, fileName, compareBeagle, compareExact, useBeagle,
                    maxBranches, maxRows, maxRows);
        }

        public static Diagnostics create(final boolean enabled,
                                         final String fileName,
                                         final boolean compareBeaglePreOrder,
                                         final boolean compareExact,
                                         final int maxBranches,
                                         final long maxRows) {
            return create(
                    enabled,
                    fileName,
                    compareBeaglePreOrder,
                    compareExact,
                    false,
                    maxBranches,
                    maxRows
            );
        }

        public static Diagnostics create(final boolean enabled,
                                         final String fileName,
                                         final boolean compareBeaglePreOrder,
                                         final boolean compareExact,
                                         final boolean useBeaglePreOrderEvidence,
                                         final int maxBranches,
                                         final long maxRows) {
            if (!enabled) {
                if (!useBeaglePreOrderEvidence) {
                    return disabled();
                }
                return new Diagnostics(
                        false,
                        null,
                        false,
                        false,
                        true,
                        maxBranches <= 0 ? Integer.MAX_VALUE : maxBranches,
                        maxRows <= 0L ? Long.MAX_VALUE : maxRows,
                        maxRows <= 0L ? Long.MAX_VALUE : maxRows
                );
            }
            final String effectiveFileName =
                    fileName == null || fileName.length() == 0 ? DEFAULT_FILE_NAME : fileName;
            return new Diagnostics(true, effectiveFileName, compareBeaglePreOrder, compareExact,
                    useBeaglePreOrderEvidence,
                    maxBranches <= 0 ? Integer.MAX_VALUE : maxBranches,
                    maxRows <= 0L ? Long.MAX_VALUE : maxRows,
                    maxRows <= 0L ? Long.MAX_VALUE : maxRows);
        }

        private void writeLine(final boolean firstProviderRow, final String line) {
            if (!enabled) {
                return;
            }
            if (fileName == null || fileName.length() == 0) {
                synchronized (System.err) {
                    if (firstProviderRow) {
                        System.err.println(HEADER);
                    }
                    System.err.println(line);
                }
                return;
            }

            final File file = new File(fileName);
            final boolean writeHeader = !file.exists() || file.length() == 0L;
            try {
                final PrintWriter out = new PrintWriter(new FileWriter(file, true));
                try {
                    if (writeHeader) {
                        out.println(HEADER);
                    }
                    out.println(line);
                } finally {
                    out.close();
                }
            } catch (IOException e) {
                throw new RuntimeException("Unable to write dependent CTMC edge-evidence diagnostics to " +
                        fileName, e);
            }
        }
    }

    private abstract static class BeaglePreOrderDelegate extends AbstractBeagleGradientDelegate {

        private final String traitName;

        private BeaglePreOrderDelegate(final String name,
                                       final Tree tree,
                                       final BeagleDataLikelihoodDelegate likelihoodDelegate,
                                       final String meaning) {
            super(name + ".dependentRewardEvidencePreOrder." + meaning,
                    tree, likelihoodDelegate);
            this.traitName = name + ".dependentRewardEvidencePreOrder." +
                    meaning + ".preOrderTouch";
        }

        @Override
        protected int getGradientLength() {
            return 0;
        }

        @Override
        protected void getNodeDerivatives(final Tree tree, final double[] first, final double[] second) {
            // This delegate only materializes preorder partials.
        }

        @Override
        protected void constructTraits(final TreeTraitProvider.Helper treeTraitHelper) {
            treeTraitHelper.addTrait(new TreeTrait.DA() {
                @Override
                public String getTraitName() {
                    return traitName;
                }

                @Override
                public Intent getIntent() {
                    return Intent.WHOLE_TREE;
                }

                @Override
                public double[] getTrait(final Tree tree, final NodeRef node) {
                    return new double[0];
                }
            });
        }
    }

    private static final class BeagleTopPreOrderDelegate extends BeaglePreOrderDelegate {

        private BeagleTopPreOrderDelegate(final String name,
                                          final Tree tree,
                                          final BeagleDataLikelihoodDelegate likelihoodDelegate) {
            super(name, tree, likelihoodDelegate, DiscretePartialsType.TOP.getMeaning());
        }

        @Override
        protected DiscretePartialsType getPreOrderType() {
            return DiscretePartialsType.TOP;
        }
    }

    private static final class BeagleBottomPreOrderDelegate extends BeaglePreOrderDelegate {

        private BeagleBottomPreOrderDelegate(final String name,
                                             final Tree tree,
                                             final BeagleDataLikelihoodDelegate likelihoodDelegate) {
            super(name, tree, likelihoodDelegate, DiscretePartialsType.BOTTOM.getMeaning());
        }

        @Override
        protected DiscretePartialsType getPreOrderType() {
            return DiscretePartialsType.BOTTOM;
        }
    }
}

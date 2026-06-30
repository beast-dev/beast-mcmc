package dr.inference.operators;

import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.UncertainSiteList;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchratemodel.BranchRateModel;
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
    private final RewardsAwareMixtureBranchRates rewardBranchRates;
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
    private final double[] diagnosticMatrix;
    private final double[] baselineIndicators;
    private final double[] baselineAtomIndices;
    private final double[] baselineCtsRewards;

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
        if (!(branchRateModel instanceof RewardsAwareMixtureBranchRates)) {
            throw new IllegalArgumentException(
                    "Dependent CTMC reward evidence requires RewardsAwareMixtureBranchRates, found " +
                            branchRateModel.getClass().getName()
            );
        }
        this.rewardBranchRates = (RewardsAwareMixtureBranchRates) branchRateModel;

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
        this.preBottomPartialsByNode = new double[tree.getNodeCount()][flattenedLength];
        this.postPartialsByNode = new double[tree.getNodeCount()][flattenedLength];
        this.postTopPartialsByNode = new double[tree.getNodeCount()][flattenedLength];
        this.beagleTopPartialsByNode = new double[tree.getNodeCount()][flattenedLength];
        this.beagleBottomPartialsByNode = new double[tree.getNodeCount()][flattenedLength];
        this.diagnosticMatrix = new double[stateCount * stateCount * categoryCount];
        final int rewardParameterCount = rewardBranchRates.getRateParameter().getDimension();
        this.baselineIndicators = new double[rewardParameterCount];
        this.baselineAtomIndices = new double[rewardParameterCount];
        this.baselineCtsRewards = new double[rewardParameterCount];

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
        baselineLogLikelihood = treeDataLikelihood.getLogLikelihood();
        snapshotRewardState(baselineIndicators, baselineAtomIndices, baselineCtsRewards);
        refreshSpectralStructures();
        fillPostPartialsForAllNodes();
        fillPostTopPartialsForAllNodes();
        fillTopPartialsFromRoot();
        fillPreBottomPartialsForAllNodes();
        beaglePreOrderAvailable = false;
        beaglePreOrderFailure = "";
        if (diagnostics.compareBeaglePreOrder || diagnostics.useBeaglePreOrderEvidence) {
            fillBeaglePreOrderPartialsIfAvailable();
        }
        prepareCount++;
    }

    @Override
    public double logEvidence(final int branchNodeNumber, final double rawReward) {
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
    }

    private double logEvidenceFromCachedMessages(final int branchNodeNumber,
                                                 final double rawReward,
                                                 final double[] topPartials,
                                                 final double[] postPartials) {
        if (branchNodeNumber < 0 || branchNodeNumber >= tree.getNodeCount()) {
            throw new IllegalArgumentException("branchNodeNumber out of range: " + branchNodeNumber);
        }
        final NodeRef node = tree.getNode(branchNodeNumber);
        if (tree.isRoot(node)) {
            throw new IllegalArgumentException("Root node has no branch: " + branchNodeNumber);
        }

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

        System.arraycopy(topPartials, 0, prePartials, 0, flattenedLength);
        System.arraycopy(postPartials, 0, this.postPartials, 0, flattenedLength);

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
    }

    private void fillPostPartialsForAllNodes() {
        for (int i = 0; i < tree.getNodeCount(); i++) {
            fillPostPartials(tree.getNode(i), postPartialsByNode[i]);
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

        fillTopPartialsBelow(root);
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

            preOrderBottomSimulation.cacheSimulatedTraits(null);
            for (int i = 0; i < tree.getNodeCount(); i++) {
                final double[] out = beagleBottomPartialsByNode[i];
                preOrderBottomDelegate.getPreorderPartials(i, DiscretePartialsType.BOTTOM, out);
                for (int j = 0; j < flattenedLength; j++) {
                    out[j] = messageMagnitude(out[j]);
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
            if (sibling != child) {
                multiplyBySiblingContribution(sibling, childTop);
            }
        }
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

        final FullParameterSnapshot current = new FullParameterSnapshot(rewardBranchRates);
        try {
            restoreRewardState(baselineIndicators, baselineAtomIndices, baselineCtsRewards);
            if (atomState >= 0) {
                rewardBranchRates.getAtomIndices().setParameterValue(parameterIndex, atomState);
                rewardBranchRates.getIndicator().setParameterValue(parameterIndex, 1.0);
            } else {
                rewardBranchRates.getRateParameter().setParameterValue(parameterIndex, rawReward);
                rewardBranchRates.getIndicator().setParameterValue(parameterIndex, 0.0);
            }
            treeDataLikelihood.makeDirty();
            return treeDataLikelihood.getLogLikelihood();
        } finally {
            current.restore(rewardBranchRates);
            treeDataLikelihood.makeDirty();
            treeDataLikelihood.getLogLikelihood();
        }
    }

    private double currentRawReward(final int branchNodeNumber) {
        final NodeRef node = tree.getNode(branchNodeNumber);
        final int parameterIndex = rewardBranchRates.getParameterIndexFromNode(node);
        final double indicator = baselineIndicators[parameterIndex];
        if (Math.abs(indicator - 1.0) <= 1.0e-9) {
            final int atomState = (int) Math.round(baselineAtomIndices[parameterIndex]);
            return rewardBranchRates.getRawRewardForAtomState(atomState);
        }
        return baselineCtsRewards[parameterIndex];
    }

    private void snapshotRewardState(final double[] indicators,
                                     final double[] atomIndices,
                                     final double[] ctsRewards) {
        final Parameter indicator = rewardBranchRates.getIndicator();
        final Parameter atoms = rewardBranchRates.getAtomIndices();
        final Parameter cts = rewardBranchRates.getRateParameter();
        for (int i = 0; i < indicators.length; i++) {
            indicators[i] = indicator.getParameterValue(i);
            atomIndices[i] = atoms.getParameterValue(i);
            ctsRewards[i] = cts.getParameterValue(i);
        }
    }

    private void restoreRewardState(final double[] indicators,
                                    final double[] atomIndices,
                                    final double[] ctsRewards) {
        final Parameter indicator = rewardBranchRates.getIndicator();
        final Parameter atoms = rewardBranchRates.getAtomIndices();
        final Parameter cts = rewardBranchRates.getRateParameter();
        for (int i = 0; i < indicators.length; i++) {
            cts.setParameterValue(i, ctsRewards[i]);
            atoms.setParameterValue(i, atomIndices[i]);
            indicator.setParameterValue(i, indicators[i]);
        }
    }

    private int matchingAtomState(final double rawReward) {
        final int stateCount = rewardBranchRates.getRewardRatesMapping().getDimension();
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

    private static final class FullParameterSnapshot {
        final double[] indicators;
        final double[] atomIndices;
        final double[] ctsRewards;

        private FullParameterSnapshot(final RewardsAwareMixtureBranchRates branchRates) {
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

        private void restore(final RewardsAwareMixtureBranchRates branchRates) {
            final Parameter indicator = branchRates.getIndicator();
            final Parameter atoms = branchRates.getAtomIndices();
            final Parameter cts = branchRates.getRateParameter();
            for (int i = 0; i < indicators.length; i++) {
                cts.setParameterValue(i, ctsRewards[i]);
                atoms.setParameterValue(i, atomIndices[i]);
                indicator.setParameterValue(i, indicators[i]);
            }
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

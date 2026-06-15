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

import java.util.List;

/**
 * BEAGLE-backed CTMC edge-evidence provider for dependent reward processes.
 *
 * The preorder buffers are generated through BEAGLE's preorder simulation path.
 * Candidate kernels are evaluated locally from the substitution-model spectral
 * representation so all continuous and atomic candidates for a branch are
 * evaluated against the same cached messages.
 */
public final class BeagleRewardDependentCtmcEdgeEvidenceProvider
        implements RewardDependentCtmcEdgeEvidenceProvider {

    private final TreeDataLikelihood treeDataLikelihood;
    private final BeagleDataLikelihoodDelegate likelihoodDelegate;
    private final Tree tree;
    private final BranchModel branchModel;
    private final RewardsAwareMixtureBranchRates rewardBranchRates;
    private final SiteRateModel siteRateModel;
    private final BeaglePreOrderDelegate preOrderDelegate;
    private final ProcessSimulation preOrderSimulation;

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
    private final double[][] postPartialsByNode;

    private EigenDecomposition[] eigenDecompositions;
    private boolean[] allRealEigen;
    private ComplexBlockKernelUtils.ComplexKernelPlan[] complexPlans;

    public BeagleRewardDependentCtmcEdgeEvidenceProvider(final TreeDataLikelihood treeDataLikelihood) {
        if (treeDataLikelihood == null) {
            throw new IllegalArgumentException("treeDataLikelihood must be non-null");
        }

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

        this.preOrderDelegate = new BeaglePreOrderDelegate(
                treeDataLikelihood.getId() == null ? "dependentRewardEvidence" : treeDataLikelihood.getId(),
                tree,
                likelihoodDelegate
        );
        this.preOrderSimulation = new ProcessSimulation(treeDataLikelihood, preOrderDelegate);

        refreshSpectralStructures();
    }

    @Override
    public void prepare() {
        refreshSpectralStructures();
        preOrderSimulation.cacheSimulatedTraits(null);
        fillPostPartialsForAllNodes();
        fillTopPartialsFromRoot();
    }

    @Override
    public double logEvidence(final int branchNodeNumber, final double rawReward) {
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

        System.arraycopy(topPartialsByNode[branchNodeNumber], 0, prePartials, 0, flattenedLength);
        System.arraycopy(postPartialsByNode[branchNodeNumber], 0, postPartials, 0, flattenedLength);

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
            java.util.Arrays.fill(topPartialsByNode[i], 0.0);
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

    private static final class BeaglePreOrderDelegate extends AbstractBeagleGradientDelegate {

        private final String traitName;

        private BeaglePreOrderDelegate(final String name,
                                       final Tree tree,
                                       final BeagleDataLikelihoodDelegate likelihoodDelegate) {
            super(name + ".dependentRewardEvidencePreOrder", tree, likelihoodDelegate);
            this.traitName = name + ".dependentRewardEvidencePreOrder.preOrderTouch";
        }

        @Override
        protected DiscretePartialsType getPreOrderType() {
            return DiscretePartialsType.BOTTOM;
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
}

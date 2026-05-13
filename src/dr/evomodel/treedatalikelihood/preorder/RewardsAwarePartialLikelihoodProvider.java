package dr.evomodel.treedatalikelihood.preorder;

import beagle.Beagle;
import dr.evolution.alignment.PatternList;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.branchmodel.RewardsAwareBranchModel;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.EvolutionaryProcessDelegate;
import dr.evomodel.treedatalikelihood.ProcessSimulation;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.Model;

import java.util.Arrays;
import java.util.List;

/**
 * PartialLikelihoodProvider
 *
 * Lightweight provider for branch-wise "below" and "above" partial likelihood vectors.
 *
 * Responsibilities:
 *  - update BEAGLE preorder partials
 *  - return postorder ("below") partials for a node/branch
 *  - construct above-branch partials for a node/branch
 *
 *
 * Assumptions:
 *  - usePreOrder = true in the BeagleDataLikelihoodDelegate
 *  - typically single-pattern discrete data when used exactly as below-partials for tips,
 *    though internal-node partial retrieval is generic
 *  - current implementation extracts only the first pattern / first rate category state vector
 *    when reading partials from BEAGLE
 */

/*
* @author Filippo Monti
 */
public final class RewardsAwarePartialLikelihoodProvider {

    private final TreeDataLikelihood treeDataLikelihood;
    private final BeagleDataLikelihoodDelegate likelihoodDelegate;
    private final RewardsAwareBranchModel rewardsAwareBranchModel;

    private final Tree tree;
    private final Beagle beagle;
    private final PatternList patternList;
    private final EvolutionaryProcessDelegate evolutionaryProcessDelegate;
    private final SiteRateModel siteRateModel;

    private final int stateCount;
    private final int patternCount;
    private final int categoryCount;
    private final int preOrderPartialOffset;

    private final ProcessSimulation processSimulation;
    private final int[] nodeToTaxonIndex;

    /**
     * Lightweight delegate used only to trigger BEAGLE preorder-partial updates.
     */
    private final class PreOrderDelegate extends ProcessSimulationDelegate.AbstractDelegate {

        private boolean substitutionProcessKnown = false;
        private double[] scratch;


        private PreOrderDelegate(final String name, final Tree tree) {
            super(name, tree);
            this.scratch = new double[stateCount * patternCount * categoryCount];
            likelihoodDelegate.addModelListener(this);
            likelihoodDelegate.addModelRestoreListener(this);
        }

        @Override
        protected void constructTraits(Helper treeTraitHelper) {

        }

        @Override
        public void simulate(final int[] operations,
                             final int operationCount,
                             final int rootNodeNumber) {

            ensureSubstitutionProcessUpToDate();
            simulateRoot(rootNodeNumber);
            beagle.updatePrePartials(operations, operationCount, Beagle.NONE);
        }

        @Override
        public void setupStatistics() {
            // no-op
        }

        @Override
        protected void simulateRoot(final int rootNumber) {
            final double[] frequencies = evolutionaryProcessDelegate.getRootStateFrequencies();

            final double[] rootPreOrderPartial = scratch;
            for (int i = 0; i < patternCount * categoryCount; i++) {
                System.arraycopy(frequencies, 0, rootPreOrderPartial, i * stateCount, stateCount);
            }

            beagle.setPartials(getPreOrderPartialIndex(rootNumber), rootPreOrderPartial);
        }

        @Override
        protected void simulateNode(int v0, int v1, int v2, int v3, int v4) {
            throw new UnsupportedOperationException("Not used with BEAGLE preorder updates");
        }

        @Override
        public int vectorizeNodeOperations(final List<NodeOperation> nodeOperations, int rootNodeNumber, final int[] operations) {
            int k = 0;

            for (NodeOperation op : nodeOperations) {
                // For preorder:
                // destination = pre-partial at left child
                // source above  = pre-partial at parent
                // source below  = post-partial at sibling
                operations[k++] = getPreOrderPartialIndex(op.getLeftChild());
                operations[k++] = Beagle.NONE;
                operations[k++] = Beagle.NONE;
                operations[k++] = getPreOrderPartialIndex(op.getNodeNumber());
                operations[k++] = evolutionaryProcessDelegate.getMatrixIndex(op.getLeftChild());
                operations[k++] = getPostOrderPartialIndex(op.getRightChild());
                operations[k++] = evolutionaryProcessDelegate.getMatrixIndex(op.getRightChild());
            }

            return nodeOperations.size();
        }

        @Override
        public int getSingleOperationSize() {
            return Beagle.OPERATION_TUPLE_SIZE;
        }

        @Override
        public void modelChangedEvent(Model model, Object object, int index) {
            substitutionProcessKnown = false;
        }

        @Override
        public void modelRestored(Model model) {
            substitutionProcessKnown = false;
        }

        private void ensureSubstitutionProcessUpToDate() {
            // This hook is left here for consistency
            if (!substitutionProcessKnown) {
                substitutionProcessKnown = true;
            }
        }
    }

    private double[] scratch1;
    private double[] scratch2;
    private double[] scratch3;

    public RewardsAwarePartialLikelihoodProvider(final TreeDataLikelihood treeDataLikelihood,
                                     final BeagleDataLikelihoodDelegate likelihoodDelegate,
                                     final RewardsAwareBranchModel rewardsAwareBranchModel) {

        this.treeDataLikelihood = treeDataLikelihood;
        this.likelihoodDelegate = likelihoodDelegate;
        this.rewardsAwareBranchModel = rewardsAwareBranchModel;

        this.tree = rewardsAwareBranchModel.getTree();
        this.beagle = likelihoodDelegate.getBeagleInstance();
        this.patternList = likelihoodDelegate.getPatternList();
        this.evolutionaryProcessDelegate = likelihoodDelegate.getEvolutionaryProcessDelegate();
        this.siteRateModel = likelihoodDelegate.getSiteRateModel();

        if (!likelihoodDelegate.isUsePreOrder()) {
            throw new IllegalArgumentException("PartialLikelihoodProvider requires usePreOrder = true");
        }

        this.stateCount = patternList.getDataType().getStateCount();
        this.patternCount = patternList.getPatternCount();
        this.categoryCount = siteRateModel.getCategoryCount();

        // pre-order partials are stored immediately after the post-order partial buffers
        this.preOrderPartialOffset = likelihoodDelegate.getPartialBufferCount();

        final ProcessSimulationDelegate delegate = new PreOrderDelegate("partialLikelihoodProvider", tree);
        this.processSimulation = new ProcessSimulation(treeDataLikelihood, delegate);

        final int nodeCount = tree.getNodeCount();
        this.nodeToTaxonIndex = new int[nodeCount];

        this.scratch1 = new double[stateCount];
        this.scratch2 = new double[stateCount];
        this.scratch3 = new double[stateCount * patternCount * categoryCount];

        Arrays.fill(nodeToTaxonIndex, -1);

        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            NodeRef tip = tree.getExternalNode(i);
            Taxon taxon = tree.getNodeTaxon(tip);
            int taxonIndex = patternList.getTaxonIndex(taxon);

            if (taxonIndex < 0) {
                throw new IllegalStateException("Taxon not found in pattern list: " + taxon);
            }

            nodeToTaxonIndex[tip.getNumber()] = taxonIndex;
        }
    }

    /**
     * Ensure preorder partials are up to date for the current tree/likelihood state.
     */
    public void updatePreOrderPartials() {
        processSimulation.cacheSimulatedTraits(tree.getRoot());
    }

    /**
     * Get the postorder ("below") partial vector at the child node of the branch.
     *
     * For tips, this returns a one-hot vector for the observed state.
     * For internal nodes, this retrieves the node partial from BEAGLE and extracts
     * the first pattern / first category slice.
     */
    public void getBelowPartialsForBranch(final NodeRef node, final double[] postPartial) {
        checkNodeNotRoot(node);
        checkVectorLength(postPartial, stateCount, "postPartial");

        if (tree.isExternal(node)) {
            final int taxonIndex = nodeToTaxonIndex[node.getNumber()];
            final int state = patternList.getPattern(0)[taxonIndex];
            Arrays.fill(postPartial, 0.0);
            postPartial[state] = 1.0;
            return;
        }

        final int postIdx = getPostOrderPartialIndex(node.getNumber());
        copyFirstPatternFirstCategoryPartials(postIdx, postPartial);
    }

    /**
     * Get the "above" partial vector for the branch leading into {@code node}.
     *
     * This is the contribution from everything in the tree except the subtree below {@code node},
     * evaluated at the parent-side state space of that branch.
     *
     * Construction:
     *  - start from preorder partial at the parent
     *  - multiply in each sibling contribution pushed up to the parent
     */
    public void getAbovePartialsForBranch(final NodeRef node, final double[] prePartial) {
        checkNodeNotRoot(node);
        checkVectorLength(prePartial, stateCount, "prePartial");

        final NodeRef parent = tree.getParent(node);
        if (parent == null) {
            throw new IllegalStateException("Non-root node has null parent: " + node.getNumber());
        }

        // Start from preorder partial at the parent
        final int preIdxParent = getPreOrderPartialIndex(parent.getNumber());
        copyFirstPatternFirstCategoryPartials(preIdxParent, prePartial);

        final double[] siblingBelow = scratch1;
        final double[] siblingUp = scratch2;

        final int childCount = tree.getChildCount(parent);
        for (int k = 0; k < childCount; k++) {
            final NodeRef sibling = tree.getChild(parent, k);
            if (sibling.getNumber() == node.getNumber()) {
                continue;
            }

            getBelowPartialsForBranch(sibling, siblingBelow);

            final double[] P = rewardsAwareBranchModel.getTransitionMatrix(sibling.getNumber());
            mulMatVecParentToChild(P, stateCount, siblingBelow, siblingUp);

            for (int i = 0; i < stateCount; i++) {
                prePartial[i] *= siblingUp[i];
            }
        }
    }
    public void getAbovePartialsForNode(final NodeRef node, final double[] prePartial) {
        checkNodeNotRoot(node);
        checkVectorLength(prePartial, stateCount, "prePartial");
        copyFirstPatternFirstCategoryPartials(getPreOrderPartialIndex(node.getNumber()), prePartial);
    }

    public void getBelowPartialsForBranch(final int branchIndex, final double[] postPartial) {
        final NodeRef node = getNodeForBranch(branchIndex);
        getBelowPartialsForBranch(node, postPartial);
    }

    public void getAbovePartialsForBranch(final int branchIndex, final double[] prePartial) {
        final NodeRef node = getNodeForBranch(branchIndex);
        getAbovePartialsForBranch(node, prePartial);
    }

    public int getStateCount() {
        return stateCount;
    }

    public Tree getTree() {
        return tree;
    }

    public TreeDataLikelihood getTreeDataLikelihood() {
        return treeDataLikelihood;
    }

    public BeagleDataLikelihoodDelegate getLikelihoodDelegate() {
        return likelihoodDelegate;
    }

    public RewardsAwareBranchModel getRewardsAwareBranchModel() {
        return rewardsAwareBranchModel;
    }

    protected int getPostOrderPartialIndex(final int nodeNumber) {
        return likelihoodDelegate.getPartialBufferIndex(nodeNumber);
    }

    protected int getPreOrderPartialIndex(final int nodeNumber) {
        return preOrderPartialOffset + nodeNumber;
    }

    private NodeRef getNodeForBranch(final int branchIndex) {
        final int nodeNum = rewardsAwareBranchModel.getNodeNumberForBranchIndex(branchIndex);
        final NodeRef node = tree.getNode(nodeNum);
        checkNodeNotRoot(node);
        return node;
    }

    private void checkNodeNotRoot(final NodeRef node) {
        if (tree.isRoot(node)) {
            throw new IllegalStateException(
                    "Requested branch corresponds to root node " + node.getNumber() + ", which is not allowed."
            );
        }
    }

    private static void checkVectorLength(final double[] x, final int expected, final String name) {
        if (x.length != expected) {
            throw new IllegalArgumentException(
                    name + " must have length " + expected + " but got " + x.length
            );
        }
    }

    /**
     * Extract the first pattern / first category state vector from a BEAGLE partial buffer.
     */
    private void copyFirstPatternFirstCategoryPartials(final int partialIndex, final double[] out) {
        final double[] all = scratch3;
        beagle.getPartials(partialIndex, Beagle.NONE, all);
        System.arraycopy(all, 0, out, 0, stateCount);
    }

    /**
     * out[i] = sum_j P[i,j] * v[j]
     * where P is row-major with rows indexed by parent state and columns by child state.
     */
    private static void mulMatVecParentToChild(final double[] P,
                                               final int dim,
                                               final double[] v,
                                               final double[] out) {
        Arrays.fill(out, 0.0);

        for (int i = 0; i < dim; i++) {
            final int rowBase = i * dim;
            double acc = 0.0;
            for (int j = 0; j < dim; j++) {
                acc += P[rowBase + j] * v[j];
            }
            out[i] = acc;
        }
    }
}
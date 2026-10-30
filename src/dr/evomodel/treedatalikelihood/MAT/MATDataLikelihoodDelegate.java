/*
 * MATDataLikelihoodDelegate.java
 *
 * Copyright © 2026, the BEAST Development Team.
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
 */

package dr.evomodel.treedatalikelihood.MAT;

import dr.evolution.alignment.PatternList;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.siteratemodel.SiteRateModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeChangedEvent;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.PreOrderSettings;
import dr.evomodel.treedatalikelihood.RateRescalingScheme;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.TreeTraversal;
import dr.evomodel.treelikelihood.PartialsRescalingScheme;
import dr.inference.model.AbstractModel;
import dr.inference.model.HastingsRatioContributor;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;

/**
 * Mutation-annotated tree (MAT) likelihood delegate, following the model of Varilly et al. (Nature,
 * 2026).
 *
 * Rather than integrating over all ancestral states, the likelihood is that of an explicit
 * mutational history: a root sequence plus a set of mutations on each branch, reconstructed by
 * parsimony ({@link MATInitializer}). The ancestral states can then be sampled jointly with the
 * topology by SPR operators carrying a {@link MutationHistoryResampler}. If the topology is changed by
 * any other operator, the MAT is brought up to date by (incremental) parsimony, which makes the
 * likelihood a deterministic function of the tree. The two should not be mixed in one analysis.
 *
 * How each site's history along a branch is represented and scored is delegated to a
 * {@link BranchHistoryModel}: by default {@link SingleHitHistoryModel} (at most one substitution per
 * site per branch, its time integrated out); alternatives allow multiple hits, summed analytically or
 * as explicit (optionally timed) paths. The root sequence has prior sum_l log pi_a.
 *
 * The topology is changed by the tree operators alone: when the likelihood is next evaluated the delegate
 * finds the subtree prune-and-regraft move from the nodes whose parent changed, has the resampler draw the
 * states the move changes (scoring the old ones under the old topology, which it keeps as the store's
 * parents and the heights of the last accepted state), and reports the log Hastings ratio of that
 * proposal to the Markov chain ({@link HastingsRatioContributor}). Without a resampler the MAT is brought up
 * to date by parsimony instead.
 *
 * The likelihood is a root term plus a term for each branch, each cached against what it depends on.
 * Only the terms of branches whose states, events or length may have changed are recomputed (the store
 * logs the first two, events from the tree the third), so an evaluation costs O(changed branches) plus a
 * sum over the branch terms.
 */
public class MATDataLikelihoodDelegate extends AbstractModel implements DataLikelihoodDelegate, HastingsRatioContributor {

    public static final String NAME = "MATDataLikelihoodDelegate";

    private final PatternList patternList;
    private final Tree tree;
    private final BranchModel branchModel;
    private final SubstitutionModel substitutionModel;
    private final FrequencyModel rootFrequencyModel;
    private final SiteRateModel siteRateModel;
    private final BranchRateModel branchRateModel;
    private final Parameter siteRates;

    private final int patternCount;
    private final int stateCount;
    private final double[] patternWeights;
    private final boolean[] constantPatterns;

    private final BranchHistoryModel historyModel;
    private final MATInitializer initializer;
    private final DefaultMutationStore mutationStore;
    private int rebuildCount = 0;
    private boolean historyResampled = false;
    private boolean mixingWarningGiven = false;

    private final double[] siteLogLikelihoods;
    private boolean siteLogLikelihoodsKnown = false;
    private double logLikelihood;
    private double storedLogLikelihood;

    // Incremental likelihood: logL = root term + sum over branches of a branch term, each cached
    // against what it depends on. modelVersion changes whenever the rate matrix, frequencies or
    // site rates may have changed.
    private long modelVersion = 0;
    private boolean modelChangedSinceStore = false;
    private long calculationCount = 0;

    // the rates, for the model version they were last read at (see ensureRates)
    private long cachedModelVersion = -1;
    private final double[] infinitesimalMatrix;
    private final double[] newInfinitesimalMatrix;
    // changes only when the rate matrix or the category rates change (not, e.g., the categories themselves)
    private long ratesVersion = 0;
    private boolean ratesKnown = false;
    private final double[] exitRates;           // q_a = sum over b != a of Q_ab
    private final double[] patternRates;        // nu_p
    private final double[] patternWeightRates;  // w_p nu_p
    private double commonSiteRate;              // the rate of every site, without per-site rates

    // the root term, with the root node, its state version and the model version it was computed with
    private double rootTerm;
    private int rootTermNode = -1;
    private long rootTermStateVersion = -1;
    private long rootTermModelVersion = -1;

    // the term of the branch above each node, with the length, state version, events and model version it
    // was computed with
    private final double[] branchTerms;
    private final double[] branchTermLengths;
    private final long[] branchTermStateVersions;
    private final BranchMutations[] branchTermMutations;
    private final long[] branchTermModelVersions;

    // each node's exit rate lambda_n = sum_p w_p nu_p q(state of p at n), with the state and model versions
    // it was computed with
    private final double[] exitRateSums;
    private final long[] exitRateStateVersions;
    private final long[] exitRateModelVersions;
    private final int[] exitRatePath;           // scratch: nodes whose exit rates are being brought up to date

    // The branch terms are only recomputed for nodes whose states, branch events or branch length may
    // have changed: the store logs the first two, events from the tree the third.
    private final IndexSet dirtyNodes;
    private boolean allNodesDirty = true;
    private long termsModelVersion = -1;
    // nodes whose terms were recomputed since storeState (the tree restores their lengths silently)
    private final IndexSet recomputedNodes;
    // whether the topology may have changed since the MAT was last checked against it
    private boolean treeChangedSinceCheck = true;

    // The subtree moves made by the operators are found and resampled at evaluation: the nodes with a tree
    // event since the step began (or all of them), the heights of the last accepted state (for the branch
    // lengths of the old topology), whether a step is in progress (between store and accept/restore), the
    // moves found in it, and the log Hastings ratio of their resampling, taken by the Markov chain.
    private MutationHistoryResampler resampler = null;
    private final IndexSet touchedNodes;
    private boolean allTouched = true;
    private final double[] knownHeights;
    private boolean inStep = false;
    private int movesThisStep = 0;
    private double pendingLogHastingsRatio = 0.0;
    private final TopologyView currentView;
    private int[] lastSubtreeMove = null;

    // Discrete site rate categories (a site rate model with more than one category), or null: each
    // variable pattern's category is sampled, and the constant patterns are summed over the categories.
    private final SiteRateCategories categories;

    public MATDataLikelihoodDelegate(Tree tree,
                                     PatternList patternList,
                                     BranchModel branchModel,
                                     SiteRateModel siteRateModel,
                                     BranchRateModel branchRateModel,
                                     Parameter siteRates) {
        this(tree, patternList, branchModel, siteRateModel, branchRateModel, siteRates, null);
    }

    /**
     * @param historyModel how histories along branches are represented (null for single hits)
     */
    public MATDataLikelihoodDelegate(Tree tree,
                                     PatternList patternList,
                                     BranchModel branchModel,
                                     SiteRateModel siteRateModel,
                                     BranchRateModel branchRateModel,
                                     Parameter siteRates,
                                     BranchHistoryModel historyModel) {
        super(NAME);
        this.tree = tree;
        this.patternList = patternList;
        this.branchModel = branchModel;
        this.siteRateModel = siteRateModel;
        this.branchRateModel = branchRateModel == null ? new DefaultBranchRateModel() : branchRateModel;
        this.siteRates = siteRates;

        if (branchModel.getSubstitutionModels().size() != 1) {
            throw new IllegalArgumentException("MAT likelihood currently requires a single substitution model across the tree");
        }
        this.substitutionModel = branchModel.getSubstitutionModels().get(0);
        this.rootFrequencyModel = branchModel.getRootFrequencyModel();

        final int categoryCount = siteRateModel == null ? 1 : siteRateModel.getCategoryCount();
        if (categoryCount > 1 && siteRates != null) {
            throw new IllegalArgumentException("MAT likelihood takes either a per-site rate parameter or a site rate " +
                    "model with more than one category, not both");
        }
        if (categoryCount > Byte.MAX_VALUE) {
            throw new IllegalArgumentException("MAT likelihood supports at most " + Byte.MAX_VALUE + " site rate categories");
        }

        this.patternCount = patternList.getPatternCount();
        this.stateCount = patternList.getDataType().getStateCount();
        if (substitutionModel.getDataType().getStateCount() != stateCount) {
            throw new IllegalArgumentException("Substitution model and pattern list have different state counts");
        }
        if (siteRates != null && siteRates.getDimension() != patternCount) {
            throw new IllegalArgumentException("Site rate parameter " + siteRates.getId() + " has dimension " +
                    siteRates.getDimension() + " but there are " + patternCount + " patterns");
        }

        this.patternWeights = patternList.getPatternWeights();

        // Each site has its own mutational history, so the patterns must not be compressed, except that
        // constant sites (which are assumed to have no mutations) may be pooled by state, as produced by
        // the constantOnly compression of SitePatterns.
        this.constantPatterns = new boolean[patternCount];
        int constantSites = 0;
        int variableSites = 0;
        for (int p = 0; p < patternCount; p++) {
            int first = patternList.getPatternState(0, p);
            boolean constant = first < stateCount;
            for (int i = 1; constant && i < patternList.getPatternLength(); i++) {
                constant = patternList.getPatternState(i, p) == first;
            }
            constantPatterns[p] = constant;
            if (constant) {
                constantSites += (int) patternWeights[p];
            } else {
                if (patternWeights[p] != 1.0) {
                    throw new IllegalArgumentException("MAT likelihood needs each variable site as its own pattern but pattern " +
                            (p + 1) + " has weight " + patternWeights[p] + ": use uncompressed patterns or constant-only " +
                            "compression (unique=\"false\" or constantOnly=\"true\" on the patterns element)");
                }
                variableSites++;
            }
        }
        Logger.getLogger("dr.evomodel").info("  MAT likelihood: " + variableSites + " variable sites, " + constantSites +
                " constant sites (assumed to have no mutations)");
        this.siteLogLikelihoods = new double[patternCount];
        int nodeCount = tree.getNodeCount();
        this.patternRates = new double[patternCount];
        this.patternWeightRates = new double[patternCount];
        this.exitRateSums = new double[nodeCount];
        this.exitRateStateVersions = new long[nodeCount];
        this.exitRateModelVersions = new long[nodeCount];
        this.exitRatePath = new int[nodeCount];
        this.branchTerms = new double[nodeCount];
        this.branchTermLengths = new double[nodeCount];
        this.branchTermStateVersions = new long[nodeCount];
        this.branchTermMutations = new BranchMutations[nodeCount];
        this.branchTermModelVersions = new long[nodeCount];
        Arrays.fill(exitRateModelVersions, -1);
        Arrays.fill(branchTermModelVersions, -1);
        this.dirtyNodes = new IndexSet(nodeCount);
        this.recomputedNodes = new IndexSet(nodeCount);
        this.touchedNodes = new IndexSet(nodeCount);
        this.knownHeights = new double[nodeCount];
        for (int n = 0; n < nodeCount; n++) {
            knownHeights[n] = tree.getNodeHeight(tree.getNode(n));
        }
        this.currentView = new TreeView();
        this.infinitesimalMatrix = new double[stateCount * stateCount];
        this.newInfinitesimalMatrix = new double[stateCount * stateCount];
        this.exitRates = new double[stateCount];

        this.historyModel = historyModel == null ? new SingleHitHistoryModel(stateCount) : historyModel;
        this.historyModel.checkCompatibility(this);
        Logger.getLogger("dr.evomodel").info("  MAT branch history model: " + this.historyModel.getName());

        initializer = new MATInitializer(patternList, tree);
        mutationStore = initializer.buildMutationStore();
        this.historyModel.initialise(mutationStore);

        this.categories = categoryCount > 1 ? new SiteRateCategories(this, siteRateModel, patternList) : null;

        addModel(branchModel);
        if (siteRateModel != null) {
            addModel(siteRateModel);
        }
        addModel(this.branchRateModel);
        if (siteRates != null) {
            addVariable(siteRates);
        }
        if (tree instanceof Model) {
            // a listener only (not a sub-model): to learn which branch lengths, and when the topology, change
            ((Model) tree).addModelListener(this);
        }
    }

    public MutationStore getMutationStore() {
        return mutationStore;
    }

    public BranchHistoryModel getHistoryModel() {
        return historyModel;
    }

    DefaultMutationStore getMutableMutationStore() {
        return mutationStore;
    }

    /**
     * Attaches the resampler that draws the ancestral states changed by a subtree move when the likelihood
     * is next evaluated (a {@link MutationHistoryResampler} attaches itself when constructed; null detaches
     * it, after which the MAT is updated by parsimony, which is not a valid MCMC move when histories are
     * sampled).
     */
    public void setResampler(MutationHistoryResampler resampler) {
        if (resampler != null && !(branchRateModel instanceof StrictClockBranchRates || branchRateModel instanceof DefaultBranchRateModel)) {
            throw new IllegalArgumentException("The MAT's mutation history resampler needs a strict clock (the rates of " +
                    "the branches of the topology before a move cannot be read from " + branchRateModel.getModelName() + ")");
        }
        this.resampler = resampler;
    }

    public MutationHistoryResampler getResampler() {
        return resampler;
    }

    /**
     * @return the last subtree move the resampler was given, as {X, S, D}: the subtree below X moved, with
     * its parent, from beside S to the branch above D (for testing)
     */
    public int[] getLastSubtreeMove() {
        return lastSubtreeMove == null ? null : lastSubtreeMove.clone();
    }

    @Override
    public double takeLogHastingsRatio() {
        final double ratio = pendingLogHastingsRatio;
        pendingLogHastingsRatio = 0.0;
        return ratio;
    }

    /**
     * @return the tree as it is
     */
    TopologyView currentView() {
        return currentView;
    }

    /**
     * Called by a {@link MutationHistoryResampler} when it has sampled new ancestral states.
     */
    void setHistoryResampled() {
        historyResampled = true;
    }

    /**
     * Called by a {@link MutationHistoryResampler} when it has changed the history without changing the tree,
     * so that the tree likelihood (which otherwise hears only of tree changes) recalculates. The branches to
     * update are found from the mutation store's change log.
     */
    void fireHistoryChanged() {
        // index -1: TreeDataLikelihood would otherwise read the index as a node number
        fireModelChanged(mutationStore, -1);
    }

    /**
     * @return the number of times the MAT has been updated by parsimony after a topology change
     */
    public int getRebuildCount() {
        return rebuildCount;
    }

    public Tree getTree() {
        return tree;
    }

    public PatternList getPatternList() {
        return patternList;
    }

    public BranchModel getBranchModel() {
        return branchModel;
    }

    public SiteRateModel getSiteRateModel() {
        return siteRateModel;
    }

    public BranchRateModel getBranchRateModel() {
        return branchRateModel;
    }

    /**
     * Makes the MAT consistent with the current topology, by parsimony if necessary.
     */
    void updateMAT() {
        final long profileStart = MATProfile.start();
        try {
            updateMATImpl();
        } finally {
            MATProfile.stop(MATProfile.Phase.LIK_UPDATE_MAT, profileStart);
        }
    }

    private void updateMATImpl() {
        if (!treeChangedSinceCheck) {
            return;
        }
        // (cleared first: the resampler is called from here)
        treeChangedSinceCheck = false;
        final int[] reparented = findReparentedNodes();
        if (reparented.length == 0) {
            assert mutationStore.hasParents(MATInitializer.getParents(tree)) : "the MAT does not match the topology";
            return;
        }
        if (resampler != null && inStep) {
            if (++movesThisStep > 1) {
                throw new IllegalStateException("The topology of " + tree.getId() + " was changed more than once in " +
                        "one MCMC step: the MAT's mutation history resampler can only resample one subtree move per step");
            }
            final int[] move = identifySubtreeMove(reparented);
            if (move == null) {
                throw new IllegalStateException("The topology of " + tree.getId() + " was changed by something other " +
                        "than a single subtree prune-and-regraft move (" + reparented.length + " nodes were reparented): " +
                        "the MAT's mutation history resampler does not support this operator");
            }
            final TopologyView oldView = new OldView(reparented);
            lastSubtreeMove = move;
            pendingLogHastingsRatio += resampler.resampleMove(move[0], move[1], move[2], oldView, currentView);
            assert mutationStore.hasParents(MATInitializer.getParents(tree)) : "the MAT does not match the topology after resampling";
        } else {
            initializer.update(mutationStore);
            historyModel.initialise(mutationStore);
            rebuildCount++;
            if (inStep && historyResampled && !mixingWarningGiven) {
                // (at the info level, as BEAST prints warnings only when verbose)
                Logger.getLogger("dr.evomodel").info("WARNING: the topology of " + tree.getId() +
                        " was changed by an operator without a mutation history resampler: the MAT has been " +
                        "updated by parsimony, which is not a valid MCMC move when histories are sampled.");
                mixingWarningGiven = true;
            }
        }
    }

    /**
     * @return the nodes whose parent in the tree differs from their parent in the store (the topology of the
     * last accepted state), found among the nodes with a tree event since the step began and their children
     */
    private int[] findReparentedNodes() {
        final IntBuffer reparented = new IntBuffer();
        if (allTouched) {
            for (int n = 0; n < tree.getNodeCount(); n++) {
                if (treeParent(n) != mutationStore.getParent(n)) {
                    reparented.add(n);
                }
            }
        } else {
            final int[] seen = new int[tree.getNodeCount()];
            for (int i = 0; i < touchedNodes.size(); i++) {
                final int t = touchedNodes.get(i);
                final NodeRef node = tree.getNode(t);
                for (int c = -1; c < tree.getChildCount(node); c++) {
                    final int n = c < 0 ? t : tree.getChild(node, c).getNumber();
                    if (seen[n] == 0) {
                        seen[n] = 1;
                        if (treeParent(n) != mutationStore.getParent(n)) {
                            reparented.add(n);
                        }
                    }
                }
            }
        }
        return reparented.toArray();
    }

    private int treeParent(int n) {
        final NodeRef parent = tree.getParent(tree.getNode(n));
        return parent == null ? -1 : parent.getNumber();
    }

    /**
     * Reads a single subtree prune-and-regraft move from the reparented nodes: the subtree below X, with its
     * parent P, was pruned from beside its sibling S (under G, P's old parent) and regrafted on the branch
     * above D (under N, D's old parent). The reparented nodes are then S, D and, unless N is G, P.
     * <p>
     * The same change of topology can have two such descriptions (when N is S, or D is G): each is a
     * legal move of a different subtree, resampling different states. The description must be a function
     * of the pair of topologies alone, so that a move and its reverse resample the same states: the one
     * with the smallest moved parent P is chosen.
     *
     * @return {X, S, D}, or null if the change is not a single subtree move
     */
    private int[] identifySubtreeMove(int[] reparented) {
        int bestP = Integer.MAX_VALUE;
        int bestS = -1;
        int bestD = -1;
        for (int s : reparented) {
            final int p = mutationStore.getParent(s);
            if (p < 0 || treeParent(s) != mutationStore.getParent(p)) {
                continue;
            }
            final int g = mutationStore.getParent(p);
            final int n = treeParent(p);
            for (int d : reparented) {
                if (d == s || treeParent(d) != p || mutationStore.getParent(d) == p || n != mutationStore.getParent(d)) {
                    continue;
                }
                final int expected = g != n ? 3 : 2;
                if (reparented.length != expected) {
                    continue;
                }
                if (g != n && !contains(reparented, p)) {
                    continue;
                }
                if (p < bestP) {
                    bestP = p;
                    bestS = s;
                    bestD = d;
                }
            }
        }
        if (bestS < 0) {
            return null;
        }
        final NodeRef pNode = tree.getNode(bestP);
        int x = -1;
        for (int i = 0; i < tree.getChildCount(pNode); i++) {
            final int c = tree.getChild(pNode, i).getNumber();
            if (c != bestD) {
                x = c;
            }
        }
        return x < 0 ? null : new int[]{x, bestS, bestD};
    }

    private static boolean contains(int[] nodes, int n) {
        for (int m : nodes) {
            if (m == n) {
                return true;
            }
        }
        return false;
    }

    /**
     * The tree as it is.
     */
    private final class TreeView implements TopologyView {

        @Override
        public int getParent(int node) {
            return treeParent(node);
        }

        @Override
        public int getChildCount(int node) {
            return tree.getChildCount(tree.getNode(node));
        }

        @Override
        public int getChild(int node, int i) {
            return tree.getChild(tree.getNode(node), i).getNumber();
        }

        @Override
        public boolean isExternal(int node) {
            return tree.isExternal(tree.getNode(node));
        }

        @Override
        public int getRoot() {
            return tree.getRoot().getNumber();
        }

        @Override
        public double getBranchLength(int node) {
            return MATDataLikelihoodDelegate.this.getBranchLength(tree.getNode(node));
        }
    }

    /**
     * The tree as it was before the subtree move that reparented the given nodes: the parents are the
     * store's (not yet synchronised with the tree), the children of the few nodes whose children changed
     * follow from them, and the branch lengths come from the heights of the last accepted state.
     */
    private final class OldView implements TopologyView {

        private final int[] reparented;
        private final int[] oldParents;
        private final int[] newParents;
        private final int[] scratch;

        OldView(int[] reparented) {
            this.reparented = reparented;
            this.oldParents = new int[reparented.length];
            this.newParents = new int[reparented.length];
            for (int i = 0; i < reparented.length; i++) {
                oldParents[i] = mutationStore.getParent(reparented[i]);
                newParents[i] = treeParent(reparented[i]);
            }
            this.scratch = new int[2 + reparented.length];
        }

        @Override
        public int getParent(int node) {
            return mutationStore.getParent(node);
        }

        /**
         * @return the old children of the node in scratch, and their number
         */
        private int oldChildren(int node) {
            int count = 0;
            final NodeRef n = tree.getNode(node);
            for (int i = 0; i < tree.getChildCount(n); i++) {
                final int c = tree.getChild(n, i).getNumber();
                if (mutationStore.getParent(c) == node) {
                    scratch[count++] = c;
                }
            }
            for (int i = 0; i < reparented.length; i++) {
                if (oldParents[i] == node && newParents[i] != node) {
                    scratch[count++] = reparented[i];
                }
            }
            return count;
        }

        @Override
        public int getChildCount(int node) {
            return oldChildren(node);
        }

        @Override
        public int getChild(int node, int i) {
            oldChildren(node);
            return scratch[i];
        }

        @Override
        public boolean isExternal(int node) {
            return tree.isExternal(tree.getNode(node));
        }

        @Override
        public int getRoot() {
            for (int i = 0; i < reparented.length; i++) {
                if (oldParents[i] < 0) {
                    return reparented[i];
                }
            }
            return tree.getRoot().getNumber();
        }

        @Override
        public double getBranchLength(int node) {
            final int parent = mutationStore.getParent(node);
            if (parent < 0) {
                return 0.0;
            }
            return (knownHeights[parent] - knownHeights[node]) * branchRateModel.getBranchRate(tree, tree.getNode(node));
        }
    }

    /**
     * The MAT tracks changes to the tree itself (see updateMAT), so the tree likelihood need not traverse
     * the tree to collect operations.
     */
    @Override
    public boolean requiresTreeTraversal() {
        return false;
    }

    @Override
    public double calculateLikelihood(List<BranchOperation> branchOperations,
                                      List<NodeOperation> nodeOperations,
                                      int rootNodeNumber) throws LikelihoodException {
        final long profileStart = MATProfile.start();
        updateMAT();
        if (categories != null) {
            categories.warnIfNoOperator();
        }
        computeLikelihood();
        calculationCount++;
        MATProfile.stop(MATProfile.Phase.LIK_CALCULATE, profileStart);
        MATProfile.count(MATProfile.Counter.EVALUATIONS, 1);
        return logLikelihood;
    }

    /**
     * Brings the rate matrix, the site rates and the history model's rates up to date.
     */
    void ensureRates() {
        final long profileStart = MATProfile.start();
        try {
            ensureRatesImpl();
        } finally {
            MATProfile.stop(MATProfile.Phase.LIK_ENSURE_RATES, profileStart);
        }
    }

    private void ensureRatesImpl() {
        if (cachedModelVersion != modelVersion) {
            substitutionModel.getInfinitesimalMatrix(newInfinitesimalMatrix);
            // (the history model's caches stay valid if the rate matrix is unchanged)
            final boolean matrixChanged = !ratesKnown || !java.util.Arrays.equals(newInfinitesimalMatrix, infinitesimalMatrix);
            if (matrixChanged) {
                System.arraycopy(newInfinitesimalMatrix, 0, infinitesimalMatrix, 0, infinitesimalMatrix.length);
                for (int a = 0; a < stateCount; a++) {
                    double q = 0.0;
                    for (int b = 0; b < stateCount; b++) {
                        if (b != a) {
                            q += infinitesimalMatrix[a * stateCount + b];
                        }
                    }
                    exitRates[a] = q;
                }
                historyModel.setRates(substitutionModel, infinitesimalMatrix, exitRates);
                ratesVersion++;
                ratesKnown = true;
            }
            if (categories != null) {
                if (categories.updateRates()) {
                    ratesVersion++;
                }
                for (int p = 0; p < patternCount; p++) {
                    patternRates[p] = categories.getPatternRate(p);
                    patternWeightRates[p] = patternWeights[p] * patternRates[p];
                }
            } else {
                double siteRateScale = getSiteRateScale();
                commonSiteRate = siteRateScale;
                for (int p = 0; p < patternCount; p++) {
                    patternRates[p] = getSiteRate(p, siteRateScale);
                    patternWeightRates[p] = patternWeights[p] * patternRates[p];
                }
            }
            cachedModelVersion = modelVersion;
        }
    }

    private void computeLikelihood() {
        ensureRates();
        if (termsModelVersion != modelVersion) {
            allNodesDirty = true;
            termsModelVersion = modelVersion;
        }
        for (int i = 0; i < mutationStore.getChangedNodeCount(); i++) {
            dirtyNodes.add(mutationStore.getChangedNode(i));
        }
        mutationStore.clearChangeLog();

        final int rootNumber = tree.getRoot().getNumber();
        long profileStart = MATProfile.start();
        if (allNodesDirty) {
            if (categories != null) {
                categories.resetTotals();
            }
            for (int n = 0; n < tree.getNodeCount(); n++) {
                if (n != rootNumber) {
                    updateBranchTerm(n);
                }
            }
            allNodesDirty = false;
        } else {
            for (int i = 0; i < dirtyNodes.size(); i++) {
                int n = dirtyNodes.get(i);
                if (n != rootNumber) {
                    updateBranchTerm(n);
                }
            }
        }
        dirtyNodes.clear();
        MATProfile.stop(MATProfile.Phase.LIK_BRANCH_TERMS, profileStart);

        profileStart = MATProfile.start();
        double total = getRootTerm(rootNumber);
        for (int n = 0; n < branchTerms.length; n++) {
            if (n != rootNumber) {
                total += branchTerms[n];
            }
        }
        if (categories != null) {
            total += categories.getConstantTerm(rootNumber, ratesVersion) + categories.getPriorTerm(modelVersion);
        }
        MATProfile.stop(MATProfile.Phase.LIK_SUM, profileStart);

        logLikelihood = total;
        siteLogLikelihoodsKnown = false;
    }

    /**
     * Recomputes the term of the branch above node n if anything it depends on has changed.
     */
    private void updateBranchTerm(int n) {
        double exitRate = historyModel.needsExitRate() ? getExitRate(n) : 0.0;
        double length = getBranchLength(tree.getNode(n));
        BranchMutations mutations = mutationStore.getBranchMutations(n);
        long stateVersion = mutationStore.getStateVersion(n);
        if (categories != null) {
            categories.updateUnchangedTerm(n, length, ratesVersion);
        }
        if (branchTermModelVersions[n] != modelVersion || branchTermLengths[n] != length ||
                branchTermStateVersions[n] != stateVersion || branchTermMutations[n] != mutations) {
            if (categories != null) {
                double[] categoryComposition = historyModel.needsComposition() ?
                        mutationStore.getCategoryComposition(n) : null;
                branchTerms[n] = historyModel.getCategoryBranchLogLikelihood(mutations, length, exitRate,
                        categoryComposition, categories.getRates(), categories.getUnchangedTerms(n),
                        patternWeights, patternRates);
            } else {
                double[] composition = historyModel.needsComposition() ? mutationStore.getComposition(n) : null;
                branchTerms[n] = historyModel.getBranchLogLikelihood(mutations, length, exitRate, composition,
                        patternWeights, patternRates);
            }
            branchTermLengths[n] = length;
            branchTermStateVersions[n] = stateVersion;
            branchTermMutations[n] = mutations;
            branchTermModelVersions[n] = modelVersion;
            MATProfile.count(MATProfile.Counter.TERMS_RECOMPUTED, 1);
        }
        recomputedNodes.add(n);
    }

    /**
     * @return lambda_n, the total exit rate of the node's sequence: from the node's composition if all
     * sites have the same rate, otherwise from its parent's and the events on its branch (bringing the
     * rates of its ancestors up to date first, from the nearest one that is)
     */
    private double getExitRate(int n) {
        if (isExitRateKnown(n)) {
            return exitRateSums[n];
        }
        if (!hasPatternRates()) {
            double[] composition = mutationStore.getComposition(n);
            double rate = 0.0;
            for (int a = 0; a < stateCount; a++) {
                rate += composition[a] * exitRates[a];
            }
            setExitRate(n, rate * commonSiteRate);
            return exitRateSums[n];
        }
        int count = 0;
        int top = n;
        exitRatePath[count++] = top;
        while (mutationStore.getParent(top) >= 0 && !isExitRateKnown(mutationStore.getParent(top))) {
            top = mutationStore.getParent(top);
            exitRatePath[count++] = top;
        }
        for (int i = count - 1; i >= 0; i--) {
            int node = exitRatePath[i];
            int parent = mutationStore.getParent(node);
            double rate;
            if (parent < 0) {
                rate = computeExitRate(node);
            } else {
                rate = exitRateSums[parent];
                BranchMutations mutations = mutationStore.getBranchMutations(node);
                for (int k = 0; k < mutations.size(); k++) {
                    int p = mutations.getPattern(k);
                    int a = mutations.getFromState(k);
                    if (mutations.isAmbiguity(k)) {
                        rate -= patternWeightRates[p] * exitRates[a];
                    } else {
                        rate += patternWeightRates[p] * (exitRates[mutations.getToState(k)] - exitRates[a]);
                    }
                }
            }
            setExitRate(node, rate);
        }
        return exitRateSums[n];
    }

    private boolean isExitRateKnown(int n) {
        return exitRateModelVersions[n] == modelVersion && exitRateStateVersions[n] == mutationStore.getStateVersion(n);
    }

    private void setExitRate(int n, double rate) {
        exitRateSums[n] = rate;
        exitRateStateVersions[n] = mutationStore.getStateVersion(n);
        exitRateModelVersions[n] = modelVersion;
    }

    private double computeExitRate(int n) {
        byte[] states = mutationStore.getStates(n);
        double rate = 0.0;
        for (int p = 0; p < patternCount; p++) {
            if (states[p] >= 0) {
                rate += patternWeightRates[p] * exitRates[states[p]];
            }
        }
        return rate;
    }

    private double getRootTerm(int rootNumber) {
        long stateVersion = mutationStore.getStateVersion(rootNumber);
        if (rootTermNode != rootNumber || rootTermStateVersion != stateVersion || rootTermModelVersion != modelVersion) {
            double[] frequencies = rootFrequencyModel.getFrequencies();
            double[] composition = mutationStore.getComposition(rootNumber);
            double term = 0.0;
            for (int a = 0; a < stateCount; a++) {
                if (composition[a] > 0.0) {
                    term += composition[a] * Math.log(frequencies[a]);
                }
            }
            rootTerm = term;
            rootTermNode = rootNumber;
            rootTermStateVersion = stateVersion;
            rootTermModelVersion = modelVersion;
        }
        return rootTerm;
    }

    /**
     * Computes the per-pattern log likelihoods in full (only needed when they are asked for), site by
     * site and branch by branch with the history model.
     */
    private void computeSiteLogLikelihoods() {
        final long profileStart = MATProfile.start();
        try {
            computeSiteLogLikelihoodsImpl();
        } finally {
            MATProfile.stop(MATProfile.Phase.LIK_SITE_LOG_LIKELIHOODS, profileStart);
        }
    }

    private void computeSiteLogLikelihoodsImpl() {
        ensureRates();
        double[] frequencies = rootFrequencyModel.getFrequencies();
        int rootNumber = tree.getRoot().getNumber();

        byte[] rootStates = mutationStore.getStates(rootNumber);
        for (int p = 0; p < patternCount; p++) {
            int a = rootStates[p];
            siteLogLikelihoods[p] = a < 0 ? 0.0 : Math.log(frequencies[a]);
        }
        int[] eventStart = new int[patternCount];
        int[] eventEnd = new int[patternCount];
        for (int n = 0; n < tree.getNodeCount(); n++) {
            if (n == rootNumber) {
                continue;
            }
            double length = getBranchLength(tree.getNode(n));
            byte[] states = mutationStore.getStates(n);
            BranchMutations mutations = mutationStore.getBranchMutations(n);
            Arrays.fill(eventStart, -1);
            int k = 0;
            while (k < mutations.getSubstitutionCount()) {
                int end = mutations.getSiteEnd(k);
                eventStart[mutations.getPattern(k)] = k;
                eventEnd[mutations.getPattern(k)] = end;
                k = end;
            }
            for (int p = 0; p < patternCount; p++) {
                if (states[p] < 0) {
                    continue;
                }
                siteLogLikelihoods[p] += eventStart[p] < 0 ?
                        historyModel.getUnchangedSiteLogLikelihood(states[p], length, patternRates[p]) :
                        historyModel.getSiteLogLikelihood(mutations, eventStart[p], eventEnd[p], length, patternRates[p]);
            }
        }
        if (categories != null) {
            categories.completeSiteLogLikelihoods(siteLogLikelihoods, frequencies, rootStates, rootNumber);
        }
        for (int p = 0; p < patternCount; p++) {
            siteLogLikelihoods[p] *= patternWeights[p];
        }
        siteLogLikelihoodsKnown = true;
    }

    double[] getPatternRates() {
        ensureRates();
        return patternRates;
    }

    double getSiteRateScale() {
        return siteRateModel == null ? 1.0 : siteRateModel.getRateForCategory(0);
    }

    /**
     * @return the relative rate of the pattern, i.e. nu_l times the site rate model's rate
     */
    double getSiteRate(int pattern, double scale) {
        return siteRates == null ? scale : scale * siteRates.getParameterValue(pattern);
    }

    boolean hasSiteRates() {
        return siteRates != null;
    }

    /**
     * @return true if patterns may have different rates (per-site rates, or site rate categories)
     */
    boolean hasPatternRates() {
        return siteRates != null || categories != null;
    }

    /**
     * @return the number of site rate categories whose membership is sampled (1 if none)
     */
    public int getCategoryCount() {
        return categories == null ? 1 : categories.getCategoryCount();
    }

    /**
     * @return the rate of each site rate category (after {@link #ensureRates})
     */
    double[] getCategoryRates() {
        return categories.getRates();
    }

    /**
     * @return the rate of the pattern (its site rate, or its category's rate; after {@link #ensureRates})
     */
    double getPatternRate(int pattern) {
        return patternRates[pattern];
    }

    /**
     * Called by a {@link SiteRateCategoryGibbsOperator} that updates the categories.
     */
    void setCategoryOperator() {
        categories.setOperator();
    }

    /**
     * Draws the site rate category of every variable pattern from its full conditional given the history
     * and the tree (a Gibbs update, see {@link SiteRateCategories#resample}).
     */
    void resampleCategories() {
        final long profileStart = MATProfile.start();
        try {
            resampleCategoriesImpl();
        } finally {
            MATProfile.stop(MATProfile.Phase.CATEGORIES_RESAMPLE, profileStart);
        }
    }

    private void resampleCategoriesImpl() {
        if (categories == null) {
            return;
        }
        updateMAT();
        ensureRates();
        categories.resample();
        modelVersion++;
        modelChangedSinceStore = true;
        // index -1: TreeDataLikelihood would otherwise read the index as a node number
        fireModelChanged(mutationStore, -1);
    }

    double[] getPatternWeights() {
        return patternWeights;
    }

    /**
     * @return true if the pattern is constant (every sequence has the same unambiguous state); constant
     * patterns are assumed to have no mutations
     */
    public boolean isConstantPattern(int pattern) {
        return constantPatterns[pattern];
    }

    SubstitutionModel getSubstitutionModel() {
        return substitutionModel;
    }

    FrequencyModel getRootFrequencyModel() {
        return rootFrequencyModel;
    }

    /**
     * @return the length of the branch above the node in expected substitutions (0 for the root)
     */
    double getBranchLength(NodeRef node) {
        NodeRef parent = tree.getParent(node);
        if (parent == null) {
            return 0.0;
        }
        return (tree.getNodeHeight(parent) - tree.getNodeHeight(node)) * branchRateModel.getBranchRate(tree, node);
    }

    @Override
    public double[] getSiteLogLikelihoods() {
        if (!siteLogLikelihoodsKnown) {
            updateMAT();
            computeSiteLogLikelihoods();
        }
        return siteLogLikelihoods.clone();
    }

    @Override
    public void makeDirty() {
        allNodesDirty = true;
        treeChangedSinceCheck = true;
    }

    @Override
    public void storeState() {
        // the heights of the state being stored: those of the nodes touched since the last store changed
        if (allTouched) {
            for (int n = 0; n < knownHeights.length; n++) {
                knownHeights[n] = tree.getNodeHeight(tree.getNode(n));
            }
        } else {
            for (int i = 0; i < touchedNodes.size(); i++) {
                final int n = touchedNodes.get(i);
                knownHeights[n] = tree.getNodeHeight(tree.getNode(n));
            }
        }
        assert heightsKnown() : "the MAT's heights of the accepted state do not match the tree";
        touchedNodes.clear();
        allTouched = false;
        inStep = true;
        movesThisStep = 0;
        pendingLogHastingsRatio = 0.0;
        recomputedNodes.clear();
        mutationStore.commit();
        storedLogLikelihood = logLikelihood;
        modelChangedSinceStore = false;
    }

    private boolean heightsKnown() {
        for (int n = 0; n < knownHeights.length; n++) {
            if (knownHeights[n] != tree.getNodeHeight(tree.getNode(n))) {
                return false;
            }
        }
        return true;
    }

    @Override
    public void restoreState() {
        inStep = false;
        pendingLogHastingsRatio = 0.0;
        // the terms recomputed since the store may have used the rejected branch lengths
        for (int i = 0; i < recomputedNodes.size(); i++) {
            dirtyNodes.add(recomputedNodes.get(i));
        }
        recomputedNodes.clear();
        mutationStore.revert();
        logLikelihood = storedLogLikelihood;
        siteLogLikelihoodsKnown = false;
        if (modelChangedSinceStore) {
            // cached terms may have been computed with the rejected model values
            modelVersion++;
        }
    }

    @Override
    protected void acceptState() {
        // the journal of changes is cleared by the next storeState(), which also takes the accepted heights
        inStep = false;
        pendingLogHastingsRatio = 0.0;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == tree) {
            // The tree likelihood hears of tree changes itself: only note which branch lengths may have
            // changed. Other objects (the height parameters, re-fired by the tree) come with a TreeChangedEvent.
            if (object instanceof TreeChangedEvent) {
                TreeChangedEvent event = (TreeChangedEvent) object;
                if (event.getParameter() == null || !event.isHeightChanged()) {
                    // (a change of a node's height parameter cannot change the topology)
                    treeChangedSinceCheck = true;
                }
                if (event.isNodeChanged() && event.getNode() != null) {
                    NodeRef node = event.getNode();
                    touchedNodes.add(node.getNumber());
                    dirtyNodes.add(node.getNumber());
                    for (int i = 0; i < tree.getChildCount(node); i++) {
                        dirtyNodes.add(tree.getChild(node, i).getNumber());
                    }
                } else {
                    allTouched = true;
                    allNodesDirty = true;
                }
            }
            return;
        }
        if (model == branchRateModel) {
            // branch rates only affect branch lengths
            allNodesDirty = true;
        } else {
            modelVersion++;
            modelChangedSinceStore = true;
        }
        // index -1: TreeDataLikelihood would otherwise read the index as a node number
        fireModelChanged(object, -1);
    }

    @Override
    @SuppressWarnings("rawtypes")
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        modelVersion++;
        modelChangedSinceStore = true;
        fireModelChanged(variable, -1);
    }

    @Override
    public int getTraitCount() {
        return 1;
    }

    @Override
    public int getTraitDim() {
        return patternCount;
    }

    @Override
    public int getPartitionCat() {
        return 1;
    }

    @Override
    public RateRescalingScheme getRateRescalingScheme() {
        return null;
    }

    @Override
    public void setCallback(TreeDataLikelihood treeDataLikelihood) {
    }

    @Override
    public int vectorizeNodeOperations(List<NodeOperation> nodeOperations, int[] operations) {
        throw new UnsupportedOperationException("MAT likelihood does not use node operations");
    }

    @Override
    public void setComputePostOrderStatisticsOnly(boolean computePostOrderStatisticsOnly) {
    }

    @Override
    public boolean providesPostOrderStatisticsOnly() {
        return false;
    }

    @Override
    public PreOrderSettings getPreOrderSettings() {
        return PreOrderSettings.getDefault();
    }

    @Override
    public boolean getPreferGPU() {
        return false;
    }

    @Override
    public boolean getUseAmbiguities() {
        return false;
    }

    @Override
    public PartialsRescalingScheme getRescalingScheme() {
        return PartialsRescalingScheme.NONE;
    }

    @Override
    public boolean getDelayRescalingUntilUnderflow() {
        return false;
    }

    @Override
    public TreeTraversal.TraversalType getOptimalTraversalType() {
        return TreeTraversal.TraversalType.POST_ORDER;
    }

    @Override
    public long getTotalCalculationCount() {
        return calculationCount;
    }

    @Override
    public String getReport() {
        return NAME + ": " + patternCount + " patterns, " + mutationStore.getMutationCount() + " mutations, " +
                mutationStore.getAmbiguityCount() + " ambiguities, rebuilt " + rebuildCount + " times";
    }
}

/*
 * CharacterTreeModel.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
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

package dr.evomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.evolution.tree.FlexibleNode;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evolution.util.Taxon;
import dr.evomodel.tree.DefaultTreeModel;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Attributable;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

/**
 * A character tree whose topology may be multifurcating and may be sampled. It holds a binary skeleton
 * (a DefaultTreeModel with unit branch lengths) and a collapse state: a keep parameter with one entry per internal node
 * of the skeleton (index = node number - S, 0 = collapsed, the root's entry held at 1) when the topology is sampled,
 * or fixed flags when a fixed Newick has a multifurcation. As a Tree it is the contracted tree, the view that is
 * logged: tips keep the skeleton's numbers, internal nodes are numbered S, S + 1, ... in post-order with children
 * ordered by their smallest state, and every branch has length 1. It fires a change only when the contracted tree
 * changes. The skeleton can only be changed by the character-tree operators: an edit that changes the root node or a
 * branch length is refused.
 *
 * @author Xiang Ji
 */
public class CharacterTreeModel extends AbstractModel implements Tree, TreeTraitProvider {

    public static final String CHARACTER_TREE_MODEL = "characterTreeModel";

    private final DefaultTreeModel skeleton;
    private final DataType dataType;
    private final Parameter keep;
    private final CharacterTreeContraction contraction;
    private final boolean variable;
    private final int stateCount;
    private final int rootNumber;
    private final int[] slotOfInternal;
    private final int[] stateOfTip;
    private final int[] minState;

    private FlexibleTree view;
    private int[] skeletonNodeOfView;
    private boolean viewKnown = false;
    private final double[] viewKeep;
    private String signature;
    private HierarchicalSubstitutionModel rateModel;
    private final Attributable.AttributeHelper attributes = new Attributable.AttributeHelper();

    private final TreeTraitProvider.Helper treeTraits = new TreeTraitProvider.Helper(HierarchicalSubstitutionModel.RATE_TRAIT,
            new TreeTrait.D() {
                public String getTraitName() {
                    return HierarchicalSubstitutionModel.RATE_TRAIT;
                }

                public Intent getIntent() {
                    return Intent.NODE;
                }

                public Double getTrait(Tree tree, NodeRef node) {
                    if (tree != CharacterTreeModel.this || tree.isExternal(node) || rateModel == null) {
                        return null;
                    }
                    ensureView();
                    return rateModel.getCladeRate(skeletonNodeOfView[node.getNumber()]);
                }
            });

    /**
     * @param newick   the character tree, possibly multifurcating; tips are the state codes of the data type
     * @param variable whether its topology may change
     * @param keep     the keep parameter of a sampled topology (its dimension is set to S - 1 and its values to 1 at
     *                 the Newick's nodes and 0 at the nodes inserted to resolve a multifurcation), or null
     */
    public CharacterTreeModel(String id, String newick, DataType dataType, boolean variable, Parameter keep) {
        super(CHARACTER_TREE_MODEL);
        setId(id);
        final CharacterTreeContraction.Built built = CharacterTreeContraction.build(id, newick, dataType, !variable, null);
        this.skeleton = built.skeleton;
        AbstractHierarchicalRateProvider.setUnitBranchLengths(skeleton);
        this.dataType = dataType;
        this.variable = variable;
        this.stateCount = skeleton.getExternalNodeCount();
        if (dataType.getStateCount() != stateCount) {
            throw new IllegalArgumentException("character tree " + id + " has " + stateCount + " tips but the data type has "
                    + dataType.getStateCount() + " states");
        }
        this.rootNumber = skeleton.getRoot().getNumber();
        this.stateOfTip = new int[stateCount];
        for (int i = 0; i < stateCount; i++) {
            stateOfTip[i] = AbstractHierarchicalRateProvider.stateOfCode(dataType, skeleton.getTaxonId(i));
        }
        this.minState = new int[skeleton.getNodeCount()];
        this.keep = keep;
        if (keep != null) {
            if (stateCount < 3) {
                throw new IllegalArgumentException("a character tree with fewer than 3 states has no non-root clade to search");
            }
            if (!variable) {
                throw new IllegalArgumentException("keep needs a variable character tree (variable=\"true\")");
            }
            keep.setDimension(stateCount - 1);
            for (int i = 0; i < stateCount - 1; i++) {
                keep.setParameterValueQuietly(i, built.collapsed[stateCount + i] ? 0.0 : 1.0);
            }
            keep.setParameterValueQuietly(rootNumber - stateCount, 1.0);
            keep.addBounds(new Parameter.DefaultBounds(1.0, 0.0, stateCount - 1));
            this.contraction = new CharacterTreeContraction(skeleton, dataType, null, keep);
            this.slotOfInternal = null;
            addVariable(keep);
        } else if (built.isMultifurcating()) {
            if (variable) {
                throw new IllegalArgumentException("the character tree " + id + " is multifurcating and variable, so it needs "
                        + "keep: fixed collapse flags belong to node numbers, which a topology move hands to other clades");
            }
            this.contraction = new CharacterTreeContraction(skeleton, dataType, built.collapsed, null);
            this.slotOfInternal = built.slotOfInternal;
        } else {
            this.contraction = null;
            this.slotOfInternal = null;
        }
        addModel(skeleton);
        skeleton.addModelRestoreListener(this);
        this.viewKeep = keep == null ? null : new double[keep.getDimension()];
        this.signature = computeSignature();
    }

    public DefaultTreeModel getSkeleton() {
        return skeleton;
    }

    public Parameter getKeep() {
        return keep;
    }

    /**
     * The collapse state, or null for a binary tree without keep (every internal node is a clade).
     */
    public CharacterTreeContraction getContraction() {
        return contraction;
    }

    public DataType getDataType() {
        return dataType;
    }

    public boolean isVariable() {
        return variable;
    }

    /**
     * For a fixed multifurcating tree: the post-order index, in the Newick as written, of the clade of each internal
     * node of the skeleton (-1 for an inserted node). Null otherwise.
     */
    public int[] getSlotOfInternal() {
        return slotOfInternal == null ? null : slotOfInternal.clone();
    }

    public boolean isCollapsed(int skeletonNode) {
        return contraction != null && contraction.isCollapsed(skeletonNode);
    }

    /**
     * n_M: the number of internal nodes of the contracted tree, root included.
     */
    public int getContractedInternalNodeCount() {
        return contraction != null ? contraction.getContractedInternalNodeCount() : skeleton.getInternalNodeCount();
    }

    /**
     * log R(M), the log of the number of binary trees that resolve the contracted tree.
     */
    public double getLogResolutionCount() {
        return contraction != null ? contraction.getLogResolutionCount() : 0.0;
    }

    /**
     * The contracted tree as a Newick string, children ordered by smallest state: equal trees give equal strings.
     */
    public String getContractedNewick() {
        ensureView();
        return signature;
    }

    private boolean keepChanged() {
        if (keep == null) {
            return false;
        }
        boolean changed = false;
        for (int i = 0; i < viewKeep.length; i++) {
            final double k = keep.getParameterValue(i);
            if (k != viewKeep[i]) {
                viewKeep[i] = k;
                changed = true;
            }
        }
        return changed;
    }

    public void setRateModel(HierarchicalSubstitutionModel model) {
        if (rateModel != null && rateModel != model) {
            throw new IllegalArgumentException("the character tree " + getId() + " already belongs to a hierarchical "
                    + "substitution model");
        }
        rateModel = model;
    }

    /**
     * Sets every node height of the skeleton to its maximum depth minus its depth, so every branch has length 1,
     * without firing events; for operators that change the topology inside their own tree edit.
     */
    public void resetUnitHeightsQuietly() {
        final int[] depth = new int[skeleton.getNodeCount()];
        final int maxDepth = depths(skeleton.getRoot(), 0, depth);
        for (int i = 0; i < skeleton.getNodeCount(); i++) {
            skeleton.setNodeHeightQuietly(skeleton.getNode(i), maxDepth - depth[i]);
        }
    }

    private int depths(NodeRef node, int d, int[] depth) {
        depth[node.getNumber()] = d;
        int max = d;
        for (int c = 0; c < skeleton.getChildCount(node); c++) {
            max = Math.max(max, depths(skeleton.getChild(node, c), d + 1, depth));
        }
        return max;
    }

    private void checkGuard() {
        if (skeleton.getRoot().getNumber() != rootNumber) {
            throw new IllegalStateException("the character tree " + getId() + " can only be changed by the character-tree "
                    + "operators: an edit moved the root node");
        }
        for (int i = 0; i < skeleton.getNodeCount(); i++) {
            final NodeRef node = skeleton.getNode(i);
            if (!skeleton.isRoot(node) && skeleton.getBranchLength(node) != 1.0) {
                throw new IllegalStateException("the character tree " + getId() + " can only be changed by the character-tree "
                        + "operators: an edit changed a branch length");
            }
        }
    }

    private int fillMinState(NodeRef node) {
        final int n = node.getNumber();
        if (skeleton.isExternal(node)) {
            minState[n] = stateOfTip[n];
        } else {
            int m = Integer.MAX_VALUE;
            for (int c = 0; c < skeleton.getChildCount(node); c++) {
                m = Math.min(m, fillMinState(skeleton.getChild(node, c)));
            }
            minState[n] = m;
        }
        return minState[n];
    }

    private List<Integer> viewChildren(int node) {
        final List<Integer> kids = new ArrayList<Integer>();
        if (contraction != null) {
            kids.addAll(contraction.contractedChildren(node));
        } else {
            final NodeRef ref = skeleton.getNode(node);
            for (int c = 0; c < skeleton.getChildCount(ref); c++) {
                kids.add(skeleton.getChild(ref, c).getNumber());
            }
        }
        Collections.sort(kids, new Comparator<Integer>() {
            public int compare(Integer a, Integer b) {
                return Integer.compare(minState[a], minState[b]);
            }
        });
        return kids;
    }

    private String computeSignature() {
        fillMinState(skeleton.getRoot());
        final StringBuilder sb = new StringBuilder();
        writeSignature(rootNumber, sb);
        return sb.append(';').toString();
    }

    private void writeSignature(int node, StringBuilder sb) {
        if (node < stateCount) {
            sb.append(skeleton.getTaxonId(node));
            return;
        }
        final List<Integer> kids = viewChildren(node);
        sb.append('(');
        for (int k = 0; k < kids.size(); k++) {
            if (k > 0) {
                sb.append(',');
            }
            writeSignature(kids.get(k), sb);
        }
        sb.append(')');
    }

    private void ensureView() {
        if (keepChanged()) {
            viewKnown = false;
        }
        if (viewKnown) {
            return;
        }
        signature = computeSignature();
        final List<Integer> skeletonNodes = new ArrayList<Integer>();
        final FlexibleNode root = buildView(rootNumber, skeletonNodes);
        final int maxDepth = setDepths(root, 0);
        setHeights(root, 0, maxDepth);
        final Map<Taxon, Integer> taxonNumbers = new HashMap<Taxon, Integer>();
        for (int i = 0; i < stateCount; i++) {
            taxonNumbers.put(skeleton.getNodeTaxon(skeleton.getExternalNode(i)), skeleton.getExternalNode(i).getNumber());
        }
        view = new FlexibleTree(root, taxonNumbers);
        view.setId(getId());
        skeletonNodeOfView = new int[view.getNodeCount()];
        for (int i = 0; i < view.getNodeCount(); i++) {
            skeletonNodeOfView[i] = (Integer) ((FlexibleNode) view.getNode(i)).getAttribute("skeletonNode");
        }
        viewKnown = true;
    }

    private FlexibleNode buildView(int node, List<Integer> skeletonNodes) {
        final FlexibleNode out = new FlexibleNode();
        out.setAttribute("skeletonNode", node);
        if (node < stateCount) {
            out.setTaxon(skeleton.getNodeTaxon(skeleton.getNode(node)));
            return out;
        }
        for (int kid : viewChildren(node)) {
            out.addChild(buildView(kid, skeletonNodes));
        }
        return out;
    }

    private int setDepths(FlexibleNode node, int d) {
        int max = d;
        for (int c = 0; c < node.getChildCount(); c++) {
            max = Math.max(max, setDepths(node.getChild(c), d + 1));
        }
        return max;
    }

    private void setHeights(FlexibleNode node, int d, int maxDepth) {
        node.setHeight(maxDepth - d);
        node.setLength(d == 0 ? 0.0 : 1.0);
        for (int c = 0; c < node.getChildCount(); c++) {
            setHeights(node.getChild(c), d + 1, maxDepth);
        }
    }

    private Tree view() {
        ensureView();
        return view;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == skeleton) {
            checkGuard();
            if (HierarchicalSubstitutionModel.changesTopology(object)) {
                contractedTreeMayHaveChanged();
            }
        }
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == keep) {
            if (keep.getParameterValue(rootNumber - stateCount) != 1.0) {
                throw new IllegalStateException("the root of the character tree " + getId() + " cannot be collapsed");
            }
            contractedTreeMayHaveChanged();
        }
    }

    private void contractedTreeMayHaveChanged() {
        viewKnown = false;
        final String now = computeSignature();
        if (!now.equals(signature)) {
            signature = now;
            fireModelChanged();
        }
    }

    @Override
    public void modelRestored(Model model) {
        if (model == skeleton) {
            viewKnown = false;
            signature = computeSignature();
        }
    }

    @Override
    protected void storeState() {
    }

    @Override
    protected void restoreState() {
        viewKnown = false;
        signature = computeSignature();
    }

    @Override
    protected void acceptState() {
    }

    @Override
    public TreeTrait[] getTreeTraits() {
        return treeTraits.getTreeTraits();
    }

    @Override
    public TreeTrait getTreeTrait(String key) {
        return treeTraits.getTreeTrait(key);
    }

    public NodeRef getRoot() {
        return view().getRoot();
    }

    public int getNodeCount() {
        return view().getNodeCount();
    }

    public NodeRef getNode(int i) {
        return view().getNode(i);
    }

    public NodeRef getInternalNode(int i) {
        return view().getInternalNode(i);
    }

    public NodeRef getExternalNode(int i) {
        return view().getExternalNode(i);
    }

    public int getExternalNodeCount() {
        return stateCount;
    }

    public int getInternalNodeCount() {
        return view().getInternalNodeCount();
    }

    public Taxon getNodeTaxon(NodeRef node) {
        return view().getNodeTaxon(node);
    }

    public boolean hasNodeHeights() {
        return true;
    }

    public double getNodeHeight(NodeRef node) {
        return view().getNodeHeight(node);
    }

    public boolean hasBranchLengths() {
        return true;
    }

    public double getBranchLength(NodeRef node) {
        return view().getBranchLength(node);
    }

    public double getNodeRate(NodeRef node) {
        return 1.0;
    }

    public Object getNodeAttribute(NodeRef node, String name) {
        return null;
    }

    public Iterator getNodeAttributeNames(NodeRef node) {
        return Collections.emptyList().iterator();
    }

    public boolean isExternal(NodeRef node) {
        return view().isExternal(node);
    }

    public boolean isRoot(NodeRef node) {
        return view().isRoot(node);
    }

    public int getChildCount(NodeRef node) {
        return view().getChildCount(node);
    }

    public NodeRef getChild(NodeRef node, int j) {
        return view().getChild(node, j);
    }

    public NodeRef getParent(NodeRef node) {
        return view().getParent(node);
    }

    public Tree getCopy() {
        return view().getCopy();
    }

    public int getTaxonCount() {
        return skeleton.getTaxonCount();
    }

    public Taxon getTaxon(int taxonIndex) {
        return skeleton.getTaxon(taxonIndex);
    }

    public String getTaxonId(int taxonIndex) {
        return skeleton.getTaxonId(taxonIndex);
    }

    public int getTaxonIndex(String id) {
        return skeleton.getTaxonIndex(id);
    }

    public int getTaxonIndex(Taxon taxon) {
        return skeleton.getTaxonIndex(taxon);
    }

    public List<Taxon> asList() {
        return skeleton.asList();
    }

    public Object getTaxonAttribute(int taxonIndex, String name) {
        return skeleton.getTaxonAttribute(taxonIndex, name);
    }

    public Iterator<Taxon> iterator() {
        return skeleton.iterator();
    }

    public Type getUnits() {
        return skeleton.getUnits();
    }

    public void setUnits(Type units) {
        skeleton.setUnits(units);
    }

    public void setAttribute(String name, Object value) {
        attributes.setAttribute(name, value);
    }

    public Object getAttribute(String name) {
        return attributes.getAttribute(name);
    }

    public Iterator<String> getAttributeNames() {
        return attributes.getAttributeNames();
    }

    @Override
    public String toString() {
        return getContractedNewick();
    }
}

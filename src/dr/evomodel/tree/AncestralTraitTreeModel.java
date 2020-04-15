/*
 * TransformedTreeModel.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.tree;

import dr.evolution.tree.*;
import dr.evolution.util.MutableTaxonListListener;
import dr.evolution.util.Taxon;
import dr.evomodel.continuous.AncestralTaxonInTree;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;

import java.util.*;

/**
 * A tree model with additional ancestral taxon nodes for peeling
 *
 * @author Marc A. Suchard
 */
public class AncestralTraitTreeModel extends AbstractModel implements MutableTreeModel, TransformableTree, Citable {

    private static final boolean DEBUG = false;

    private final int ancestorCount;

    private final int treeExternalCount;
    private final int treeInternalCount;

    private final int externalCount;
    private final int internalCount;

    private int extraInternal;
    
    private ShadowNode[] nodes;
    private ShadowNode[] storedNodes;

    private ShadowNode root;
    private int storedRootNumber;

    @Override
    public NodeRef getOriginalNode(NodeRef transformedNode) {

        assert (transformedNode != null);

        checkShadowTree();
        ShadowNode node = (ShadowNode) transformedNode;
        return node.getOriginalNode();
    }

    @Override
    public NodeRef getTransformedNode(NodeRef originalNode) {
        return new BasicNode(mapOriginalToShadowNumber(originalNode.getNumber()));
    }

    public class ShadowNode implements NodeRef {

        private int number = -1;
        private int originalNumber = -1;

        private AncestralTaxonInTree ancestor = null;

        private ShadowNode child0 = null;
        private ShadowNode child1 = null;
        private ShadowNode parent = null;

        private boolean used = false;

        private ShadowNode() { /* Do nothing */ }

        private ShadowNode(int number, NodeRef originalNode, AncestralTaxonInTree ancestor) {
            this.number = number;
            this.originalNumber = originalNode != null ?
                    originalNode.getNumber() :
                    -1;
            this.ancestor = ancestor;
            this.used = true;
        }

        private void adoptValues(ShadowNode donor, ShadowNode[] copy) {
            this.number = donor.number;
            this.originalNumber = donor.originalNumber;
            this.ancestor = donor.ancestor;
            this.used = donor.used;

            if (donor.child0 != null) {
                this.child0 = copy[donor.child0.getNumber()];
            } else {
                this.child0 = null;
            }

            if (donor.child1 != null) {
                this.child1 = copy[donor.child1.getNumber()];
            } else {
                this.child1 = null;
            }

            if (donor.parent != null) {
                this.parent = copy[donor.parent.getNumber()];
            } else {
                this.parent = null;
            }
        }

        @Override
        public int getNumber() { return number; }

        @Override
        public void setNumber(int n) {
            throw new RuntimeException("Node number is not modifiable");
        }

        private int getOriginalNumber() { return originalNumber; }

        private NodeRef getOriginalNode() {
            return originalNumber >= 0 ?
                    treeModel.getNode(originalNumber) :
                    null;
        }

        private NodeRef getChild(int i) {
            if (i == 0) {
                return child0;
            } else if (i == 1) {
                return child1;
            } else {
                throw new IllegalArgumentException("Binary trees only!");
            }
        }

        private boolean isExternal() { return child0 == null && child1 == null; }

        public String toString() {
            int pa = parent != null ? parent.getNumber() : -1;
            int c0 = child0 != null ? child0.getNumber() : -1;
            int c1 = child1 != null ? child1.getNumber() : -1;
            String anc = ancestor != null ?  ancestor.getTaxon().getId() : "-1";
            String u = used ? "true" : "false";
            double height = getNodeHeight(this);
            boolean ex = isExternal();
            int cnt = getChildCount(this);

            return "node " + number + " " + pa + " " + c0 + " " + c1 + " : " + originalNumber + " " + anc + " " + u +
                    " " + height + " " + ex + " " + cnt;
        }

        private boolean isUsed() {
            return used;
        }

        private void setUnused() {
            this.used = false;
        }
    }

    public static String toString(ShadowNode[] nodes, int root) {
        StringBuilder sb = new StringBuilder();
        for (ShadowNode node : nodes) {
            if (node != null) {
                sb.append(node.toString()).append("\n");
            } else {
                sb.append("null\n");
            }
        }
        sb.append("root = ").append(root);
        return sb.toString();
    }

    public AncestralTraitTreeModel(String id, MutableTreeModel tree, List<AncestralTaxonInTree> ancestors) {
        super(id);
        this.treeModel = tree;
        this.ancestors = ancestors;

        ancestorCount = ancestors.size();

        treeExternalCount = treeModel.getExternalNodeCount();
        treeInternalCount = treeModel.getInternalNodeCount();

        externalCount = treeExternalCount + ancestorCount;
        internalCount = treeInternalCount + ancestorCount;

        addModel(tree);

        int index = 0;
        for (AncestralTaxonInTree ancestor : ancestors) {
            addRestrictedPartials(ancestor, index);
            ++index;
        }

        nodes = new ShadowNode[externalCount + internalCount];
        for (int i = 0; i < nodes.length; ++i) {
            nodes[i] = new ShadowNode();
        }
    }

    public Tree getOriginalTree() {
        return treeModel;
    }

    private void checkShadowTree() {
        if (!validShadowTree) {
            buildShadowTree();
            validShadowTree = true;
        }
    }

    private void buildShadowTree() {

        if (DEBUG) {
            System.out.println("ATTM.bST");
        }

        setupClamps();

        for (ShadowNode node : nodes) { // TODO Only need to set extra nodes
            node.setUnused();
        }

        extraInternal = 0;
        root = buildRecursivelyShadowTree(treeModel.getRoot(), null);
    }

    private void storeNode(ShadowNode node) {
        nodes[node.getNumber()] = node;
    }

    private static void sortByTime(List<AncestralTaxonInTree> ancestors) {
        ancestors.sort((lhs, rhs) -> -Double.compare(lhs.getHeight(), rhs.getHeight()));
    }

    private ShadowNode buildRecursivelyShadowTree(NodeRef originalNode,
                                                  ShadowNode parentNode) {
        final int originalNumber = originalNode.getNumber();
        final int newNumber = mapOriginalToShadowNumber(originalNumber);

        ShadowNode newNode = new ShadowNode(newNumber, originalNode, null);
        newNode.parent = parentNode;
        storeNode(newNode);

        NodeRef originalChild0 = treeModel.getChild(originalNode, 0);
        NodeRef originalChild1 = treeModel.getChild(originalNode, 1);

        ShadowNode recurse0 = newNode;
        ShadowNode recurse1 = newNode;

        if (nodeToClampMap.containsKey(originalNode.getNumber())) {

            List<AncestralTaxonInTree> ancestors = nodeToClampMap.get(originalNode.getNumber());

            if (ancestors.size() > 1) {
                sortByTime(ancestors);
            }

            // Add tips
            for (AncestralTaxonInTree ancestor : ancestors) {
                final int newTipNumber = treeExternalCount + ancestor.getIndex();

                ShadowNode newTipNode = new ShadowNode(newTipNumber, null, ancestor);

                ShadowNode newInternalNode = new ShadowNode(externalCount + treeInternalCount + extraInternal,
                        null, ancestor);

                if (ancestor.getPathChildNumber() == 0) {

                    recurse0.child0 = newInternalNode;
                    newInternalNode.parent = recurse0;

                    newInternalNode.child1 = newTipNode;
                    newTipNode.parent = newInternalNode;

                    recurse0 = newInternalNode; // recurse0.child0 is free

                } else {

                    recurse1.child1 = newInternalNode;
                    newInternalNode.parent = recurse1;

                    newInternalNode.child0 = newTipNode;
                    newTipNode.parent = newInternalNode;

                    recurse1 = newInternalNode; // recurse1.child1 is free
                }

                storeNode(newTipNode);
                storeNode(newInternalNode);

                ++extraInternal;
            }
        }

        if (!treeModel.isExternal(originalNode)) {
            recurse0.child0 = buildRecursivelyShadowTree(originalChild0, recurse0);
            recurse1.child1 = buildRecursivelyShadowTree(originalChild1, recurse1);
        }

        return newNode;
    }

    private int mapOriginalToShadowNumber(int originalNumber) {

        assert (originalNumber >= 0);
        assert (originalNumber < treeExternalCount + treeInternalCount);

        int newNumber = originalNumber;
        if (originalNumber >= treeExternalCount) {
            newNumber += ancestorCount;
        }
        return newNumber;
    }

    public String toString() {
        return TreeUtils.newick(this);
    }

    public boolean isVariable() {
        return treeModel instanceof AbstractModel && ((AbstractModel) treeModel).isVariable();
    }

    public double getNodeHeight(NodeRef iNode) {
        assert (iNode != null);

        checkShadowTree();

        double height;

        ShadowNode node = (ShadowNode) iNode;
        int originalNumber = node.getOriginalNumber();
        if (originalNumber >= 0) {
            height = treeModel.getNodeHeight(node.getOriginalNode());
        } else {

            final AncestralTaxonInTree ancestor = node.ancestor;

            if (ancestor.isOnAncestralPath()) { // TODO Refactor into subclasses
                double ancestorHeight = ancestor.getHeight();
                if (node.isExternal()) {
                    ancestorHeight -= ancestor.getPseudoBranchLength();
                }
                double rootHeight = treeModel.getNodeHeight(treeModel.getRoot());
                ancestorHeight = Math.min(rootHeight, ancestorHeight);
                height = ancestorHeight;
            } else { // Below is the original
                if (node.isExternal()) {
                    height = treeModel.getNodeHeight(node.parent.parent.getOriginalNode())
                            - ancestor.getPseudoBranchLength();
                } else {
                    height = treeModel.getNodeHeight(node.parent.getOriginalNode());
                }
            }
        }

        return height;
    }

    private static final boolean FIX_BRANCH_LENGTH = true;

    public double getBranchLength(NodeRef iNode) {

        if (FIX_BRANCH_LENGTH) {
            throw new RuntimeException("Not yet implemented");
        }

        System.err.println("DIE");

        assert (iNode != null);

        checkShadowTree();

        ShadowNode node = (ShadowNode) iNode;
        if (!node.isUsed()) {
            return 0.0;
        }

        int originalNumber = node.getOriginalNumber();
        if (originalNumber >= 0) {
            return treeModel.getBranchLength(node.getOriginalNode());
        } else {
            if (node.isExternal()) {
                return node.ancestor.getPseudoBranchLength();
            } else {
                return 0;
            }
        }
    }

    private void storeNodeStructure() {

        final int length = nodes.length;

        // Initialize once
        if (storedNodes == null) {
            storedNodes = new ShadowNode[length];
            for (int i = 0; i < length; ++i) {
                storedNodes[i] = new ShadowNode();
            }
        }

        // Copy
        for (int i = 0; i < length; ++i) {
            storedNodes[i].adoptValues(nodes[i], storedNodes);
        }
    }
    
    protected void storeState() {

        assert (nodes != null);

        savedValidShadowTree = validShadowTree;

        if (validShadowTree) {
            storeNodeStructure();
            storedRootNumber = root.getNumber();
        }
    }

    /**
     * Restore the stored state
     */
    protected void restoreState() {

        assert (storedNodes != null);

        validShadowTree = savedValidShadowTree;

        if (validShadowTree) {
            ShadowNode[] tmp = nodes;
            nodes = storedNodes;
            storedNodes = tmp;

            root = nodes[storedRootNumber];
        }
    }

    private class RemappedTreeChangeEvent implements TreeChangedEvent {

        final private TreeChangedEvent event;
        final private NodeRef node;

        private RemappedTreeChangeEvent(TreeChangedEvent event, NodeRef node) {
            this.event = event;
            this.node = node;
        }

        @Override public int getIndex() { return event.getIndex(); }

        @Override public NodeRef getNode() { return node; }

        @Override public Parameter getParameter() { return event.getParameter(); }

        @Override public boolean isNodeChanged() { return event.isNodeChanged(); }

        @Override public boolean isTreeChanged() { return event.isTreeChanged(); }

        @Override public boolean isNodeParameterChanged() { return event.isNodeParameterChanged(); }

        @Override public boolean isHeightChanged() { return event.isHeightChanged(); }
    }

    protected void handleModelChangedEvent(Model model, Object object, int index) {

        if (DEBUG) {
            System.out.println("ATTM.hMCE hit with " + object.getClass().getCanonicalName());
        }

        if (model == treeModel) {

            if (object instanceof TreeChangedEvent) {
                final TreeChangedEvent treeChangedEvent = (TreeChangedEvent) object;

                if (treeChangedEvent.isTreeChanged()) {

                    if (DEBUG) {
                        System.out.println("\tATTM.hMCE invalidate");
                    }

                    validShadowTree = false;
                    fireModelChanged(new TreeChangedEvent.WholeTree());

                } else if (treeChangedEvent.isNodeChanged()) {

                    if (DEBUG) {
                        System.out.println("\tATTM.hMCE nodeChange");
                    }

                    final NodeRef originalNode = treeChangedEvent.getNode();
                    if (originalNode != null) { // Remap
                        ShadowNode shadow = nodes[mapOriginalToShadowNumber(originalNode.getNumber())];

                        if (DEBUG) {
                            System.out.println("\tATTM.hMCE isRoot? " + isRoot(shadow) + " " + treeModel.isRoot(originalNode));
                        }

                        if (!shadow.isExternal()) {
                            if (shadow.child0.ancestor != null) {
                                /* If there is an ancestor (and zero-branch-length internal node),
                                 * we apparently need to hit the both the parent of the ancestor and the
                                 * zero-branch-length internal node.  MAS is unsure why we must hit the
                                 * ancestor parent.
                                 */

                                if (shadow.parent != null) {
                                    fireModelChanged(new RemappedTreeChangeEvent(treeChangedEvent, shadow.parent), index);
                                }

                                fireModelChanged(new RemappedTreeChangeEvent(treeChangedEvent, shadow.child0), index);
                                shadow = shadow.child1;
                            }
                        }

                        object = new RemappedTreeChangeEvent(treeChangedEvent, shadow);

                    }
                    fireModelChanged(object, index);

                } else {
                    throw new IllegalArgumentException("TreeModel should not generate other events");
                }

            } else if (object instanceof Parameter) {

                if (hasAncestralPathTaxa) {
                    validShadowTree = false;
                    fireModelChanged(new TreeChangedEvent.WholeTree());
                }
            } else {
                throw new IllegalArgumentException("TreeModel should not generate other objects");
            }

        } else if (model instanceof AncestralTaxonInTree && ancestors.contains(model)) {

            if (hasAncestralPathTaxa) {
                validShadowTree = false;
            }
            fireModelChanged(new TreeChangedEvent.WholeTree());

        } else {
            throw new IllegalArgumentException("Illegal model");
        }

        if (DEBUG) {
            System.out.println("ATTM.hMCE end");
        }
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        // Do nothing; no variables
    }

    // Delegate the rest to treeModel

    public NodeRef getRoot() {
        checkShadowTree();
        return root;
    }

    public int getNodeCount() {
        return externalCount + internalCount;
    }

    public NodeRef getNode(int i) {
        assert (i >= 0 && i < externalCount + internalCount);
        
        checkShadowTree();
        return nodes[i];
    }

    public NodeRef getInternalNode(int i) {
        assert (i >= 0 && i < internalCount);

        checkShadowTree();
        return nodes[i + externalCount];
    }

    public NodeRef getExternalNode(int i) {
        assert (i >= 0 && i < externalCount);

        checkShadowTree();
        return nodes[i];
    }

    public int getExternalNodeCount() {
        return externalCount;
    }

    public int getInternalNodeCount() {
        return internalCount;
    }

    public Taxon getNodeTaxon(NodeRef node) {
        assert (node != null);

        checkShadowTree();

        int originalNumber = ((ShadowNode) node).getOriginalNumber();
        if (originalNumber >= 0) {
            return treeModel.getNodeTaxon(treeModel.getNode(originalNumber));
        } else {
            return getTaxonByTreeIndex(node.getNumber());
        }
    }

    public boolean hasNodeHeights() {
        return treeModel.hasNodeHeights();
    }

    public boolean hasBranchLengths() {
        return treeModel.hasBranchLengths();
    }

    public double getNodeRate(NodeRef node) {
        throw new RuntimeException("Not yet implemented");
    }

    public Object getNodeAttribute(NodeRef node, String name) {
        throw new RuntimeException("Not yet implemented");
    }

    public Iterator getNodeAttributeNames(NodeRef node) {
        throw new RuntimeException("Not yet implemented");
    }

    public boolean isExternal(NodeRef node) {
        assert (node != null);

        checkShadowTree();
        return ((ShadowNode) node).isExternal();
    }

    public boolean isRoot(NodeRef node) {
        assert (node != null);

        checkShadowTree();
        return node == root;
    }

    public int getChildCount(NodeRef node) {
        assert (node != null);

        checkShadowTree();
        if (((ShadowNode) node).isExternal()) {
            return 0;
        } else {
            return 2;
        }
    }

    public NodeRef getChild(NodeRef node, int j) {
        assert (node != null);
        
        checkShadowTree();
        return ((ShadowNode) node).getChild(j);
    }

    public NodeRef getParent(NodeRef node) {
        assert (node != null);
        
        checkShadowTree();
        return ((ShadowNode) node).parent;
    }

    public Tree getCopy() {
        throw new RuntimeException("Not yet implemented");
    }

    protected void acceptState() {
        // Do nothing
    }

    private final MutableTreeModel treeModel;
    private final List<AncestralTaxonInTree> ancestors;

    public double[] getMultivariateNodeTrait(NodeRef node, String name) {
        throw new RuntimeException("Not yet implemented");
    }

    public void setMultivariateTrait(NodeRef n, String name, double[] value) {
        throw new RuntimeException("Not yet implemented");
    }

    public boolean beginTreeEdit() {
        throw new RuntimeException("Not yet implemented");
    }

    public void endTreeEdit() {
        throw new RuntimeException("Not yet implemented");
    }

    public void addChild(NodeRef parent, NodeRef child) {
        throw new RuntimeException("Not yet implemented");
    }

    public void removeChild(NodeRef parent, NodeRef child) {
        throw new RuntimeException("Not yet implemented");
    }

    public void replaceChild(NodeRef node, NodeRef child, NodeRef newChild) {
        throw new RuntimeException("Not yet implemented");
    }

    public void setRoot(NodeRef root) {
        throw new RuntimeException("Not yet implemented");
    }

    public void setNodeHeight(NodeRef node, double height) {
        throw new RuntimeException("Not yet implemented");
    }

    public void setNodeRate(NodeRef node, double height) {
        throw new RuntimeException("Not yet implemented");
    }

    public void setBranchLength(NodeRef node, double length) {
        throw new RuntimeException("Not yet implemented");
    }

    public void setNodeAttribute(NodeRef node, String name, Object value) {
        throw new RuntimeException("Not yet implemented");
    }

    public void addMutableTreeListener(MutableTreeListener listener) {
        throw new RuntimeException("Not yet implemented");
    }

    public void setAttribute(String name, Object value) {
        throw new RuntimeException("Not yet implemented");
    }

    public Object getAttribute(String name) {
        throw new RuntimeException("Not yet implemented");
    }

    public Iterator<String> getAttributeNames() {
        throw new RuntimeException("Not yet implemented");
    }

    public int addTaxon(Taxon taxon) {
        throw new RuntimeException("Not yet implemented");
    }

    public boolean removeTaxon(Taxon taxon) {
        throw new RuntimeException("Not yet implemented");
    }

    public void setTaxonId(int taxonIndex, String id) {
        treeModel.setTaxonId(taxonIndex, id);
        throw new RuntimeException("Not yet implemented");
    }

    public void setTaxonAttribute(int taxonIndex, String name, Object value) {
        throw new RuntimeException("Not yet implemented");
    }

    public void addMutableTaxonListListener(MutableTaxonListListener listener) {
        throw new RuntimeException("Not yet implemented");
    }

    public int getTaxonCount() {
        return treeModel.getTaxonCount() + ancestorCount;
    }

    public Taxon getTaxon(int taxonIndex) {
        Taxon taxon;
        if (taxonIndex < treeExternalCount) {
            taxon = treeModel.getTaxon(taxonIndex);
        } else {
            taxon = getTaxonByTreeIndex(taxonIndex);
        }
        return taxon;
    }

    public String getTaxonId(int taxonIndex) {
        String result;
        if (taxonIndex < treeExternalCount) {
            result = treeModel.getTaxonId(taxonIndex);
        } else {
            Taxon taxon = getTaxonByTreeIndex(taxonIndex);
            result = taxon.getId();
        }
        return result;
    }

    public int getTaxonIndex(String id) {
        for (int i = 0, n = getTaxonCount(); i < n; i++) {
            if (getTaxonId(i).equals(id)) return i;
        }
        return -1;
    }

    public int getTaxonIndex(Taxon taxon) {

        int index = treeModel.getTaxonIndex(taxon);
        if (index != -1) {
            return index;
        } else {
            throw new RuntimeException("Not yet implemented");
        }
    }

    public List<Taxon> asList() {
        List<Taxon> taxa = new ArrayList<>();
        for (int i = 0, n = getTaxonCount(); i < n; i++) {
            taxa.add(getTaxon(i));
        }
        return taxa;
    }

    public Object getTaxonAttribute(int taxonIndex, String name) {
        if (taxonIndex < treeModel.getExternalNodeCount()) {
            return treeModel.getTaxonAttribute(taxonIndex, name);
        } else {
            Taxon taxon = getTaxonByTreeIndex(taxonIndex);
            if (taxon != null) {
                return taxon.getAttribute(name);
            }
            return null;
        }
    }

    public Iterator<Taxon> iterator() {
        throw new RuntimeException("Not yet implemented");
    }

    public Type getUnits() {
        return treeModel.getUnits();
    }

    public void setUnits(Type units) {
        treeModel.setUnits(units);
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TRAIT_MODELS;
    }

    @Override
    public String getDescription() {
        return "Bayesian estimation of Pagel's lambda";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CommonCitations.VRANCKEN_2015_SIMULTANEOUSLY);
    }

    private Taxon getTaxonByTreeIndex(int index) {
        return ancestors.get(index - treeExternalCount).getTaxon();
    }

    private void setupClamps() {
        nodeToClampMap.clear();
        recursiveSetupMrcaClamps(treeModel, treeModel.getRoot(), new BitSet(), clampList, nodeToClampMap);
        setupAncestralPathClamps(treeModel, clampList, nodeToClampMap);
    }

    private static void addAncestralTaxonToMap(Map<Integer, List<AncestralTaxonInTree>> map,
                                               int index,
                                               AncestralTaxonInTree ancestor) {

        final List<AncestralTaxonInTree> list = map.containsKey(index) ? map.get(index) : new ArrayList<>();
        list.add(ancestor);
        map.put(index, list);
    }

    private void setupAncestralPathClamps(Tree tree,
                                          Map<BitSet, AncestralTaxonInTree> clampList,
                                          Map<Integer, List<AncestralTaxonInTree>> nodeToClampMap) {

        hasAncestralPathTaxa = false;

        for (int i = 0; i < tree.getExternalNodeCount(); ++i) {
            NodeRef node = tree.getExternalNode(i);
            BitSet tip = new BitSet();
            tip.set(node.getNumber());

            if (clampList.containsKey(tip)) {
                AncestralTaxonInTree partials = clampList.get(tip);

                double pathHeight = partials.getHeight();
                assert (pathHeight > 0.0);

                NodeRef parent = tree.getParent(node);
                double parentHeight = tree.getNodeHeight(parent);
                boolean isChild0 = tree.getChild(parent, 0) == node;

                while (parentHeight < pathHeight && parent != tree.getRoot()) {
                    node = parent;
                    parent = tree.getParent(parent);
                    parentHeight = tree.getNodeHeight(parent);
                    isChild0 = tree.getChild(parent, 0) == node;
                }

                partials.setNode(parent, isChild0 ? 0 : 1);
                addAncestralTaxonToMap(nodeToClampMap, parent.getNumber(), partials);
                hasAncestralPathTaxa = true;
            }
        }
    }

    private static void recursiveSetupMrcaClamps(Tree tree, NodeRef node,
                                                 BitSet tips,
                                                 Map<BitSet, AncestralTaxonInTree> clampList,
                                                 Map<Integer, List<AncestralTaxonInTree>> nodeToClampMap) {

        if (tree.isExternal(node)) {
            tips.set(node.getNumber());
        } else {
            for (int i = 0; i < tree.getChildCount(node); i++) {
                NodeRef child = tree.getChild(node, i);

                BitSet childTips = new BitSet();
                recursiveSetupMrcaClamps(tree, child, childTips, clampList, nodeToClampMap);
                tips.or(childTips);
            }

            if (clampList.containsKey(tips)) {
                AncestralTaxonInTree partials = clampList.get(tips);
                partials.setNode(node);
                addAncestralTaxonToMap(nodeToClampMap, node.getNumber(), partials);
            }
        }
    }

    private void addRestrictedPartials(AncestralTaxonInTree nodeClamp, int index) {

        clampList.put(nodeClamp.getTipBitSet(), nodeClamp);
        addModel(nodeClamp);
        nodeClamp.setIndex(index);
    }

    final private Map<BitSet, AncestralTaxonInTree> clampList = new HashMap<>();
    final private Map<Integer, List<AncestralTaxonInTree>> nodeToClampMap = new HashMap<>();

    private boolean hasAncestralPathTaxa = false;

    private boolean validShadowTree = false;
    private boolean savedValidShadowTree;
}

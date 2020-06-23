/*
 * CompatibilityStatistic.java
 *
 * Copyright (c) 2002-2020 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.Taxon;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.*;

/**
 * Tests whether a tree is compatible with a less-resolved constraints tree. The constraints tree does not need to
 * include every taxa in the test tree, but every taxa in the constraints tree must be in the test tree. Taxa not in the
 * constraints tree are ignored in the comparison.
 * *
 *
 * @author Andrew Rambaut
 */
public class ConstraintsTreeLikelihood extends AbstractModelLikelihood {

    public ConstraintsTreeLikelihood(String name, Tree targetTree, Tree constraintsTree) throws TreeUtils.MissingTaxonException {

        super(name);

        for (int i = 0; i < constraintsTree.getTaxonCount(); i++) {
            String id = constraintsTree.getTaxonId(i);
            if (targetTree.getTaxonIndex(id) == -1) {
                throw new TreeUtils.MissingTaxonException(constraintsTree.getTaxon(i));
            }
            Taxon taxon = targetTree.getTaxon(targetTree.getTaxonIndex(id));
            for (int j = 0; j <targetTree.getExternalNodeCount() ; j++) {
                NodeRef tip = targetTree.getExternalNode(j);
                if(targetTree.getNodeTaxon(tip).equals(taxon)){
                    constrainedTips.add(tip.getNumber());
                    break;
                }
            }
        }

        setupClades(constraintsTree, constraintsTree.getRoot(), targetTree);

        updateNode = new boolean[targetTree.getNodeCount()];
        targetTreeNodeCladeMap = new BitSet[targetTree.getNodeCount()];
        lastUpdated = new ArrayList<>();

        for (int i = 0; i < updateNode.length; i++) {
            updateNode[i] = true;
        }

        if (targetTree instanceof TreeModel) {
            addModel((TreeModel)targetTree);
        }
        this.targetTree = targetTree;
    }

    private void updateAllNodes(){
        for (int i = 0; i < updateNode.length; i++) {
            updateNode[i] = true;
            targetClades.clear();
        }
        likelihoodKnown = false;

    }

    private void updateNodeAndAncestors(NodeRef node){
        while(node!=null){
            int nodeIndex = node.getNumber();
            updateNode[nodeIndex] = true;
            //Sometimes this has already been removed - why?
            targetClades.remove(targetTreeNodeCladeMap[nodeIndex]);
            node = targetTree.getParent(node);

            lastUpdated.add(nodeIndex);
        }
        likelihoodKnown = false;
    }


    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
//        makeDirty();
            if (object instanceof TreeChangedEvent) {

                if (((TreeChangedEvent) object).isNodeChanged()) {
                    // If a node event occurs the node and its two child nodes
                    // are flagged for updating (this will result in everything
                    // above being updated as well. Node events occur when a node
                    // is added to a branch, removed from a branch or its height or
                    // rate changes.
                    updateNodeAndAncestors(((TreeChangedEvent) object).getNode());

                } else if (((TreeChangedEvent) object).isTreeChanged()) {
                    // Full tree events result in a complete updating of the tree likelihood
                    // This event type is now used for EmpiricalTreeDistributions.
//                    System.err.println("Full tree update event - these events currently aren't used\n" +
//                            "so either this is in error or a new feature is using them so remove this message.");
                    updateAllNodes();
                } else {
                    // Other event types are ignored (probably trait changes).
                    //System.err.println("Another tree event has occured (possibly a trait change).");
                }
            }
    }

    @Override
    protected void storeState() {
        storedLikelihoodKnown = likelihoodKnown;
        storedLogLikelihood = logLikelihood;
        lastUpdated = new ArrayList<>();
    }

    @Override
    protected void restoreState() {
        likelihoodKnown = storedLikelihoodKnown;
        logLikelihood = storedLogLikelihood;
        for (int i :lastUpdated){
            updateNode[i] = true;
            targetClades.remove(targetTreeNodeCladeMap[i]);
        }
    }

    @Override
    protected void acceptState() {
        // do nothing

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        // do nothing
    }

    @Override
    public Model getModel() {
        return this;
    }

    @Override
    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = isCompatible() ? 0.0 : Double.NEGATIVE_INFINITY;
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    @Override
    public void makeDirty() {
        likelihoodKnown = false;
        updateAllNodes();
    }

    /**
     * Returns true if all the clades in the constraints tree are present in the target tree
     * @return
     */
    private boolean isCompatible() {
        getClades(targetTree.getRoot());
        for (BitSet clade: constraintsClades) {
            if (!targetClades.contains(clade)) {
                return false;
            }
        }

        return true;
    }

    /**
     * Compiles the set of clades from the constraints tree using the numbers from the target tree.
     * @param constraintsTree
     * @param node
     * @param targetTree
     * @return
     */
    private BitSet setupClades(Tree constraintsTree, NodeRef node, Tree targetTree) {

            BitSet clade = new BitSet();
            if (constraintsTree.isExternal(node)) {
                    String taxonId = constraintsTree.getNodeTaxon(node).getId();
                    clade.set(targetTree.getTaxonIndex(taxonId));
            } else {
                for (int i = 0; i < constraintsTree.getChildCount(node); i++) {
                    NodeRef child = constraintsTree.getChild(node, i);
                    clade.or(setupClades(constraintsTree, child, targetTree));
                }
                constraintsClades.add(clade);
            }
            return clade;
    }

    /**
     * Compiles the set of clades defined by bitsets on the tip numbers of the target Tree.
     * @param node - current node in traversal
     * @return BitSet - returned for recursion.
     */
    private BitSet getClades(NodeRef node) {
        int nodeIndex = node.getNumber();
        boolean isExternal = targetTree.isExternal(node);

        if(updateNode[nodeIndex]){
            BitSet clade = new BitSet();
            if (isExternal) {
                if (constrainedTips.contains(node.getNumber())) {
                    clade.set(node.getNumber());
                }
            } else {
                for (int i = 0; i < targetTree.getChildCount(node); i++) {
                    NodeRef child = targetTree.getChild(node, i);
                    clade.or(getClades(child));
                }
            }
            updateNode[nodeIndex]=false;
            targetTreeNodeCladeMap[nodeIndex]=clade;
            targetClades.add(clade);

        }
            return targetTreeNodeCladeMap[nodeIndex];
    }

    private final Tree targetTree;
    private final Set<BitSet> constraintsClades = new HashSet<>();
    private final Set<Integer> constrainedTips = new HashSet<>();
    private boolean likelihoodKnown = false;
    private boolean storedLikelihoodKnown=false;
    private double logLikelihood = 0.0;
    private double storedLogLikelihood = 0.0;
    private boolean[] updateNode;
    private List<Integer> lastUpdated;
    private Set<BitSet> targetClades = new HashSet<>();

    // With targetClade Set (above) we are storing the clades of the target tree twice
    // I don't see a way around this. the Set makes the comparison with the constraintsClades fast and
    // the array means we can keep a map between the node and it's clade.
    private  BitSet[] targetTreeNodeCladeMap;


}
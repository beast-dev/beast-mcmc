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
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.BitSet;
import java.util.HashSet;
import java.util.Set;

/**
 * Tests whether a tree is compatible with a less-resolved constraints tree
 * *
 *
 * @author Andrew Rambaut
 */
public class ConstraintsTreeLikelihood extends AbstractModelLikelihood {

    public ConstraintsTreeLikelihood(String name, Tree targetTree, Tree constraintsTree) throws TreeUtils.MissingTaxonException {

        super(name);

        for (int i = 0; i < targetTree.getTaxonCount(); i++) {
            String id = targetTree.getTaxonId(i);
            if (constraintsTree.getTaxonIndex(id) == -1) {
                throw new TreeUtils.MissingTaxonException(targetTree.getTaxon(i));
            }
        }

        getClades(constraintsTree, constraintsTree.getRoot(), targetTree, constraintsClades);

        if (targetTree instanceof TreeModel) {
            addModel((TreeModel)targetTree);
        }
        this.targetTree = targetTree;
    }


    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        makeDirty();
    }

    @Override
    protected void storeState() {
        assert(likelihoodKnown) : "Likelihood not known before a store";
        storedLogLikelihood = logLikelihood;
    }

    @Override
    protected void restoreState() {
        logLikelihood = storedLogLikelihood;
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
    }

    /**
     * Returns true if all the clades in the constraints tree are present in the target tree
     * @return
     */
    private boolean isCompatible() {
        Set<BitSet> targetClades = new HashSet<>();
        getClades(targetTree, targetTree.getRoot(), null, targetClades);

        for (BitSet clade: constraintsClades) {
            if (!targetClades.contains(clade)) {
                return false;
            }
        }

        return true;
    }

    /**
     * Compiles a set of clades defined by bitsets on the tip numbers. If reference tree is
     * given then it uses it to make sure the numbers refer to the tips in the reference tree.
     * @param tree
     * @param node
     * @param referenceTree
     * @param clades
     * @return
     */
    private BitSet getClades(Tree tree, NodeRef node, Tree referenceTree, Set<BitSet> clades) {
        BitSet clade = new BitSet();

        if (tree.isExternal(node)) {
            if (referenceTree != null) {
                String taxonId = tree.getNodeTaxon(node).getId();
                clade.set(referenceTree.getTaxonIndex(taxonId));
            }else{
                clade.set(node.getNumber());
            }


        } else {
            for (int i = 0; i < tree.getChildCount(node); i++) {
                NodeRef child = tree.getChild(node, i);
                clade.or(getClades(tree, child, referenceTree, clades));
            }

            clades.add(clade);
        }
        return clade;
    }

    private final Tree targetTree;
    private final Set<BitSet> constraintsClades = new HashSet<>();
    private boolean likelihoodKnown = false;
    private double logLikelihood = 0.0;
    private double storedLogLikelihood = 0.0;
}
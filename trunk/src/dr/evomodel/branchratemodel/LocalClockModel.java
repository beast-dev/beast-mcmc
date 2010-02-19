/*
 * LocalClockModel.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.branchratemodel.LocalClockModelParser;
import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

import java.util.BitSet;
import java.util.HashMap;
import java.util.Map;

/**
 * @author Andrew Rambaut
 * @version $Id: LocalClockModel.java,v 1.1 2005/04/05 09:27:48 rambaut Exp $
 */
public class LocalClockModel extends AbstractModel implements BranchRateModel {

    private TreeModel treeModel;
    protected Map<Integer, LocalClock> localTipClocks = new HashMap<Integer, LocalClock>();
    protected Map<BitSet, LocalClock> localCladeClocks = new HashMap<BitSet, LocalClock>();
    private boolean updateNodeClocks = true;
    private Map<NodeRef, LocalClock> nodeClockMap = new HashMap<NodeRef, LocalClock>();
    private final Parameter globalRateParameter;

    public LocalClockModel(TreeModel treeModel, Parameter globalRateParameter) {

        super(LocalClockModelParser.LOCAL_CLOCK_MODEL);
        this.treeModel = treeModel;

        addModel(treeModel);

        this.globalRateParameter = globalRateParameter;
        addVariable(globalRateParameter);
    }

    public void addExternalBranchClock(TaxonList taxonList, Parameter rateParameter, boolean relative) throws Tree.MissingTaxonException {
        BitSet tips = getTipsForTaxa(treeModel, taxonList);
        LocalClock clock = new LocalClock(rateParameter, relative, tips);
        for (int i = tips.nextSetBit(0); i >= 0; i = tips.nextSetBit(i + 1)) {
            localTipClocks.put(i, clock);
        }
        addVariable(rateParameter);
    }

    public void addCladeClock(TaxonList taxonList, Parameter rateParameter, boolean relative, boolean includeStem, boolean excludeClade) throws Tree.MissingTaxonException {
        BitSet tips = getTipsForTaxa(treeModel, taxonList);
        LocalClock clock = new LocalClock(rateParameter, relative, tips, includeStem, excludeClade);
        localCladeClocks.put(tips, clock);
        addVariable(rateParameter);
    }

    /**
     * @param tree the tree
     * @param taxa the taxa
     * @return A bitset with the node numbers set.
     * @throws dr.evolution.tree.Tree.MissingTaxonException
     *          if a taxon in taxa is not contained in the tree
     */
    private BitSet getTipsForTaxa(Tree tree, TaxonList taxa) throws Tree.MissingTaxonException {

        BitSet tips = new BitSet();

        for (int i = 0; i < taxa.getTaxonCount(); i++) {

            Taxon taxon = taxa.getTaxon(i);
            boolean found = false;
            for (int j = 0; j < tree.getExternalNodeCount(); j++) {

                NodeRef node = tree.getExternalNode(j);
                if (tree.getNodeTaxon(node).getId().equals(taxon.getId())) {
                    tips.set(node.getNumber());
                    found = true;
                    break;
                }
            }

            if (!found) {
                throw new Tree.MissingTaxonException(taxon);
            }
        }

        return tips;
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        updateNodeClocks = true;
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        fireModelChanged();
    }

    protected void storeState() {
    }

    protected void restoreState() {
        updateNodeClocks = true;
    }

    protected void acceptState() {
    }


    // BranchRateModel implementation

    public double getBranchRate(Tree tree, NodeRef node) {

        if (tree.isRoot(node)) {
            throw new IllegalArgumentException("root node doesn't have a rate!");
        }

        if (updateNodeClocks) {
            nodeClockMap.clear();
            setupRateParameters(tree, tree.getRoot(), new BitSet());

            updateNodeClocks = false;
        }

        double rate = globalRateParameter.getParameterValue(0);

        LocalClock localClock = nodeClockMap.get(node);
        if (localClock != null) {
            if (localClock.isRelativeRate()) {
                rate *= localClock.getRateParameter().getParameterValue(0);
            } else {
                rate = localClock.getRateParameter().getParameterValue(0);
            }
        }

        return rate;
    }

    public String getBranchAttributeLabel() {
        return RATE;
    }

    public String getAttributeForBranch(Tree tree, NodeRef node) {
        return Double.toString(getBranchRate(tree, node));
    }

    private void setupRateParameters(Tree tree, NodeRef node, BitSet tips) {
        LocalClock clock = null;

        if (tree.isExternal(node)) {
            tips.set(node.getNumber());
            clock = localTipClocks.get(node.getNumber());
        } else {
            for (int i = 0; i < tree.getChildCount(node); i++) {
                NodeRef child = tree.getChild(node, i);
                BitSet childTips = new BitSet();
                setupRateParameters(tree, child, childTips);

                tips.or(childTips);
            }
            clock = localCladeClocks.get(tips);
        }

        if (clock != null) {
            setNodeClock(tree, node, clock, clock.includeStem(), clock.excludeClade());
        }
    }

    private void setNodeClock(Tree tree, NodeRef node, LocalClock localClock, boolean includeStem, boolean excludeClade) {

        if (!tree.isExternal(node) && !excludeClade) {
            for (int i = 0; i < tree.getChildCount(node); i++) {
                NodeRef child = tree.getChild(node, i);
                setNodeClock(tree, child, localClock, true, true);
            }
        }

        if (includeStem && !nodeClockMap.containsKey(node)) {
            nodeClockMap.put(node, localClock);
        }
    }

    private class LocalClock {
        LocalClock(Parameter rateParameter, boolean relativeRate, BitSet tips) {
            this.rateParameter = rateParameter;
            this.relativeRate = relativeRate;
            this.tips = tips;
            this.isClade = false;
            this.includeStem = true;
            this.excludeClade = true;
        }

        LocalClock(Parameter rateParameter, boolean relativeRate, BitSet tips, boolean includeStem, boolean excludeClade) {
            this.rateParameter = rateParameter;
            this.relativeRate = relativeRate;
            this.tips = tips;
            this.isClade = true;
            this.includeStem = includeStem;
            this.excludeClade = excludeClade;
        }

        boolean includeStem() {
            return this.includeStem;
        }

        boolean excludeClade() {
            return excludeClade;
        }

        boolean isClade() {
            return this.isClade;
        }

        boolean isRelativeRate() {
            return relativeRate;
        }

        Parameter getRateParameter() {
            return this.rateParameter;
        }

        private final Parameter rateParameter;
        private final boolean relativeRate;
        private final BitSet tips;
        private final boolean isClade;
        private final boolean includeStem;
        private final boolean excludeClade;
    }

}
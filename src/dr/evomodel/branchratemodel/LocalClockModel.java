/*
 * LocalClockModel.java
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

package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.branchratemodel.LocalClockModelParser;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.*;

/**
 * @author Andrew Rambaut
 * @version $Id: LocalClockModel.java,v 1.1 2005/04/05 09:27:48 rambaut Exp $
 */
public class LocalClockModel extends AbstractBranchRateModel implements Citable {

    private TreeModel treeModel;
    protected Map<Integer, LocalClock> localTipClocks = new HashMap<Integer, LocalClock>();
    protected Map<BitSet, LocalClock> localCladeClocks = new HashMap<BitSet, LocalClock>();
    protected LocalClock trunkClock = null;

    private boolean updateNodeClocks = true;
    private Map<NodeRef, LocalClock> nodeClockMap = new HashMap<NodeRef, LocalClock>();
    private final Parameter globalRateParameter;

    public LocalClockModel(TreeModel treeModel, Parameter globalRateParameter) {

        super(LocalClockModelParser.LOCAL_CLOCK_MODEL);
        this.treeModel = treeModel;

        addModel(treeModel);

        this.globalRateParameter = globalRateParameter;
        addVariable(globalRateParameter);

        // add the super class' tree traits (just the rate)
        helper.addTrait(this);

        updateNodeClocks = true;

    }

    public void addExternalBranchClock(TaxonList taxonList, Parameter rateParameter, boolean isRelativeRate) throws TreeUtils.MissingTaxonException {
        Set<Integer> tips = TreeUtils.getTipsForTaxa(treeModel, taxonList);
        LocalClock clock = new LocalClock(rateParameter, isRelativeRate, tips, ClockType.EXTERNAL);
        for (int i : tips) {
            localTipClocks.put(i, clock);
        }
        addVariable(rateParameter);
    }

    public void addCladeClock(TaxonList taxonList, Parameter rateParameter, boolean isRelativeRate, double stemProportion, boolean excludeClade) throws TreeUtils.MissingTaxonException {
        Set<Integer> tips = TreeUtils.getTipsForTaxa(treeModel, taxonList);
        BitSet tipBitSet = TreeUtils.getTipsBitSetForTaxa(treeModel, taxonList);
        LocalClock clock = new LocalClock(rateParameter, isRelativeRate, tips, stemProportion, excludeClade);
        localCladeClocks.put(tipBitSet, clock);
        addVariable(rateParameter);
    }

    public void addTrunkClock(TaxonList taxonList, Parameter rateParameter, Parameter indexParameter, boolean isRelativeRate) throws TreeUtils.MissingTaxonException {
        if (trunkClock != null) {
            throw new RuntimeException("Trunk already defined for this LocalClockModel");
        }

        List<Integer> tipList = new ArrayList<Integer>(TreeUtils.getTipsForTaxa(treeModel, taxonList));
        trunkClock = new LocalClock(rateParameter, indexParameter, isRelativeRate, tipList, ClockType.TRUNK);
        addVariable(rateParameter);
        if (indexParameter != null) {
            addVariable(indexParameter);
        }

        helper.addTrait("trunk", new TreeTrait.S() {
            //            @Override
            public String getTraitName() {
                return "trunk";
            }

            //            @Override
            public Intent getIntent() {
                return Intent.BRANCH;
            }

            //            @Override
            public String getTrait(Tree tree, NodeRef node) {
                setupNodeClocks(tree);
                if (nodeClockMap.get(node) == trunkClock) {
                    return "T";
                }
                return "B";
            }
        });
    }

    public void handleModelChangedEvent(Model model, Object object, int index) {
        updateNodeClocks = true;
        fireModelChanged();
    }

    protected final void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (trunkClock != null && variable == trunkClock.indexParameter) {
            updateNodeClocks = true;
        }
        fireModelChanged();
    }

    protected void storeState() {
    }

    protected void restoreState() {
        updateNodeClocks = true;
    }

    protected void acceptState() {
    }

    // TreeTraitProvider overrides

    public TreeTrait[] getTreeTraits() {
        return helper.getTreeTraits();
    }

    public TreeTrait getTreeTrait(String key) {
        return helper.getTreeTrait(key);
    }

    // BranchRateModel implementation

    public double getBranchRate(final Tree tree, final NodeRef node) {

        if (tree.isRoot(node)) {
            throw new IllegalArgumentException("root node doesn't have a rate!");
        }

        setupNodeClocks(tree);

        double rate = globalRateParameter.getParameterValue(0);

        LocalClock parentClock = nodeClockMap.get(tree.getParent(node));
        LocalClock localClock = nodeClockMap.get(node);
        if (localClock != null) {
            double parentRate = rate;
            double stemProportion = 1.0;

            if (localClock != parentClock) {
                // this is the branch where the rate switch occurs
                if (parentClock != null) {
                    if (parentClock.isRelativeRate()) {
                        parentRate *= localClock.getRateParameter().getParameterValue(0);
                    } else {
                        parentRate = localClock.getRateParameter().getParameterValue(0);
                    }
                }
                stemProportion = localClock.getStemProportion();
            }

            if (localClock.isRelativeRate()) {
                rate *= localClock.getRateParameter().getParameterValue(0);
            } else {
                rate = localClock.getRateParameter().getParameterValue(0);
            }

            rate = (rate * stemProportion) + (parentRate * (1.0 - stemProportion));
        }

        return rate;
    }

    /**
     * Set up the map from node to clock.
     * @param tree
     */
    private void setupNodeClocks(final Tree tree) {
        if (updateNodeClocks) {
            nodeClockMap.clear();
            setupRateParameters(tree, tree.getRoot(), new BitSet());

            if (trunkClock != null) {
                // backbone will overwrite other local clocks
                setupTrunkRates(tree, tree.getRoot());
            }

            updateNodeClocks = false;
        }

    }

    /**
     * Traverse the tree getting bitsets for each node based on the tips below. If this
     * bitset is in the localCladeClocks map then set that clock for all the nodes below.
     * Pre-order traversal so shallowest clades are set first.
     * @param tree
     * @param node
     * @param tips
     */
    private void setupRateParameters(Tree tree, NodeRef node, BitSet tips) {
        LocalClock clock;

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
            setNodeClock(tree, node, clock, clock.excludeClade());
        }
    }

    private boolean setupTrunkRates(Tree tree, NodeRef node) {
        LocalClock clock = null;

        if (tree.isExternal(node)) {
            if (trunkClock.indexParameter != null) {
                if (trunkClock.tipList.get((int) trunkClock.indexParameter.getParameterValue(0)) == node.getNumber()) {
                    clock = trunkClock;
                }
            } else if (trunkClock.tipList.contains(node.getNumber())) {
                clock = trunkClock;
            }
        } else {
            for (int i = 0; i < tree.getChildCount(node); i++) {
                NodeRef child = tree.getChild(node, i);
                if (setupTrunkRates(tree, child)) {
                    // if any of the desendents are back bone then this node is too
                    clock = trunkClock;
                }
            }
        }

        if (clock != null) {
            setNodeClock(tree, node, clock, clock.excludeClade());
            return true;
        }

        return false;
    }

    /**
     * Traverse down the clade, associating all the nodes with the specified clock.
     * @param tree
     * @param node
     * @param localClock
     * @param excludeClade
     */
    private void setNodeClock(Tree tree, NodeRef node, LocalClock localClock, boolean excludeClade) {

        if (!tree.isExternal(node) && !excludeClade) {
            for (int i = 0; i < tree.getChildCount(node); i++) {
                NodeRef child = tree.getChild(node, i);
                setNodeClock(tree, child, localClock, false);
            }
        }

        if (!nodeClockMap.containsKey(node)) {
            nodeClockMap.put(node, localClock);
        }
    }

    enum ClockType {
        CLADE,
        TRUNK,
        EXTERNAL
    }

    private class LocalClock {

        LocalClock(Parameter rateParameter, boolean isRelativeRate, Set<Integer> tipSet, ClockType type) {
            this.rateParameter = rateParameter;
            this.indexParameter = null;
            this.isRelativeRate = isRelativeRate;
            this.tips = tipSet;
            this.tipList = null;
            this.type = type;
            this.stemProportion = 1.0;
            this.excludeClade = true;
        }

        LocalClock(Parameter rateParameter, Parameter indexParameter, boolean isRelativeRate, List<Integer> tipList, ClockType type) {
            this.rateParameter = rateParameter;
            this.indexParameter = indexParameter;
            this.isRelativeRate = isRelativeRate;
            this.tips = null;
            this.tipList = tipList;
            this.type = type;
            this.stemProportion = 1.0;
            this.excludeClade = true;
        }

        LocalClock(Parameter rateParameter, boolean isRelativeRate, Set<Integer> tips, double stemProportion, boolean excludeClade) {
            this.rateParameter = rateParameter;
            this.indexParameter = null;
            this.isRelativeRate = isRelativeRate;
            this.tips = tips;
            this.tipList = null;
            this.type = ClockType.CLADE;
            this.stemProportion = stemProportion;
            this.excludeClade = excludeClade;
        }

        double getStemProportion() {
            return this.stemProportion;
        }

        boolean excludeClade() {
            return excludeClade;
        }

        ClockType getType() {
            return this.type;
        }

        boolean isRelativeRate() {
            return isRelativeRate;
        }

        Parameter getRateParameter() {
            return this.rateParameter;
        }

        private final Parameter rateParameter;
        private final Parameter indexParameter;
        private final boolean isRelativeRate;
        private final Set<Integer> tips;
        private final List<Integer> tipList;
        private final ClockType type;
        private final double stemProportion;
        private final boolean excludeClade;
    }

    private final Helper helper = new Helper();

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.MOLECULAR_CLOCK;
    }

    @Override
    public String getDescription() {
        return "Local clock model";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CITATION);
    }

    public static Citation CITATION = new Citation(
            new Author[]{
                    new Author("AD", "Yoder"),
                    new Author("Z", "Yang")
            },
            "Estimation of Primate Speciation Dates Using Local Molecular Clocks",
            2000,
            "Mol Biol Evol",
            17, 1081, 1090
    );
}
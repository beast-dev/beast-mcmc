/*
 * LocalClockModel.java
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

package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.TaxonList;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.branchratemodel.LocalClockModelParser;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.*;

/**
 * @author Andrew Rambaut
 */
public class LocalClockModel extends AbstractBranchRateModel implements Citable, DifferentiableBranchRates {

    private TreeModel treeModel;
    protected Map<Integer, LocalClock> localTipClocks = new HashMap<Integer, LocalClock>();
    protected Map<BitSet, LocalClock> localCladeClocks = new HashMap<BitSet, LocalClock>();
    protected LocalClock trunkClock = null;

    private boolean updateNodeClocks = true;
    private Map<NodeRef, LocalClock> nodeClockMap = new HashMap<NodeRef, LocalClock>();
    private final Parameter globalRateParameter;
    private final BranchRateModel globalBranchRates;

    public LocalClockModel(TreeModel treeModel, Parameter globalRateParameter) {

        super(LocalClockModelParser.LOCAL_CLOCK_MODEL);
        this.treeModel = treeModel;

        addModel(treeModel);

        this.globalRateParameter = globalRateParameter;
        this.globalBranchRates = null;

        addVariable(globalRateParameter);

        // add the super class' tree traits (just the rate)
        helper.addTrait(this);

        updateNodeClocks = true;

    }

    public LocalClockModel(TreeModel treeModel, BranchRateModel globalBranchRates) {

        super(LocalClockModelParser.LOCAL_CLOCK_MODEL);
        this.treeModel = treeModel;

        addModel(treeModel);

        this.globalRateParameter = null;
        this.globalBranchRates = globalBranchRates;

        addModel(globalBranchRates);

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
        clockOrder.add(clock);
        addVariable(rateParameter);
    }

    public void addExternalBranchClock(TaxonList taxonList, BranchRateModel branchRates, boolean isRelativeRate) throws TreeUtils.MissingTaxonException {
        Set<Integer> tips = TreeUtils.getTipsForTaxa(treeModel, taxonList);
        LocalClock clock = new LocalClock(branchRates, isRelativeRate, tips, ClockType.EXTERNAL);
        for (int i : tips) {
            localTipClocks.put(i, clock);
        }
        clockOrder.add(clock);
        addModel(branchRates);
    }

    public void addCladeClock(TaxonList taxonList, Parameter rateParameter, boolean isRelativeRate, Parameter stemParameter, boolean stemAsTime, boolean excludeClade) throws TreeUtils.MissingTaxonException {
        Set<Integer> tips = TreeUtils.getTipsForTaxa(treeModel, taxonList);
        BitSet tipBitSet = TreeUtils.getTipsBitSetForTaxa(treeModel, taxonList);
        LocalClock clock = new LocalClock(rateParameter, isRelativeRate, tips, stemParameter, stemAsTime, excludeClade);
        localCladeClocks.put(tipBitSet, clock);
        clockOrder.add(clock);
        addVariable(rateParameter);
        if (stemParameter != null) {
            addVariable(stemParameter);
        }
    }

    public void addCladeClock(TaxonList taxonList, BranchRateModel branchRates, boolean isRelativeRate, Parameter stemParameter, boolean stemAsTime, boolean excludeClade) throws TreeUtils.MissingTaxonException {
        Set<Integer> tips = TreeUtils.getTipsForTaxa(treeModel, taxonList);
        BitSet tipBitSet = TreeUtils.getTipsBitSetForTaxa(treeModel, taxonList);
        LocalClock clock = new LocalClock(branchRates, isRelativeRate, tips, stemParameter, stemAsTime, excludeClade);
        localCladeClocks.put(tipBitSet, clock);
        clockOrder.add(clock);
        addModel(branchRates);
        if (stemParameter != null) {
            addVariable(stemParameter);
        }
    }

    public void addTrunkClock(TaxonList taxonList, Parameter rateParameter, Parameter indexParameter, boolean isRelativeRate) throws TreeUtils.MissingTaxonException {
        if (trunkClock != null) {
            throw new RuntimeException("Trunk already defined for this LocalClockModel");
        }

        List<Integer> tipList = new ArrayList<Integer>(TreeUtils.getTipsForTaxa(treeModel, taxonList));
        trunkClock = new LocalClock(rateParameter, indexParameter, isRelativeRate, tipList, ClockType.TRUNK);
        clockOrder.add(trunkClock);
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

    public void addTrunkClock(TaxonList taxonList, BranchRateModel branchRates, Parameter indexParameter, boolean isRelativeRate) throws TreeUtils.MissingTaxonException {
        if (trunkClock != null) {
            throw new RuntimeException("Trunk already defined for this LocalClockModel");
        }

        List<Integer> tipList = new ArrayList<Integer>(TreeUtils.getTipsForTaxa(treeModel, taxonList));
        trunkClock = new LocalClock(branchRates, indexParameter, isRelativeRate, tipList, ClockType.TRUNK);
        clockOrder.add(trunkClock);
        addModel(branchRates);
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

        double rate = (globalRateParameter != null ?
                globalRateParameter.getParameterValue(0) :
                globalBranchRates.getBranchRate(tree, node));

        LocalClock parentClock = nodeClockMap.get(tree.getParent(node));
        LocalClock localClock = nodeClockMap.get(node);
        if (localClock != null) {
            double parentRate = rate;
            double stemProportion = 1.0;
            double stemTime = 0.0;

            if (localClock != parentClock) {
                // this is the branch where the rate switch occurs
                if (parentClock != null) {
                    if (parentClock.isRelativeRate()) {
                        parentRate *= localClock.getBranchRate(tree, tree.getParent(node));
                    } else {
                        parentRate = localClock.getBranchRate(tree, tree.getParent(node));
                    }
                }
                if (localClock.stemAsTime) {
                    // this could be greater than 1 in which case bad things might happen
                    stemTime = localClock.getStemValue();
                    stemProportion = stemTime / tree.getBranchLength(node);

                    if (stemProportion > 1.0) {
                        // it should be ensured that this never happens.
                        throw new IllegalArgumentException("A stem proportion for a local clock is > 1.0");
                    }
                } else {
                    stemProportion = localClock.getStemValue();
                    stemTime = tree.getBranchLength(node) * stemProportion;
                }
            }

            if (localClock.isRelativeRate()) {
                rate *= localClock.getBranchRate(tree, node);
            } else {
                rate = localClock.getBranchRate(tree, node);
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

    // DifferentiableBranchRates: the gradient is wrt the global rate (if a parameter) and the rate parameter of each
    // local clock, stacked in that order; stems and the trunk index are held fixed

    @Override
    public double getBranchRateDifferential(Tree tree, NodeRef node) {
        return 1.0; // the per-branch gradient is wrt the branch rate itself; updateGradientLogDensity chains it
    }

    @Override
    public double getBranchRateSecondDifferential(Tree tree, NodeRef node) {
        return 0.0;
    }

    @Override
    public Parameter getRateParameter() {
        if (rateParameters == null) {
            List<Parameter> parameters = new ArrayList<Parameter>();
            if (globalRateParameter == null) {
                throw new UnsupportedOperationException("Local clock gradients need a global rate parameter, " +
                        "not a branch rate model");
            }
            parameters.add(globalRateParameter);
            for (LocalClock clock : clockOrder) {
                if (clock.rateParameter == null) {
                    throw new UnsupportedOperationException("Local clock gradients need rate parameters, " +
                            "not branch rate models");
                }
                if (!parameters.contains(clock.rateParameter)) {
                    parameters.add(clock.rateParameter);
                }
            }
            rateParameters = parameters;
            compoundRateParameter = new CompoundParameter("localClockRates",
                    parameters.toArray(new Parameter[0]));
        }
        return compoundRateParameter;
    }

    @Override
    public int getParameterIndexFromNode(NodeRef node) {
        throw new UnsupportedOperationException("A local clock rate is shared by many branches");
    }

    @Override
    public ArbitraryBranchRates.BranchRateTransform getTransform() {
        throw new RuntimeException("Not yet implemented!");
    }

    /**
     * @param gradient d log L / d rate_e, indexed by branch (node number with the root skipped)
     * @return d log L / d theta for theta in getRateParameter()
     */
    @Override
    public double[] updateGradientLogDensity(double[] gradient, double[] value, int from, int to) {
        getRateParameter();
        setupNodeClocks(treeModel);
        double[] result = new double[rateParameters.size()];
        int v = 0;
        for (int i = 0; i < treeModel.getNodeCount(); ++i) {
            final NodeRef node = treeModel.getNode(i);
            if (!treeModel.isRoot(node)) {
                addBranchRateDerivatives(treeModel, node, gradient[v], result);
                v++;
            }
        }
        return result;
    }

    @Override
    public double[] updateDiagonalHessianLogDensity(double[] diagonalHessian, double[] gradient, double[] value, int from, int to) {
        throw new RuntimeException("Not yet implemented!");
    }

    /*
     * Adds weight * d rate / d theta to result, following getBranchRate:
     * rate = s * R + (1 - s) * P, with R the clock's rate (relative: g * r, else r), P the parent-side rate and s the
     * stem proportion (1 off the switching branch).
     */
    private void addBranchRateDerivatives(Tree tree, NodeRef node, double weight, double[] result) {
        final int global = 0;
        final double g = globalRateParameter.getParameterValue(0);

        final LocalClock localClock = nodeClockMap.get(node);
        if (localClock == null) {
            result[global] += weight;
            return;
        }
        final LocalClock parentClock = nodeClockMap.get(tree.getParent(node));
        final int local = rateParameters.indexOf(localClock.rateParameter);
        final double r = localClock.rateParameter.getParameterValue(0);

        double stemProportion = 1.0;
        // the parent-side rate as getBranchRate computes it (it reads the local clock's rate)
        double dParentRateDg = 1.0;
        double dParentRateDr = 0.0;
        if (localClock != parentClock) {
            if (parentClock != null) {
                if (parentClock.isRelativeRate()) {
                    dParentRateDg = r;
                    dParentRateDr = g;
                } else {
                    dParentRateDg = 0.0;
                    dParentRateDr = 1.0;
                }
            }
            if (localClock.stemAsTime) {
                stemProportion = localClock.getStemValue() / tree.getBranchLength(node);
            } else {
                stemProportion = localClock.getStemValue();
            }
        }

        final double dRateDg = localClock.isRelativeRate() ? r : 0.0;
        final double dRateDr = localClock.isRelativeRate() ? g : 1.0;

        result[global] += weight * (stemProportion * dRateDg + (1.0 - stemProportion) * dParentRateDg);
        result[local] += weight * (stemProportion * dRateDr + (1.0 - stemProportion) * dParentRateDr);
    }

    private final List<LocalClock> clockOrder = new ArrayList<LocalClock>();
    private List<Parameter> rateParameters = null;
    private CompoundParameter compoundRateParameter = null;

    enum ClockType {
        CLADE,
        TRUNK,
        EXTERNAL
    }

    private class LocalClock {

        LocalClock(Parameter rateParameter, boolean isRelativeRate, Set<Integer> tipSet, ClockType type) {
            this.rateParameter = rateParameter;
            this.branchRates = null;
            this.indexParameter = null;
            this.isRelativeRate = isRelativeRate;
            this.tips = tipSet;
            this.tipList = null;
            this.type = type;
            this.stemParameter = null;
            this.stemAsTime = false;
            this.excludeClade = true;
        }

        LocalClock(BranchRateModel branchRates, boolean isRelativeRate, Set<Integer> tipSet, ClockType type) {
            this.rateParameter = null;
            this.branchRates = branchRates;
            this.indexParameter = null;
            this.isRelativeRate = isRelativeRate;
            this.tips = tipSet;
            this.tipList = null;
            this.type = type;
            this.stemParameter = null;
            this.stemAsTime = false;
            this.excludeClade = true;
        }

        LocalClock(Parameter rateParameter, boolean isRelativeRate, Set<Integer> tips, Parameter stemParameter, boolean stemAsTime, boolean excludeClade) {
            this.rateParameter = rateParameter;
            this.branchRates = null;
            this.indexParameter = null;
            this.isRelativeRate = isRelativeRate;
            this.tips = tips;
            this.tipList = null;
            this.type = ClockType.CLADE;
            this.stemParameter = stemParameter;
            this.stemAsTime = stemAsTime;
            this.excludeClade = excludeClade;
        }

        LocalClock(BranchRateModel branchRates, boolean isRelativeRate, Set<Integer> tips, Parameter stemParameter, boolean stemAsTime, boolean excludeClade) {
            this.rateParameter = null;
            this.branchRates = branchRates;
            this.indexParameter = null;
            this.isRelativeRate = isRelativeRate;
            this.tips = tips;
            this.tipList = null;
            this.type = ClockType.CLADE;
            this.stemParameter = stemParameter;
            this.stemAsTime = stemAsTime;
            this.excludeClade = excludeClade;
        }

        LocalClock(Parameter rateParameter, Parameter indexParameter, boolean isRelativeRate, List<Integer> tipList, ClockType type) {
            this.rateParameter = rateParameter;
            this.branchRates = null;
            this.indexParameter = indexParameter;
            this.isRelativeRate = isRelativeRate;
            this.tips = null;
            this.tipList = tipList;
            this.type = type;
            this.stemParameter = null;
            this.stemAsTime = false;
            this.excludeClade = true;
        }

        LocalClock(BranchRateModel branchRates, Parameter indexParameter, boolean isRelativeRate, List<Integer> tipList, ClockType type) {
            this.rateParameter = null;
            this.branchRates = branchRates;
            this.indexParameter = indexParameter;
            this.isRelativeRate = isRelativeRate;
            this.tips = null;
            this.tipList = tipList;
            this.type = type;
            this.stemParameter = null;
            this.stemAsTime = false;
            this.excludeClade = true;
        }

        double getStemValue() {
            if (stemParameter != null) {
                return stemParameter.getParameterValue(0);
            } else {
                // if no parameter then default to 0 stem 
                return 0.0;
            }
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

        double getBranchRate(Tree tree, NodeRef node) {
            if (rateParameter != null) {
                return this.rateParameter.getParameterValue(0);
            } else {
                return this.branchRates.getBranchRate(tree, node);
            }
        }

        private final Parameter rateParameter;
        private final BranchRateModel branchRates;
        private final Parameter indexParameter;
        private final boolean isRelativeRate;
        private final Set<Integer> tips;
        private final List<Integer> tipList;
        private final ClockType type;
        private final Parameter stemParameter;
        private final boolean stemAsTime;
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
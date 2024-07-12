/*
 * BastaLikelihood.java
 *
 * Copyright (c) 2002-2023 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.evomodel.coalescent.basta;

import dr.evolution.alignment.PatternList;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.bigfasttree.BestSignalsFromBigFastTreeIntervals;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.*;
import dr.util.Citable;
import dr.util.Citation;
import dr.xml.Reportable;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

import static dr.evomodel.coalescent.basta.ProcessOnCoalescentIntervalDelegate.*;

/**
 * @author Guy Baele
 * @author Yucai Shao
 * @author Marc A. Suchard
 */

public class BastaLikelihood extends AbstractModelLikelihood implements
        TreeTraitProvider, Citable, Profileable, Reportable {

    private static final boolean COUNT_TOTAL_OPERATIONS = true;

    private final BastaLikelihoodDelegate likelihoodDelegate;

    private final Tree tree;
    private final PatternList patternList;
    private final SubstitutionModel substitutionModel;
    private final Parameter popSizeParameter;
    private final BranchRateModel branchRateModel;
    private final int stateCount;

    private final Helper treeTraits = new Helper();

    private final CoalescentIntervalTraversal treeTraversalDelegate;
    private final BestSignalsFromBigFastTreeIntervals treeIntervals;

    private double logLikelihood;
    private double storedLogLikelihood;
    protected boolean likelihoodKnown;

    private boolean populationSizesKnown;
    private boolean treeIntervalsKnown;
    private boolean transitionMatricesKnown;

    public BastaLikelihood(String name,
                           Tree treeModel,
                           PatternList patternList,
                           SubstitutionModel substitutionModel,
                           Parameter popSizeParameter,
                           BranchRateModel branchRateModel,
                           BastaLikelihoodDelegate likelihoodDelegate,
                           int numberSubIntervals,
                           boolean useAmbiguities) {

        super(name);

        assert likelihoodDelegate != null;
        assert treeModel != null;
        assert branchRateModel != null;
        assert patternList.getPatternCount() == 1;
        assert useAmbiguities;

        if (!(branchRateModel instanceof StrictClockBranchRates)) {
            throw new RuntimeException("Not yet implemented");
        }

        final Logger logger = Logger.getLogger("dr.evomodel");

        logger.info("\nUsing BastaLikelihood");

        this.patternList = patternList;

        this.likelihoodDelegate = likelihoodDelegate;
        addModel(likelihoodDelegate);

        this.tree = treeModel;

        this.branchRateModel = branchRateModel;
        addModel(branchRateModel);

        this.substitutionModel = substitutionModel;
        addModel(substitutionModel);

        this.popSizeParameter = popSizeParameter;
        addVariable(popSizeParameter);

        this.stateCount = substitutionModel.getDataType().getStateCount();

        if (tree instanceof TreeModel) {
            treeIntervals = new BestSignalsFromBigFastTreeIntervals((TreeModel) treeModel);
            addModel(treeIntervals);
        } else {
            throw new RuntimeException("Not yet implemented");
        }

        treeTraversalDelegate = new CoalescentIntervalTraversal(treeModel, treeIntervals, branchRateModel, numberSubIntervals);

        setTipData();

        likelihoodKnown = false;
        populationSizesKnown = false;
        treeIntervalsKnown = false;
        transitionMatricesKnown = false;
    }

    public CoalescentIntervalTraversal getTraversalDelegate() { return treeTraversalDelegate; }

    public SubstitutionModel getSubstitutionModel() { return substitutionModel; } // TODO generify for multiple models (e.g. epochs)

    private void setTipData() {

        int[] data = patternList.getPattern(0);

        for (int i = 0; i < tree.getExternalNodeCount(); ++i) {
            NodeRef node = tree.getExternalNode(i);

            int index = patternList.getTaxonIndex(tree.getNodeTaxon(node).getId());
            int datum = data[index];

            if (datum >= stateCount) {
                throw new RuntimeException("Not yet implemented");
            }

            double[] partials = new double[stateCount];
            partials[datum] = 1.0;

            likelihoodDelegate.setPartials(node.getNumber(), partials);
        }
    }

    @Override
    public final Model getModel() {
        return this;
    }

    @Override @SuppressWarnings("Duplicates")
    public final double getLogLikelihood() {
        if (COUNT_TOTAL_OPERATIONS) totalGetLogLikelihoodCount++;

        if (!likelihoodKnown) {
            long startTime;
            if (COUNT_TOTAL_OPERATIONS) {
                totalCalculateLikelihoodCount++;
                startTime = System.nanoTime();
            }

            logLikelihood = calculateLogLikelihood();

            if (COUNT_TOTAL_OPERATIONS) {
                long endTime = System.nanoTime();
                totalLikelihoodTime += (endTime - startTime) / 1000000;
            }

            likelihoodKnown = true;
        }

        return logLikelihood;
    }

    @Override
    public final void makeDirty() {
        if (COUNT_TOTAL_OPERATIONS) totalMakeDirtyCount++;

        likelihoodKnown = false;
        treeIntervalsKnown = false;
        populationSizesKnown = false;
        transitionMatricesKnown = false;

        likelihoodDelegate.makeDirty();
        updateAllNodes();
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (variable == popSizeParameter) {
            populationSizesKnown = false;
            likelihoodKnown = false;
        } else {
            throw new RuntimeException("Not yet implemented");
        }
    }

    @Override
    protected final void handleModelChangedEvent(Model model, Object object, int index) {

        if (model == treeIntervals) {
            treeIntervalsKnown = false;
            transitionMatricesKnown = false;
        } else if (model == branchRateModel) {
            treeIntervalsKnown = false; // TODO should not be necessary
            transitionMatricesKnown = false;
        } else if (model == substitutionModel) {
            transitionMatricesKnown = false;
        } else {
            throw new RuntimeException("Not yet implemented");
        }

        if (COUNT_TOTAL_OPERATIONS) totalModelChangedCount++;

        likelihoodKnown = false;
        fireModelChanged();
    }

    @Override
    protected final void storeState() {
        assert (likelihoodKnown) : "the likelihood should always be known at this point in the cycle";
        assert (populationSizesKnown);
        assert (treeIntervalsKnown);
        assert (transitionMatricesKnown);

        storedLogLikelihood = logLikelihood;
    }

    @Override
    protected final void restoreState() {
        logLikelihood = storedLogLikelihood;

        likelihoodKnown = true;
        populationSizesKnown = true;
        treeIntervalsKnown = true;
        transitionMatricesKnown = true;
    }

    @Override
    protected void acceptState() { } // nothing to do

    private double calculateLogLikelihood() {

        if (!transitionMatricesKnown) {
            // update eigen-decomposition
            likelihoodDelegate.updateEigenDecomposition(0, substitutionModel.getEigenDecomposition(), false); // TODO do conditionally and double-buffer
        }

        if (!populationSizesKnown) {
            // update population sizes
            likelihoodDelegate.updatePopulationSizes(0, popSizeParameter.getParameterValues(), false); // TODO do conditionally and double-buffer
        }

        if (!treeIntervalsKnown) {
            // update operations on tree
            treeTraversalDelegate.dispatchTreeTraversalCollectBranchAndNodeOperations();
        }

        final List<BranchIntervalOperation> branchOperations =
                treeTraversalDelegate.getBranchIntervalOperations();
        final List<TransitionMatrixOperation> matrixOperations =
                transitionMatricesKnown ? NO_OPT :
                treeTraversalDelegate.getMatrixOperations();
        final List<Integer> intervalStarts = treeTraversalDelegate.getIntervalStarts();

        if (COUNT_TOTAL_OPERATIONS) {
            totalPropagationCount += branchOperations.size();
            totalMatrixUpdateCount += matrixOperations.size();
            totalIntervalReductionCount += treeTraversalDelegate.getCoalescentIntervalCount();
        }

        final NodeRef root = tree.getRoot();
        double logL = likelihoodDelegate.calculateLikelihood(branchOperations, matrixOperations,
                intervalStarts, root.getNumber());

        // after traverse all nodes and patterns have been updated --
        //so change flags to reflect this.
        setAllNodesUpdated();

        treeIntervalsKnown = true;
        populationSizesKnown = true;
        transitionMatricesKnown = true;

        return logL;
    }

    public double[] getGradientLogDensity(StructuredCoalescentLikelihoodGradient wrt) {

        final List<BranchIntervalOperation> branchOperations =
                treeTraversalDelegate.getBranchIntervalOperations();
        final List<TransitionMatrixOperation> matrixOperations =
                transitionMatricesKnown ? NO_OPT :
                        treeTraversalDelegate.getMatrixOperations();
        final List<Integer> intervalStarts = treeTraversalDelegate.getIntervalStarts();

        final NodeRef root = tree.getRoot();

        calculateLogLikelihood(); // TODO Only execute if necessary

        double[] gradient = likelihoodDelegate.calculateGradientGeneric(branchOperations, matrixOperations, intervalStarts,
                root.getNumber(), wrt);

        return wrt.chainRule(gradient);
    }

    private void setAllNodesUpdated() {
        treeTraversalDelegate.setAllNodesUpdated();
    }

    /**
     * Set update flag for a node only
     */
    protected void updateNode(NodeRef node) {
        if (COUNT_TOTAL_OPERATIONS) totalRateUpdateSingleCount++;

        treeTraversalDelegate.updateNode(node);
        likelihoodKnown = false;
    }

    protected void updateAllNodes() {
        if (COUNT_TOTAL_OPERATIONS) totalRateUpdateAllCount++;

        treeTraversalDelegate.updateAllNodes();
        likelihoodKnown = false;
    }

    @Override
    public String getReport() {
        StringBuilder sb = new StringBuilder();

        double logL = getLogLikelihood();

        String delegateString = likelihoodDelegate.getReport();
        if (delegateString != null) {
            sb.append(delegateString);
        }

        sb.append(getClass().getName()).append("(").append(logL).append(")");

        if (COUNT_TOTAL_OPERATIONS)
            sb.append(
                    "\n  propagation operations = ").append(totalPropagationCount).append(
                    "\n  matrix updates = ").append(totalMatrixUpdateCount).append(
                    "\n  interval operations = ").append(totalIntervalReductionCount).append(
                    "\n  model changes = ").append(totalModelChangedCount).append(
                    "\n  make dirties = ").append(totalMakeDirtyCount).append(
                    "\n  calculate likelihoods = ").append(totalCalculateLikelihoodCount).append(
                    "\n  get likelihoods = ").append(totalGetLogLikelihoodCount).append(
                    "\n  all rate updates = ").append(totalRateUpdateAllCount).append(
                    "\n  partial rate updates = ").append(totalRateUpdateSingleCount).append(
                    "\n  average likelihood time = ").append(totalLikelihoodTime / totalCalculateLikelihoodCount);

        return sb.toString();
    }

    @Override
    public TreeTrait[] getTreeTraits() {
        return treeTraits.getTreeTraits();
    }

    @Override
    public TreeTrait getTreeTrait(String key) {
        return treeTraits.getTreeTrait(key);
    }

    public void addTrait(TreeTrait trait) {
        treeTraits.addTrait(trait);
    }

    public void addTraits(TreeTrait[] traits) {
        treeTraits.addTraits(traits);
    }

    @Override
    public Citation.Category getCategory() { return Citation.Category.TREE_PRIORS; }

    @Override
    public String getDescription() {
        if (likelihoodDelegate instanceof Citable) {
            return ((Citable)likelihoodDelegate).getDescription();
        } else {
            return null;
        }
    }

    public BastaLikelihoodDelegate getLikelihoodDelegate() {  return likelihoodDelegate; }

    @Override
    public List<Citation> getCitations() {
        if (likelihoodDelegate instanceof Citable) {
            return ((Citable)likelihoodDelegate).getCitations();
        } else {
            return new ArrayList<>();
        }
    }

    @Override
    public long getTotalCalculationCount() {
        return likelihoodDelegate.getTotalCalculationCount();
    }

    private final List<TransitionMatrixOperation> NO_OPT = new ArrayList<>();

    private int totalPropagationCount = 0;
    private int totalMatrixUpdateCount = 0;
    private int totalIntervalReductionCount = 0;
    private int totalGetLogLikelihoodCount = 0;
    private int totalModelChangedCount = 0;
    private int totalMakeDirtyCount = 0;
    private int totalCalculateLikelihoodCount = 0;
    private int totalRateUpdateAllCount = 0;
    private int totalRateUpdateSingleCount = 0;

    private long totalLikelihoodTime = 0;

    public Parameter getPopSizes() {return popSizeParameter;}
}

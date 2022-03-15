/*
 * StructuredCoalescentLikelihood.java
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

package dr.evomodel.coalescent.basta;

import dr.evolution.alignment.PatternList;
import dr.evolution.coalescent.IntervalList;
import dr.evolution.coalescent.IntervalType;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneralDataType;
import dr.evolution.tree.*;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evomodel.bigfasttree.BigFastTreeIntervals;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.coalescent.TreeIntervals;
import dr.evomodel.substmodel.GeneralSubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.tree.TreeChangedEvent;
import dr.evomodel.treelikelihood.AncestralStateTraitProvider;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.math.MathUtils;
import dr.util.*;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * @author Guy Baele
 *
 * Implementation of BASTA: Bayesian structured coalescent approximation.
 * Original paper: Nicola De Maio, Chieh-Hsi Wu, Kathleen O'Reilly and Daniel Wilson
 * "New routes to phylogeography: a Bayesian structured coalescent approximation".
 * PLOS Genetics 11, e1005421; doi: 10.1371/journal.pgen.1005421
 */
public class StructuredCoalescentLikelihood extends AbstractModelLikelihood implements Units, Citable {

    //TODO: the likelihood class should have minimum functionality, i.e. if the likelihood needs to be recomputed
    //then compute it by calling a StructuredCoalescentModel calculateLikelihood method; if not, then simply return it
    //TODO: create StructuredCoalescentModel class
    //TODO: have StructuredCoalescentLikelihood listen to both BeagleMatrixExponentiationDelegate and
    //the StructuredCoalescentModel classes?

    private static final boolean DEBUG = true;
    private static final boolean MATRIX_DEBUG = false;
    private static final boolean UPDATE_DEBUG = false;

    //private static final boolean USE_BEAGLE = false;
    private static final boolean ASSOC_MULTIPLICATION = true;
    private static final boolean USE_TRANSPOSE = false;

    private static final boolean MINIMUM_EVALUATION = false;

    public StructuredCoalescentLikelihood(Tree tree, BranchRateModel branchRateModel, Parameter popSizes, PatternList patternList,
                                          DataType dataType, String tag, GeneralSubstitutionModel generalSubstitutionModel, int subIntervals,
                                          TaxonList includeSubtree, List<TaxonList> excludeSubtrees, boolean useMAP) throws TreeUtils.MissingTaxonException {

        super(StructuredCoalescentLikelihoodParser.STRUCTURED_COALESCENT);

        this.treeModel = (TreeModel)tree;
        this.patternList = patternList;
        this.dataType = dataType;
        this.useMAP = useMAP;

        if (tree instanceof TreeModel) {
            addModel((TreeModel) tree);
        }

        if (tree instanceof TreeModel) {
            this.intervals = new BigFastTreeIntervals((TreeModel) tree);
        } else {
            throw new IllegalArgumentException("Please provide a TreeModel for the structured coalescent model.");
        }

        this.popSizes = popSizes;
        addVariable(this.popSizes);

        if (branchRateModel != null) {
            this.branchRateModel = branchRateModel;
        } else {
            this.branchRateModel = new DefaultBranchRateModel();
        }
        addModel(this.branchRateModel);

        this.generalSubstitutionModel = generalSubstitutionModel;
        addModel(this.generalSubstitutionModel);

        this.demes = generalSubstitutionModel.getDataType().getStateCount();
        this.startExpected = new double[this.demes];
        this.endExpected = new double[this.demes];

        int nodeCount = treeModel.getNodeCount();

        this.intervalStartProbs = new Double[intervals.getIntervalCount()+1][];
        this.intervalEndProbs = new Double[intervals.getIntervalCount()+1][];
        this.intervalStartSquareProbs = new Double[intervals.getIntervalCount()+1][];
        this.intervalEndSquareProbs = new Double[intervals.getIntervalCount()+1][];
        this.coalescentLeftProbs = new Double[intervals.getIntervalCount()+1][];
        this.coalescentRightProbs = new Double[intervals.getIntervalCount()+1][];
        //keep track of currently active lineage probabilities and copy to the lists above
        this.activeLineages = new ArrayList<Double[]>();
        //TODO figure out if I still need this
        //this.currentLineageProbs = new ArrayList<Double[]>();

        this.activeLineageList = new ProbDist[nodeCount];

        this.addedLineages = new boolean[nodeCount];
        for (int i = 0; i < addedLineages.length; i++) {
            addedLineages[i] = false;
        }
        this.addedLength = 0;

        this.nodeProbDist = new ProbDist[nodeCount];
        for (int i = 0; i < this.nodeProbDist.length; i++) {
            this.nodeProbDist[i] = new ProbDist(demes);
        }
        //set the pattern and starting lineage probability for each external node
        for (int i = 0; i < treeModel.getExternalNodeCount(); i++) {
            NodeRef refNode = treeModel.getExternalNode(i);
            this.nodeProbDist[refNode.getNumber()].patternIndex = patternList.getPattern(0)[patternList.getTaxonIndex(treeModel.getNodeTaxon(refNode).getId())];
            this.nodeProbDist[refNode.getNumber()].node = refNode;
            this.nodeProbDist[refNode.getNumber()].intervalType = IntervalType.SAMPLE;
            this.nodeProbDist[refNode.getNumber()].leftChild = null;
            this.nodeProbDist[refNode.getNumber()].rightChild = null;
        }

        this.maxCoalescentIntervals = treeModel.getTaxonCount() * 2 - 2;
        this.currentCoalescentInterval = 0;
        this.migrationMatrices = new double[maxCoalescentIntervals][this.demes*this.demes];
        this.storedMigrationMatrices = new double[maxCoalescentIntervals][this.demes*this.demes];

        this.times = new ArrayList<ComparableDouble>();
        this.children = new ArrayList<Integer>();
        this.nodes = new ArrayList<NodeRef>();

        this.storedTimes = new ArrayList<ComparableDouble>();
        this.storedChildren = new ArrayList<Integer>();
        this.storedNodes = new ArrayList<NodeRef>();

        this.treeModelUpdateFired = false;
        this.rateChanged = false;

        this.likelihoodKnown = false;

    }

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    public final Model getModel() {
        return this;
    }

    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
            //updateAllDensities(false);
        }
        return logLikelihood;
    }

    /**
     * Calculates the log structured coalescent density across all coalescent intervals.
     */
    public double calculateLogLikelihood() {

        //TODO write if-else clause that determines whether a full likelihood recalculation is required
        //TODO if not, write a method that only recomputes part of the probability distributions
        computeProbabilityDistributions();

        return calculateLogLikelihood(0);

    }

    /**
     * Calculates the log structured coalescent density from a given interval number up to the root.
     */
    public double calculateLogLikelihood(int startingInterval) {

        double logL = Double.NEGATIVE_INFINITY;

        //TODO only iterate over a subset of the intervals
        int intervalCount = intervals.getIntervalCount();
        for (int i = 0; i < intervalCount; i++) {
            double intervalLength = intervals.getIntervalTime(i);

            //TODO this could be funky for zero-length intervals
            if (intervalLength != 0.0) {
                //get all the information from the various lists
                Double[] lineageStartCount = intervalStartProbs[i];
                Double[] lineageStartCountSquare = intervalStartSquareProbs[i];
                Double[] lineageEndCount = intervalEndProbs[i];
                Double[] lineageEndCountSquare = intervalEndSquareProbs[i];

                double halfLength = intervalLength/2.0;
                if (halfLength != 0.0) {
                    double densityOne = 0.0;
                    double densityTwo = 0.0;
                    for (int j = 0; j < demes; j++) {
                        densityOne += (((lineageStartCount[j]*lineageStartCount[j])-lineageStartCountSquare[j])/(2.0*popSizes.getParameterValue(j)));
                        densityTwo += (((lineageEndCount[j]*lineageEndCount[j])-lineageEndCountSquare[j])/(2.0*popSizes.getParameterValue(j)));
                    }
                    logL += -halfLength*densityOne;
                    logL += -halfLength*densityTwo;
                }

                if (intervals.getIntervalType(i) == IntervalType.COALESCENT) {
                    Double[] leftProbs = new Double[demes];
                    Double[] rightProbs = new Double[demes];
                    for (int j = 0 ; j < demes; j++) {
                        leftProbs[j] = coalescentLeftProbs[i][j];
                        rightProbs[j] = coalescentRightProbs[i][j];
                    }

                } else {
                    //do nothing
                }

            } else {
                //do nothing
            }
        }

        return logL;

    }

    /**
     * Compute all required probability distributions for calculating the overall structured coalescent density.
     * This methods computes equations 11 and 12 from the BASTA manuscript.
     */
    private void computeProbabilityDistributions() {

        int intervalCount = intervals.getIntervalCount();
        System.out.println("number of intervals: " + intervalCount);

        //first sampling event is not considered so take this into account
        NodeRef node = intervals.getSamplingNode(0);
        System.out.println("sampling node: "+ treeModel.getNodeTaxon(node).getId());
        Double[] lineageCount = new Double[demes];
        for (int i = 0; i < demes; i++) {
            lineageCount[i] = 0.0;
        }
        //set start probabilities of first sampling interval
        lineageCount[patternList.getPattern(0)[patternList.getTaxonIndex(treeModel.getNodeTaxon(node).getId())]] = 1.0;
        this.intervalStartProbs[0] = Arrays.copyOf(lineageCount, demes);
        this.intervalStartSquareProbs[0] = Arrays.copyOf(lineageCount, demes);
        //add to list of active lineages
        this.activeLineages.add(lineageCount);

        //each interval is defined by its ending event, i.e. does the interval end in a sampling or coalescent event
        for (int i = 0; i < intervalCount; i++) {

            System.out.println("interval type: " + intervals.getIntervalType(i));

            if (intervals.getIntervalType(i) == IntervalType.COALESCENT) {
                //get coalescent node from interval
                node = intervals.getCoalescentNode(i);

                /*for (int j = 0; j < treeModel.getChildCount(node); j++) {
                    System.out.println("coalescent child node: " + treeModel.getChild(node, j));
                }*/

                //TODO add information to coalescentLeftProbs and coalescentRightProbs



            } else if (intervals.getIntervalType(i) == IntervalType.SAMPLE) {
                //check for zero-length interval
                if (intervals.getInterval(i) == 0.0) {
                    //multiple samples at same sampling time
                    System.out.println("zero-length interval");
                    node = intervals.getSamplingNode(i+1);
                    System.out.println("sampling node: "+ treeModel.getNodeTaxon(node).getId());


                } else {
                    //one sample at a time

                    //compute end probabilities of first sampling interval
                    double intervalLength = intervals.getIntervalTime(0);
                    //matrix exponentiation to compute end interval probabilities; equation 11
                    incrementActiveLineages(this.activeLineages, intervalLength, 0);
                    Double[] lineageCountSquare = new Double[demes];
                    for (int j = 0; j < demes; j++) {
                        lineageCountSquare[j] = lineageCount[j]*lineageCount[j];
                    }
                    //computing Sd (see supplementary)
                    Double[] sumActiveLineages = new Double[demes];
                    for (int j = 0; j < demes; j++) {
                        sumActiveLineages[j] = 0.0;
                    }
                    for (Double[] lineage : activeLineages) {
                        for (int j = 0; j < demes; j++) {
                            sumActiveLineages[j] += lineage[j];
                        }
                    }
                    this.intervalEndProbs[i] = Arrays.copyOf(sumActiveLineages, demes);
                    this.intervalEndSquareProbs[i] = Arrays.copyOf(lineageCountSquare, demes);

                    //get the node number of the sampling node
                    node = intervals.getSamplingNode(i+1);
                    System.out.println("sampling node: "+ treeModel.getNodeTaxon(node).getId());

                    //TODO add to the active lineages for next interval iteration



                }
            }
        }

    }

    private static void transpose(double[] matrix, int dim) {
        for (int i = 0; i  < dim; ++i) {
            for (int j = i + 1; j < dim; ++j) {
                final int ij = i * dim + j;
                final int ji = j * dim + i;
                double tmp = matrix[ij];
                matrix[ij] = matrix[ji];
                matrix[ji] = tmp;
            }
        }
    }

    private void incrementActiveLineages(ArrayList<Double[]> lineageCount, double intervalLength, int interval) {
        final double branchRate;
        synchronized (branchRateModel) {
            branchRate = branchRateModel.getBranchRate(treeModel, treeModel.getRoot());
        }

        if (!matricesKnown) {
            generalSubstitutionModel.getTransitionProbabilities(branchRate * intervalLength, migrationMatrices[interval]);
            if (USE_TRANSPOSE) {
                transpose(migrationMatrices[interval], demes);
            }

            if (MATRIX_DEBUG) {
                System.out.println("-----------");
                double[] matrix = new double[demes*demes];
                generalSubstitutionModel.getInfinitesimalMatrix(matrix);
                for (int i = 0; i < demes; i++) {
                    for (int j = 0; j < demes; j++) {
                        System.out.print(matrix[i*demes+j] + " ");
                    }
                    System.out.println();
                }
                System.out.println("-----------");
                System.out.println("Matrix exponentiation (t=" + interval + ") is: ");
                for (int i = 0; i < demes * demes; i++) {
                    System.out.print(migrationMatrices[interval][i] + " ");
                    if ((i + 1) % demes == 0) {
                        System.out.println();
                    }
                }
                System.out.println("-----------");
            }
        }

    }

    /**
     * When a lineage/branch has not been fully processed/computed towards the log likelihood, increase
     * all the branch lengths of those lineages still active.
     *
     * @param increment
     */
    private void incrementActiveLineages(double increment) {
        if (DEBUG) {
            System.out.println("Incrementing active lineages by " + increment);
        }

        final double branchRate;
        synchronized (branchRateModel) {
            branchRate = branchRateModel.getBranchRate(treeModel, treeModel.getRoot());
        }
        if (!matricesKnown) {
            generalSubstitutionModel.getTransitionProbabilities(branchRate * increment, migrationMatrices[this.currentCoalescentInterval]);
            if (USE_TRANSPOSE) {
                transpose(migrationMatrices[this.currentCoalescentInterval], demes);
            }

            if (MATRIX_DEBUG) {
                System.out.println("-----------");
                double[] matrix = new double[demes*demes];
                generalSubstitutionModel.getInfinitesimalMatrix(matrix);
                for (int i = 0; i < demes; i++) {
                    for (int j = 0; j < demes; j++) {
                        System.out.print(matrix[i*demes+j] + " ");
                    }
                    System.out.println();
                }
                System.out.println("-----------");
                System.out.println("Matrix exponentiation (t=" + increment + ") is: ");
                for (int i = 0; i < demes * demes; i++) {
                    System.out.print(migrationMatrices[this.currentCoalescentInterval][i] + " ");
                    if ((i + 1) % demes == 0) {
                        System.out.println();
                    }
                }
                System.out.println("-----------");
            }
        }

        //TODO evaluate later if this should be done in parallel
        //double start = System.currentTimeMillis();
        /*for (ProbDist pd : activeLineageList) {
            pd.incrementIntervalLength(increment, migrationMatrices[this.currentCoalescentInterval]);
            if (DEBUG) {
                System.out.println("  " + pd);
            }
        }*/
        for (int i = 0; i < addedLength; i++) {
            activeLineageList[i].incrementIntervalLength(increment, migrationMatrices[this.currentCoalescentInterval]);
        }
        //double end = System.currentTimeMillis();
        this.currentCoalescentInterval++;
    }

    // **************************************************************
    // ModelListener IMPLEMENTATION
    // **************************************************************

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (DEBUG) {
            System.out.println("handleModelChangedEvent: " + model.getModelName() + ", " + object + " (class " + object.getClass() + ")");
        }
        if (model == treeModel) {
            if (object instanceof TreeChangedEvent) {
                // If a node event occurs the node and its two child nodes
                // are flagged for updating (this will result in everything
                // above being updated as well. Node events occur when a node
                // is added to a branch, removed from a branch or its height or
                // rate changes.
                if (((TreeChangedEvent) object).isNodeChanged()) {
                    //System.out.println("isNodeChanged: " + ((TreeChangedEvent) object).getNode().getNumber());
                    //System.out.println(treeModel.getNodeHeight(((TreeChangedEvent) object).getNode());
                    //double changeHeight = treeModel.getNodeHeight(((TreeChangedEvent) object).getNode());
                    //TODO use what's in the current times variable to decide which ProbDist to update?
                    //TODO NOT SUFFICIENT: sketch out an example

                    if (MINIMUM_EVALUATION) {
                        NodeRef node = ((TreeChangedEvent) object).getNode();
                        int nodeNumber = ((TreeChangedEvent) object).getNode().getNumber();
                        this.nodeProbDist[nodeNumber].needsUpdate = true;
                        while (treeModel.getParent(node) != null && (treeModel.getParent(node) != treeModel.getRoot())) {
                            node = treeModel.getParent(node);
                            this.nodeProbDist[node.getNumber()].needsUpdate = true;
                        }
                    /*for (ProbDist pd : this.nodeProbDist) {
                        if (pd.needsUpdate) {
                            System.out.println("update: " + pd.node.getNumber() + " (total: " + this.nodeProbDist.length + ")");
                        }
                    }*/
                    } else {
                        //updateAllDensities(true);
                    }

                    //print out the tree to check
                    System.out.println(treeModel);
                    //check both times and storedTimes
                    /*System.out.println("times");
                    for (int i = 0; i < times.size(); i++) {
                        System.out.print(times.get(i) + " ");
                    }
                    System.out.println();
                    System.out.println("storedTimes");
                    for (int i = 0; i < storedTimes.size(); i++) {
                        System.out.print(storedTimes.get(i) + " ");
                    }
                    System.out.println();*/

                    //TODO call collectTimes? probably not as we need the old ones
                    //times.clear();
                    //children.clear();
                    //nodes.clear();
                    //collectAllTimes(this.treeModel, this.treeModel.getRoot(), nodes, times, children);

                    //check both times and storedTimes
                    //System.out.println("times");
                    /*for (int i = 0; i < times.size(); i++) {
                        System.out.print(times.get(i) + " ");
                    }
                    System.out.println();
                    System.out.println("storedTimes");
                    for (int i = 0; i < storedTimes.size(); i++) {
                        System.out.print(storedTimes.get(i) + " ");
                    }
                    System.out.println();*/

                    //TODO uncomment this
                    /*double minHeight = Math.min(changeHeight, );

                    for (ProbDist pd : this.nodeProbDist) {
                        if (treeModel.getNodeHeight(pd.node) >= minHeight) {
                            pd.needsUpdate = true;
                        }
                    }*/

                    //TODO give updateTransitionProbabilities more responsibility?
                } else if (((TreeChangedEvent) object).isHeightChanged()) {
                    //System.out.println("isHeightChanged: " + ((TreeChangedEvent) object).getNode().getNumber());

                } else if (((TreeChangedEvent) object).isTreeChanged()) {
                    //System.out.println("isTreeChanged");
                    // Full tree events result in a complete updating of the tree likelihood
                    // This event type is now used for EmpiricalTreeDistributions.
                    System.err.println("Full tree update event - these events currently aren't used\n" +
                            "so either this is in error or a new feature is using them so remove this message.");
                    //updateAllDensities(true);
                } else {
                    //System.out.println("else: ?");
                    // Other event types are ignored (probably trait changes).
                    System.err.println("Another tree event has occurred (possibly a trait change).");
                }
            }
            //for all the nodes that are older than the event, set needsUpdate to true
            //then trigger a recalculation that makes use of an adjusted traverseTree method (that checks whether
            //or not the ProbDist needs to be updated
            likelihoodKnown = false;
            //TODO not all matrices will have to be recomputed all the time
            matricesKnown = false;
            treeModelUpdateFired = true;
            areStatesRedrawn = false;
        } else if (model == branchRateModel) {
            matricesKnown = false;
            //the following to accommodate events stemming from the upDownOperator
            this.rateChanged = true;
            /*likelihoodKnown = false;
            for (ProbDist pd : this.nodeProbDist) {
                pd.needsUpdate = true;
            }*/
            //updateAllDensities();
            likelihoodKnown = false;
            areStatesRedrawn = false;
        } else if (model == generalSubstitutionModel) {
            matricesKnown = false;
            //TODO is this necessary? turns out it is to avoid store/restore issues but why??
            this.rateChanged = true;
            /*likelihoodKnown = false;
            for (ProbDist pd : this.nodeProbDist) {
                pd.needsUpdate = true;
            }*/
            //updateAllDensities();
            likelihoodKnown = false;
            areStatesRedrawn = false;
        } else {
            throw new RuntimeException("Unknown handleModelChangedEvent source, exiting.");
        }
    }

    // **************************************************************
    // VariableListener IMPLEMENTATION
    // **************************************************************


    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        if (DEBUG) {
            System.out.println("handleVariableChangedEvent: " + variable.getVariableName());
        }
        //if even one of the popSizes has changed, the whole density needs to be recomputed
        //hence we should use the adaptive multivariate transition kernel on all popSizes
        likelihoodKnown = false;
        areStatesRedrawn = false;
        //a change in one of the popSizes does not affect matrix exponentiation
        matricesKnown = true;
    }

    protected void storeState() {
        //super.storeState();
        for (int i = 0; i < this.finalCoalescentInterval; i++) {
            System.arraycopy(this.migrationMatrices[i],0, this.storedMigrationMatrices[i], 0, demes*demes);
        }
        storedLikelihoodKnown = likelihoodKnown;
        storedLogLikelihood = logLikelihood;

        if (areStatesRedrawn) {
            for (int i = 0; i < reconstructedStates.length; i++) {
                System.arraycopy(reconstructedStates[i], 0, storedReconstructedStates[i], 0, reconstructedStates[i].length);
            }
        }

        storedAreStatesRedrawn = areStatesRedrawn;
        //TODO look into the line below
        //storedJointLogLikelihood = jointLogLikelihood;
    }

    protected void restoreState() {
        //super.restoreState();
        double[] tmp;
        for (int i = 0; i < this.finalCoalescentInterval; i++) {
            tmp = this.migrationMatrices[i];
            this.migrationMatrices[i] = this.storedMigrationMatrices[i];
            this.storedMigrationMatrices[i] = tmp;
        }
        likelihoodKnown = storedLikelihoodKnown;
        logLikelihood = storedLogLikelihood;

        int[][] temp = reconstructedStates;
        reconstructedStates = storedReconstructedStates;
        storedReconstructedStates = temp;

        areStatesRedrawn = storedAreStatesRedrawn;
        //TODO look into the line below
        //jointLogLikelihood = storedJointLogLikelihood;
    }

    @Override
    protected void acceptState() {
        // do nothing
    }

    public void makeDirty() {
        //updateAllDensities(true);
        likelihoodKnown = false;
        matricesKnown = false;
        areStatesRedrawn = false;
    }

    public TreeModel getTreeModel() {
        return this.treeModel;
    }

    /**
     * Private class that allows for objects that hold the computed probability distribution of lineages among demes
     */
    private class ProbDist {

        //lineage probability distribution at start of interval
        private double[] startLineageProbs;
        //lineage probability distribution at end of interval
        private double[] endLineageProbs;
        //private double[] expectedLineages;
        private double intervalLength;

        //this node variable serves a double purpose as it's either a tip node in the case of a sampling event
        //or an internal node in the case of a coalescent event
        private NodeRef node;

        private int patternIndex;

        //keep track of whether the interval length has already been incremented for matrix exponentiation
        private boolean incremented = false;

        //internal boolean that indicates if this density needs to be recomputed
        private boolean needsUpdate = true;

        //coalescent or sample
        private IntervalType intervalType;

        //child nodes only serve a purpose when dealing with a coalescent event
        private NodeRef leftChild = null;
        private NodeRef rightChild = null;

        public ProbDist(int nDemes) {
            this.startLineageProbs = new double[nDemes];
            this.endLineageProbs = new double[nDemes];
        }

        public ProbDist(int nDemes, double distance, NodeRef node) {
            this.startLineageProbs = new double[nDemes];
            this.endLineageProbs = new double[nDemes];
            this.intervalLength = distance;
            this.node = node;
        }

        public ProbDist(int nDemes, double distance, NodeRef node, NodeRef leftChild, NodeRef rightChild) {
            this.startLineageProbs = new double[nDemes];
            this.endLineageProbs = new double[nDemes];
            this.intervalLength = distance;
            this.node = node;
            this.leftChild = leftChild;
            this.rightChild = rightChild;
        }

        public void update(double distance) {
            this.intervalLength = distance;
            for (int i = 0; i < startLineageProbs.length; i++) {
                this.startLineageProbs[i] = 0.0;
                this.endLineageProbs[i] = 0.0;
            }
            this.startLineageProbs[this.patternIndex] = 1.0;
            this.endLineageProbs[this.patternIndex] = 1.0;
        }

        public void update(double distance, NodeRef node, IntervalType intervalType, NodeRef leftChild, NodeRef rightChild) {
            this.intervalLength = distance;
            this.node = node;
            this.intervalType = intervalType;
            this.leftChild = leftChild;
            this.rightChild = rightChild;
            this.incremented = false;
            for (int i = 0; i < startLineageProbs.length; i++) {
                this.startLineageProbs[i] = 0.0;
                this.endLineageProbs[i] = 0.0;
            }
            this.startLineageProbs[this.patternIndex] = 1.0;
            this.endLineageProbs[this.patternIndex] = 1.0;
        }

        //compute the probability distribution of lineages among demes for a coalescent event
        public double computeCoalescedLineage(ProbDist leftProbDist, ProbDist rightProbDist) {
            double sum = 0.0;
            double[] sumComponents = new double[demes];
            for (int i = 0; i < demes; i++) {
                sumComponents[i] = (leftProbDist.endLineageProbs[i] * rightProbDist.endLineageProbs[i])/popSizes.getParameterValue(i);
                sum += sumComponents[i];
            }

            if (DEBUG) {
                System.out.println("coalescent end lineage prob 0 = " + leftProbDist.endLineageProbs[0]);
                System.out.println("coalescent end lineage prob 1 = " + leftProbDist.endLineageProbs[1]);
            }
            for (int i = 0; i < demes; i++) {
                this.startLineageProbs[i] = sumComponents[i]/sum;
                if (DEBUG) {
                    System.out.println("coalescent start lineage prob " + i + " = " + (sumComponents[i]/sum));
                }
            }
            if (DEBUG) {
                System.out.println("E_i = " + Math.log(sum));
            }
            //new interval (and new node) so set its length to zero
            this.intervalLength = 0.0;
            return Math.log(sum);
        }

        //compute the end probability densities and expected numbers of lineages
        public void computeEndLineageDensities(double lineageLength, double[] migrationMatrix) {
            //TODO this should be possible in parallel for each lineage within the same coalescent interval
            //TODO will need to collect all the tasks first (i.e. collectEndLineageDensities) and then start in parallel
            for (int k = 0; k < demes; k++) {
                    /*double value = 0.0;
                    for (int l = 0; l < demes; l++) {
                        value += this.startLineageProbs[l] * migrationMatrix[l*demes+k];
                    }
                    this.endLineageProbs[k] = value;*/
                this.endLineageProbs[k] = USE_TRANSPOSE ?
                        rdot(demes, startLineageProbs, 0, 1, migrationMatrix, k * demes, 1) :
                        rdot(demes, startLineageProbs, 0, 1, migrationMatrix, k, demes);
            }
            //}

        }

        //TODO this should really become 1 call to computeEndLineageDensities
        public void incrementIntervalLength(double increment, double[] migrationMatrix) {
            this.intervalLength += increment;
            if (incremented) {
                if (ASSOC_MULTIPLICATION) {
                    System.arraycopy(this.endLineageProbs, 0, this.startLineageProbs, 0, demes);
                    //for (int i = 0; i < demes; i++) {
                    //    startLineageProbs[i] = endLineageProbs[i];
                    //}
                    computeEndLineageDensities(increment, migrationMatrix);
                } else {
                    throw new RuntimeException("Only incremental matrix exponentiation allowed for performance reasons.");
                }
            } else {
                if (DEBUG) {
                    System.out.println("this.intervalLength computeEndLineageDensities: " + this.intervalLength + " ; increment = " + increment);
                }
                computeEndLineageDensities(this.intervalLength, migrationMatrix);
            }
            this.incremented = true;
        }

        public IntervalType getIntervalType() {
            return this.intervalType;
        }

        public NodeRef getLeftChild() {
            return this.leftChild;
        }

        public NodeRef getRightChild() {
            return this.rightChild;
        }

        @Override
        public String toString() {
            String output = "Node " + this.node + " ; length = " + this.intervalLength + " S(";
            for (int i = 0; i < startLineageProbs.length; i++) {
                output += startLineageProbs[i] + " ";
            }
            output += ") E(";
            for (int i = 0; i < endLineageProbs.length; i++) {
                output += endLineageProbs[i] + " ";
            }
            output += ")";
            return output;
        }

    }

    // **************************************************************
    // Method from JavaBlas v1.2.4 (www.jblas.org)
    // **************************************************************

    /** Compute scalar product between dx and dy. */
    public static double rdot(int n, double[] dx, int dxIdx, int incx, double[] dy, int dyIdx, int incy) {
        double s = 0.0;
        if (incx == 1 && incy == 1 && dxIdx == 0 && dyIdx == 0) {
            for (int i = 0; i < n; i++)
                s += dx[i] * dy[i];
        }
        else {
            for (int c = 0, xi = dxIdx, yi = dyIdx; c < n; c++, xi += incx, yi += incy) {
                s += dx[xi] * dy[yi];
            }
        }
        return s;
    }

    // **************************************************************
    // Units IMPLEMENTATION
    // **************************************************************

    /**
     * Sets the units these coalescent intervals are measured in.
     */
    public final void setUnits(Type u) {
        treeModel.setUnits(u);
    }

    /**
     * Returns the units these coalescent intervals are measured in.
     */
    public final Type getUnits() {
        return treeModel.getUnits();
    }

    // **************************************************************
    // Citable IMPLEMENTATION
    // **************************************************************

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.TREE_PRIORS;
    }

    @Override
    public String getDescription() {
        return "Bayesian structured coalescent approximation";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CITATION);
    }

    public static Citation CITATION = new Citation(
            new Author[]{
                    new Author("Nicola", "De Maio"),
                    new Author("Chieh-Hsi", "Wu"),
                    new Author("Kathleen", "O'Reilly"),
                    new Author("Daniel", "Wilson")
            },
            "New routes to phylogeography: a Bayesian structured coalescent approximation",
            2015,
            "PLOS Genetics",
            11, "e1005421",
            "10.1371/journal.pgen.1005421"
    );

    // ****************************************************************
    // Private and protected stuff
    // ****************************************************************

    protected double logLikelihood;
    protected double storedLogLikelihood;
    protected boolean likelihoodKnown = false;
    protected boolean storedLikelihoodKnown = false;

    //the tree model
    private TreeModel treeModel;

    //the branch rate model
    private BranchRateModel branchRateModel;

    //the population sizes
    private Parameter popSizes;

    //the discrete trait data
    private PatternList patternList;

    //objects for reconstructing ancestral states
    private final DataType dataType;
    private boolean useMAP;
    private int[][] reconstructedStates;
    private int[][] storedReconstructedStates;
    protected boolean areStatesRedrawn = false;
    protected boolean storedAreStatesRedrawn = false;

    //expected starting lineage counts
    private double[] startExpected;

    //expected ending lineage counts
    private double[] endExpected;

    //array with a probability distribution for each node
    private ProbDist[] nodeProbDist;

    //elements for constructing the coalescent intervals and times
    private ArrayList<ComparableDouble> times;
    private ArrayList<Integer> children;
    private ArrayList<NodeRef> nodes;
    private int[] indices;

    //stored versions of the coalescent intervals and times
    private ArrayList<ComparableDouble> storedTimes;
    private ArrayList<Integer> storedChildren;
    private ArrayList<NodeRef> storedNodes;
    private int[] storedIndices;

    //list of currently active lineages
    //private ArrayList<ProbDist> activeLineageList;
    //temporary list of lineages
    //private ArrayList<ProbDist> tempLineageList;

    //list of currently active lineages
    private ProbDist[] activeLineageList;
    //temporary list of lineages
    //private ProbDist[] tempLineageList;
    private boolean[] addedLineages;
    private int addedLength;

    private BigFastTreeIntervals intervals;

    //probability densities at the start and end of each coalescent interval
    private Double[][] intervalStartProbs;
    private Double[][] intervalEndProbs;
    private Double[][] intervalStartSquareProbs;
    private Double[][] intervalEndSquareProbs;
    //keep track of information to compute structured coalescent density at coalescent events
    private Double[][] coalescentLeftProbs;
    private Double[][] coalescentRightProbs;
    //keep track of (varying) list of active lineages
    private ArrayList<Double[]> activeLineages;
    //private Double[][] currentLineageProbs;

    //the migration model
    private GeneralSubstitutionModel generalSubstitutionModel;

    //number of demes for the structured coalescent model
    private int demes;

    //variables that allow to use storeState and restoreState
    private int maxCoalescentIntervals;
    private int currentCoalescentInterval;
    private double[][] migrationMatrices;
    private int finalCoalescentInterval;
    private double[][] storedMigrationMatrices;

    private boolean matricesKnown;
    private boolean rateChanged;
    private boolean treeModelUpdateFired;

}

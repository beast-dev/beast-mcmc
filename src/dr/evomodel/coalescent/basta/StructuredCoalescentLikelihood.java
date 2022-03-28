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
import dr.evolution.coalescent.IntervalType;
import dr.evolution.datatype.Codons;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneralDataType;
import dr.evolution.tree.*;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
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
public class StructuredCoalescentLikelihood extends AbstractModelLikelihood implements Units, Citable, AncestralStateTraitProvider, TreeTraitProvider {

    //TODO: the likelihood class should have minimum functionality, i.e. if the likelihood needs to be recomputed
    //then compute it by calling a StructuredCoalescentModel calculateLikelihood method; if not, then simply return it
    //TODO: create StructuredCoalescentModel class
    //TODO: have StructuredCoalescentLikelihood listen to both BeagleMatrixExponentiationDelegate and
    //the StructuredCoalescentModel classes?

    private static final boolean DEBUG = false;
    private static final boolean MATRIX_DEBUG = false;
    private static final boolean UPDATE_DEBUG = false;

    //private static final boolean USE_BEAGLE = false;
    private static final boolean ASSOC_MULTIPLICATION = true;
    private static final boolean USE_TRANSPOSE = true;

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
            System.out.println((TreeModel) tree);
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

        reconstructedStates = new int[treeModel.getNodeCount()][patternList.getPatternCount()];
        storedReconstructedStates = new int[treeModel.getNodeCount()][patternList.getPatternCount()];

        treeTraits.addTrait(new TreeTrait.IA() {
            public String getTraitName() {
                return tag;
            }

            public Intent getIntent() {
                return Intent.NODE;
            }

            public Class getTraitClass() {
                return int[].class;
            }

            public int[] getTrait(Tree tree, NodeRef node) {
                return getStatesForNode(tree, node);
            }

            public String getTraitString(Tree tree, NodeRef node) {
                return formattedState(getStatesForNode(tree, node), dataType);
            }
        });

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
            updateAllDensities(false);
        }
        //System.out.println("total time by incrementing active lineages = " + this.timeIncrementActiveLineages/1000.0);
        return logLikelihood;
    }

    /**
     * Calculates the log likelihood of this set of coalescent intervals,
     * given a demographic model.
     */
    public double calculateLogLikelihood() {
        areStatesRedrawn = false;

        if (DEBUG) {
            System.out.println("\nStructuredCoalescentLikelihood.calculateLogLikelihood():");
            System.out.println("Coalescent Intervals list of length: " + (this.treeModel.getNodeCount() - 1));
            System.out.println("Tree: " + treeModel.getNewick() + "\n");

            System.out.println("Pattern information: " + patternList.getPatternCount() + " pattern(s) of type: " + patternList.getDataType());
            for (int i = 0; i < patternList.getTaxonCount(); i++) {
                System.out.println("  Taxon " + patternList.getTaxon(i) + ": " + patternList.getPattern(0)[i]);
            }
            System.out.println();
        }

        logLikelihood = traverseTree(treeModel, treeModel.getRoot(), patternList);
        redrawAncestralStates();

        return logLikelihood;
    }

    //based on the traverseTree method in OldAbstractCoalescentLikelihood
    private double traverseTree(Tree tree, NodeRef root, PatternList patternList) {
        double MULTIFURCATION_LIMIT = 1e-9;

        if (treeModelUpdateFired && !rateChanged) {
            if (DEBUG) {
                System.out.println("treeModelUpdateFired");
                System.out.println("current coalescent intervals:");
                for (int i = 1; i < times.size(); i++) {
                    System.out.println((i-1) + ": " + (times.get(indices[i]).doubleValue() - times.get(indices[i-1]).doubleValue()));
                }
            }
            updateTransitionProbabilities();
            treeModelUpdateFired = false;
        } else {
            times.clear();
            children.clear();
            nodes.clear();
            collectAllTimes(tree, root, nodes, times, children);
            indices = new int[times.size()];

            HeapSort.sort(times, indices);
        }

        double lnL = 0.0;

        // start is the time of the first tip
        double start = times.get(indices[0]).doubleValue();

        int i = 0;
        int j = 0;
        while (i < times.size()) {

            int lineagesRemoved = 0;
            int lineagesAdded = 0;

            final double finish = times.get(indices[i]).doubleValue();
            double next = finish;

            if (DEBUG) {
                System.out.println("start = " + start + " ; finish = " + finish);
            }

            double intervalLength = finish - start;

            while (Math.abs(next - finish) < MULTIFURCATION_LIMIT) {
                final int children = this.children.get(indices[i]);
                if (children == 0) {
                    lineagesAdded += 1;
                } else {
                    lineagesRemoved += (children - 1);
                }
                i += 1;
                if (i == times.size()) break;

                next = times.get(indices[i]).doubleValue();
            }

            if (lineagesAdded > 0) {
                if (DEBUG){
                    System.out.println("lineagesAdded = " + lineagesAdded);
                }
                if ((intervalLength) > MULTIFURCATION_LIMIT) {
                    incrementActiveLineages(finish-start);
                    if (DEBUG) {
                        System.out.println("start = " + start + " ; finish = " + finish + " ; intervalCount > 0");
                    }
                    while (Math.abs(treeModel.getNodeHeight(nodes.get(indices[j])) - finish) < MULTIFURCATION_LIMIT) {
                        NodeRef refNode = nodes.get(indices[j]);
                        if (treeModel.isExternal(refNode)) {
                            //ProbDist newProbDist = new ProbDist(demes, 0.0, refNode);
                            //nodeProbDist[refNode.getNumber()].update( 0.0, refNode, IntervalType.SAMPLE);
                            if (nodeProbDist[refNode.getNumber()].needsUpdate) {
                                nodeProbDist[refNode.getNumber()].update(0.0);
                            }
                            //newProbDist.setIntervalType(IntervalType.SAMPLE);
                            //newProbDist.startLineageProbs[patternList.getPattern(0)[patternList.getTaxonIndex(treeModel.getNodeTaxon(newProbDist.node).getId())]] =  1.0;
                            //newProbDist.copyLineageDensities();
                            //nodeProbDist[refNode.getNumber()].startLineageProbs[patternList.getPattern(0)[patternList.getTaxonIndex(treeModel.getNodeTaxon(refNode).getId())]] =  1.0;
                            //nodeProbDist[refNode.getNumber()].copyLineageDensities();
                            if (DEBUG) {
                                System.out.println(nodeProbDist[refNode.getNumber()]);
                            }
                            //temporary list required to keep accurate track of expected lineage counts
                            //tempLineageList.add(nodeProbDist[refNode.getNumber()]);
                            addedLineages[refNode.getNumber()] = true;
                        } else {
                            if (treeModel.getChildCount(refNode) > 2) {
                                throw new RuntimeException("Structured coalescent currently only allows strictly bifurcating trees.");
                            }
                            NodeRef leftChild = treeModel.getChild(refNode, 0);
                            NodeRef rightChild = treeModel.getChild(refNode, 1);
                            //ProbDist newProbDist = new ProbDist(demes, intervalLength, refNode, leftChild, rightChild);
                            //newProbDist.setIntervalType(IntervalType.COALESCENT);
                            if (nodeProbDist[refNode.getNumber()].needsUpdate) {
                                nodeProbDist[refNode.getNumber()].update(intervalLength, refNode, IntervalType.COALESCENT, leftChild, rightChild);
                            }
                            //temporary list required to keep accurate track of expected lineage counts
                            //tempLineageList.add(nodeProbDist[refNode.getNumber()]);
                            addedLineages[refNode.getNumber()] = true;
                        }
                        j++;
                        if (j >= indices.length) {
                            j = 0;
                            break;
                        }
                    }
                } else {
                    //very first sampling event
                    if (DEBUG) {
                        System.out.println("first sampling event");
                    }
                    while (Math.abs(treeModel.getNodeHeight(nodes.get(indices[j])) - start) < MULTIFURCATION_LIMIT) {
                        NodeRef refNode = nodes.get(indices[j]);
                        if (treeModel.isExternal(refNode)) {
                            //ProbDist newProbDist = new ProbDist(demes, 0.0, refNode);
                            //nodeProbDist[refNode.getNumber()].update( 0.0, refNode, IntervalType.SAMPLE);
                            if (nodeProbDist[refNode.getNumber()].needsUpdate) {
                                nodeProbDist[refNode.getNumber()].update(0.0);
                            }
                            //newProbDist.setIntervalType(IntervalType.SAMPLE);
                            //newProbDist.startLineageProbs[patternList.getPattern(0)[patternList.getTaxonIndex(treeModel.getNodeTaxon(newProbDist.node).getId())]] =  1.0;
                            //newProbDist.copyLineageDensities();
                            //nodeProbDist[refNode.getNumber()].startLineageProbs[patternList.getPattern(0)[patternList.getTaxonIndex(treeModel.getNodeTaxon(refNode).getId())]] =  1.0;
                            //nodeProbDist[refNode.getNumber()].copyLineageDensities();
                            if (DEBUG) {
                                System.out.println(nodeProbDist[refNode.getNumber()]);
                            }
                            //temporary list required to keep accurate track of expected lineage counts
                            //tempLineageList.add(nodeProbDist[refNode.getNumber()]);
                            addedLineages[refNode.getNumber()] = true;
                        } else {
                            throw new RuntimeException("First interval cannot be a coalescent event.");
                        }
                        //System.out.println("** (j=" + j + ") " + treeModel.getNodeHeight(nodes.get(indices[j])));
                        j++;
                        //System.out.println("** (j=" + j + ") ");
                        if (j >= indices.length) {
                            j = 0;
                            break;
                        }
                    }
                }
                start = finish;
            }

            if (lineagesRemoved > 0) {
                if (DEBUG) {
                    System.out.println("lineagesRemoved = " + lineagesRemoved);
                    System.out.println("start = " + start + " ; finish = " + finish);
                }
                incrementActiveLineages(finish-start);
                while (Math.abs(treeModel.getNodeHeight(nodes.get(indices[j])) - finish) < MULTIFURCATION_LIMIT) {
                    NodeRef refNode = nodes.get(indices[j]);
                    if (treeModel.isExternal(refNode)) {
                        //ProbDist newProbDist = new ProbDist(demes, intervalLength, refNode);
                        //nodeProbDist[refNode.getNumber()].update( 0.0, refNode, IntervalType.SAMPLE);
                        if (nodeProbDist[refNode.getNumber()].needsUpdate) {
                            nodeProbDist[refNode.getNumber()].update(0.0);
                        }
                        //newProbDist.setIntervalType(IntervalType.SAMPLE);
                        //newProbDist.startLineageProbs[patternList.getPattern(0)[patternList.getTaxonIndex(treeModel.getNodeTaxon(newProbDist.node).getId())]] = 1.0;
                        //newProbDist.copyLineageDensities();
                        //nodeProbDist[refNode.getNumber()].startLineageProbs[patternList.getPattern(0)[patternList.getTaxonIndex(treeModel.getNodeTaxon(refNode).getId())]] = 1.0;
                        //nodeProbDist[refNode.getNumber()].copyLineageDensities();
                        //temporary list required to keep accurate track of expected lineage counts
                        //tempLineageList.add(nodeProbDist[refNode.getNumber()]);
                        addedLineages[refNode.getNumber()] = true;
                    } else {
                        if (treeModel.getChildCount(refNode) > 2) {
                            throw new RuntimeException("Structured coalescent currently only allows strictly bifurcating trees.");
                        }
                        NodeRef leftChild = treeModel.getChild(refNode, 0);
                        NodeRef rightChild = treeModel.getChild(refNode, 1);
                        ProbDist leftProbDist = nodeProbDist[leftChild.getNumber()];
                        ProbDist rightProbDist = nodeProbDist[rightChild.getNumber()];
                        //tempLineageList.remove(leftProbDist);
                        //tempLineageList.remove(rightProbDist);
                        addedLineages[leftChild.getNumber()] = false;
                        addedLineages[rightChild.getNumber()] = false;
                        /*if (DEBUG) {
                            System.out.println("currently active lineages: " + activeLineageList.size());
                        }*/
                        //ProbDist newProbDist = new ProbDist(demes, intervalLength, refNode, leftChild, rightChild);
                        //newProbDist.setIntervalType(IntervalType.COALESCENT);
                        if (nodeProbDist[refNode.getNumber()].needsUpdate) {
                            nodeProbDist[refNode.getNumber()].update(intervalLength, refNode, IntervalType.COALESCENT, leftChild, rightChild);
                        }
                        //lnL += newProbDist.computeCoalescedLineage(leftProbDist, rightProbDist);
                        lnL += nodeProbDist[refNode.getNumber()].computeCoalescedLineage(leftProbDist, rightProbDist);
                        if (!treeModel.isRoot(refNode)) {
                            //temporary list required to keep accurate track of expected lineage counts
                            //tempLineageList.add(nodeProbDist[refNode.getNumber()]);
                            addedLineages[refNode.getNumber()] = true;
                        } else {
                            if (DEBUG) {
                                System.out.println("ROOT");
                            }
                        }
                    }
                    //System.out.println("** (j=" + j + ") " + treeModel.getNodeHeight(nodes.get(indices[j])));
                    j++;
                    //System.out.println("** (j=" + j + ") ");
                    if (j >= indices.length) {
                        break;
                    }
                }
                start = finish;
            }

            //compute expected lineage counts here
            if (!(finish == 0.0)) {
                if (DEBUG) {
                    System.out.println("Computing expected lineage counts: " + intervalLength);
                }
                //computeExpectedLineageCounts();
                //and then compute the log likelihood
                lnL += computeLogLikelihood(intervalLength);
                if (DEBUG) {
                    System.out.println("Computing (log) likelihood contribution; logP = " + lnL);
                }
            }

            //and then compute the log likelihood
            /*if (!(finish == 0.0)) {
                lnL += computeLogLikelihood(intervalLength);
                if (DEBUG) {
                    System.out.println("Computing (log) likelihood contribution; logP = " + lnL);
                }
            }*/

            //update the list of active lineages for next iteration
            /*activeLineageList.clear();
            for (ProbDist pd : tempLineageList) {
                activeLineageList.add(pd);
            }*/
            this.addedLength = 0;
            for (int k = 0; k < addedLineages.length; k++) {
                if (addedLineages[k]) {
                    activeLineageList[addedLength] = nodeProbDist[k];
                    addedLength++;
                }
            }
            /*if (DEBUG) {
                System.out.println("currently active lineages: " + activeLineageList.size() + "\n");
            }*/

        }

        if (DEBUG) {
            System.out.println("Structured coalescent lnL = " + lnL);
        }

        this.finalCoalescentInterval = this.currentCoalescentInterval;
        //System.out.println("finalCoalescentInterval = " + finalCoalescentInterval);
        this.currentCoalescentInterval = 0;
        this.matricesKnown = true;
        this.rateChanged = false;

        return lnL;
    }

    /**
     * Compute the (log) likelihood over all currently active lineages.
     *
     * @return log likelihood contribution over the currently active lineages
     */
    private double computeLogLikelihood(double intervalLength) {

        double intervalOne = 0.0;
        double intervalTwo = 0.0;

        if (DEBUG) {
            System.out.println("interval length = " + intervalLength);
            //System.out.println("active lineage list length = " + activeLineageList.size());
        }

        /*for (ProbDist pd : activeLineageList) {
            for (int i = 0; i < demes; i++) {
                startProbs[i] += pd.startLineageProbs[i]*pd.startLineageProbs[i];
                endProbs[i] += pd.endLineageProbs[i]*pd.endLineageProbs[i];
            }
        }*/

        double[] startProbs = new double[demes];
        double[] endProbs = new double[demes];
        for (int i = 0; i < this.demes; i++) {
            this.startExpected[i] = 0.0;
            this.endExpected[i] = 0.0;
        }
        for (int j = 0; j < addedLength; j++) {
            ProbDist pd = activeLineageList[j];
            for (int i = 0; i < this.demes; i++) {
                startProbs[i] += pd.startLineageProbs[i]*pd.startLineageProbs[i];
                endProbs[i] += pd.endLineageProbs[i]*pd.endLineageProbs[i];
                this.startExpected[i] += pd.startLineageProbs[i];
                this.endExpected[i] += pd.endLineageProbs[i];
            }
        }

        //TODO these two products can be stored inside of the ProbDist so as to only recompute when absolutely necessary
        /*for (int i = 0; i < addedLength; i++) {
            ProbDist pd = activeLineageList[i];
            for (int j = 0; j < demes; j++) {
                startProbs[j] += pd.startLineageProbs[j]*pd.startLineageProbs[j];
                endProbs[j] += pd.endLineageProbs[j]*pd.endLineageProbs[j];
            }
        }*/

        for (int i = 0; i < demes; i++) {
            intervalOne += (startExpected[i]*startExpected[i] - startProbs[i])/popSizes.getParameterValue(i);
            intervalTwo += (endExpected[i]*endExpected[i] - endProbs[i])/popSizes.getParameterValue(i);
        }

        intervalOne *= -intervalLength/4.0;
        intervalTwo *= -intervalLength/4.0;

        if (DEBUG) {
            System.out.println("interval 1 (log) likelihood = " + intervalOne);
            System.out.println("interval 2 (log) likelihood = " + intervalTwo);
            System.out.println("total (log) likelihood = " + (intervalOne + intervalTwo));
        }

        /*for (ProbDist pd : this.nodeProbDist) {
            pd.needsUpdate = false;
        }*/

        return intervalOne + intervalTwo;
    }

    /**
     * Iterate over all the currently active lineages and compute the expected lineage counts.
     */
    private void computeExpectedLineageCounts() {
        double[] start, end;
        for (int i = 0; i < this.demes; i++) {
            this.startExpected[i] = 0.0;
            this.endExpected[i] = 0.0;
        }
        /*for (ProbDist pd : this.activeLineageList) {
            start = pd.startLineageProbs;
            end = pd.endLineageProbs;
            for (int i = 0; i < this.demes; i++) {
                this.startExpected[i] += start[i];
                this.endExpected[i] += end[i];
            }
        }*/
        for (int j = 0; j < addedLength; j++) {
            ProbDist pd = activeLineageList[j];
            start = pd.startLineageProbs;
            end = pd.endLineageProbs;
            for (int i = 0; i < this.demes; i++) {
                this.startExpected[i] += start[i];
                this.endExpected[i] += end[i];
            }
        }

        if (DEBUG) {
            System.out.print("  E_start(");
            for (int i = 0; i < this.demes; i++) {
                System.out.print(this.startExpected[i] + " ");
            }
            System.out.print(")\n  E_end(");
            for (int i = 0; i < this.demes; i++) {
                System.out.print(this.endExpected[i] + " ");
            }
            System.out.println(")");
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

    /**
     * Extract coalescent times and tip information into ArrayList times from tree.
     * Upon return times contain the time of each node in the subtree below top, and at the corresponding index
     * of children is the descendant count for that time.
     *
     * @param top          the node to start from
     * @param tree         given tree
     * @param times        array to fill with times
     * @param children     array to fill with descendents count
     */
    private void collectAllTimes(Tree tree, NodeRef top, ArrayList<NodeRef> nodes,
                                        ArrayList<ComparableDouble> times, ArrayList<Integer> children) {

        times.add(new ComparableDouble(tree.getNodeHeight(top)));
        nodes.add(top);
        children.add(tree.getChildCount(top));

        for (int i = 0; i < tree.getChildCount(top); i++) {
            NodeRef child = tree.getChild(top, i);
            collectAllTimes(tree, child, nodes, times, children);
        }
    }

    /**
     * compare the two times lists and only recompute those matrices that correspond to different coalescent
     * interval lengths
     */
    private void updateTransitionProbabilities() {

        if (MATRIX_DEBUG) {
            int index = 0;
            for (int i = 1; i < storedTimes.size(); i++) {
                double length = storedTimes.get(storedIndices[i]).doubleValue() - storedTimes.get(storedIndices[i-1]).doubleValue();
                if (length != 0.0) {
                    System.out.println("----------- START ----------- updateTransitionProbabilities()");
                    System.out.println("Matrix exponentiation (t=" + length + ") is: ");
                    for (int j = 0; j < demes * demes; j++) {
                        System.out.print(migrationMatrices[index][j] + " ");
                        if ((j + 1) % demes == 0) {
                            System.out.println();
                        }
                    }
                    System.out.println("-----------");
                    index++;
                }
            }
        }

        double[] oldLengths = new double[times.size()-1];
        for (int i = 1; i < times.size(); i++) {
            oldLengths[i-1] = times.get(indices[i]).doubleValue() - times.get(indices[i-1]).doubleValue();
        }

        if (DEBUG) {
            System.out.println("old coalescent intervals:");
            for (int i = 1; i < times.size(); i++) {
                System.out.println((i - 1) + ": " + (times.get(indices[i]).doubleValue() - times.get(indices[i - 1]).doubleValue()));
            }
        }

        storedTimes.clear();
        storedNodes.clear();
        storedChildren.clear();

        collectAllTimes(treeModel, treeModel.getRoot(), storedNodes, storedTimes, storedChildren);

        storedIndices = new int[storedTimes.size()];
        HeapSort.sort(storedTimes, storedIndices);

        if (DEBUG) {
            System.out.println("new coalescent intervals:");
            for (int i = 1; i < storedTimes.size(); i++) {
                System.out.println((i - 1) + ": " + (storedTimes.get(storedIndices[i]).doubleValue() - storedTimes.get(storedIndices[i - 1]).doubleValue()));
            }
        }

        double[] newLengths = new double[storedTimes.size()-1];
        for (int i = 1; i < storedTimes.size(); i++) {
            newLengths[i-1] = storedTimes.get(storedIndices[i]).doubleValue() - storedTimes.get(storedIndices[i-1]).doubleValue();
        }

        final double branchRate;
        synchronized (branchRateModel) {
            branchRate = branchRateModel.getBranchRate(treeModel, treeModel.getRoot());
        }

        if (UPDATE_DEBUG) {
            for (int i = 0; i < times.size(); i++) {
                System.out.print(times.get(i) + " ");
            }
            System.out.println();
            for (int i = 0; i < storedTimes.size(); i++) {
                System.out.print(storedTimes.get(i) + " ");
            }
            System.out.println();
        }

        int matrixIndex = 0;

        //take into account interval length of 0.0
        //if the number of intervals would increase, currently throw an error but full evaluation would be better
        if (oldLengths.length == newLengths.length) {
            for (int i = 0; i < oldLengths.length; i++) {
                if ((oldLengths[i] != newLengths[i]) && (newLengths[i] != 0.0)) {
                    if (UPDATE_DEBUG) {
                        System.out.println("updating transition probability matrix for coalescent interval length " + newLengths[i] + " at index " + matrixIndex);
                    }
                    generalSubstitutionModel.getTransitionProbabilities(branchRate * newLengths[i], migrationMatrices[matrixIndex]);
                    if (USE_TRANSPOSE) {
                        transpose(migrationMatrices[matrixIndex], demes);
                    }
                    matrixIndex++;
                } else {
                    if (MATRIX_DEBUG) {
                        System.out.println("no matrix update required");
                    }
                    if (newLengths[i] != 0.0) {
                        matrixIndex++;
                    }
                }
            }
        } else {
            throw new RuntimeException("Number of coalescent intervals has increased?");
        }

        if (MATRIX_DEBUG) {
            int index = 0;
            for (int i = 1; i < storedTimes.size(); i++) {
                double length = storedTimes.get(storedIndices[i]).doubleValue() - storedTimes.get(storedIndices[i-1]).doubleValue();
                if (length != 0.0) {
                    System.out.println("----------- END ----------- updateTransitionProbabilities()");
                    System.out.println("Matrix exponentiation (t=" + length + ") is: ");
                    for (int j = 0; j < demes * demes; j++) {
                        System.out.print(migrationMatrices[index][j] + " ");
                        if ((j + 1) % demes == 0) {
                            System.out.println();
                        }
                    }
                    System.out.println("-----------");
                    index++;
                }
            }
        }

        times = storedTimes;
        nodes = storedNodes;
        children = storedChildren;
        indices = storedIndices;

        matricesKnown = true;

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

                    //print out the tree to check
                    //System.out.println(treeModel);
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
                    updateAllDensities(true);
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
        updateAllDensities(true);
        likelihoodKnown = false;
        matricesKnown = false;
        areStatesRedrawn = false;
    }

    protected void updateAllDensities(boolean status) {
        for (ProbDist pd : this.nodeProbDist) {
            pd.needsUpdate = status;
        }
        //likelihoodKnown = false;
    }

    public void redrawAncestralStates() {
        //jointLogLikelihood = 0;
        traverseSample(treeModel, treeModel.getRoot());
        areStatesRedrawn = true;
    }

    //TODO DTA employs the parent state for this; need this here as well?
    public void traverseSample(TreeModel tree, NodeRef node) {
        //System.out.println("traverseSample: " + node.toString());

        //int[] state = new int[patternList.getPatternCount()];

        if (!tree.isExternal(node)) {
            //internal node including the root
            double[] liks = nodeProbDist[node.getNumber()].startLineageProbs;
            /*System.out.println("length: " + liks.length);
            for (int i = 0; i < liks.length; i++) {
                System.out.println(liks[i]);
            }*/
            for (int j = 0; j < patternList.getPatternCount(); j++) {
                reconstructedStates[node.getNumber()][j] = drawChoice(liks);
            }
        } else {
            //external node, currently not allowing ambiguities
            double[] liks = nodeProbDist[node.getNumber()].endLineageProbs;
            /*System.out.println("length: " + liks.length);
            for (int i = 0; i < liks.length; i++) {
                System.out.println(liks[i]);
            }*/
            for (int j = 0; j < patternList.getPatternCount(); j++) {
                reconstructedStates[node.getNumber()][j] = drawChoice(liks);
            }
        }

        //System.exit(0);

        NodeRef child1 = tree.getChild(node, 0);
        if (child1 != null) {
            traverseSample(tree, child1);
        }

        NodeRef child2 = tree.getChild(node, 1);
        if (child2 != null) {
            traverseSample(tree, child2);
        }

    }

    private int drawChoice(double[] measure) {
        if (useMAP) {
            double max = measure[0];
            int choice = 0;
            for (int i = 1; i < measure.length; i++) {
                if (measure[i] > max) {
                    max = measure[i];
                    choice = i;
                }
            }
            return choice;
        } else {
            return MathUtils.randomChoicePDF(measure);
        }
    }

    @Override
    public TreeModel getTreeModel() {
        return this.treeModel;
    }

    protected TreeTraitProvider.Helper treeTraits = new TreeTraitProvider.Helper();

    public TreeTrait getTreeTrait(String key) {
        return treeTraits.getTreeTrait(key);
    }

    public String formattedState(int[] state) {
        return formattedState(state, dataType);
    }

    private static String formattedState(int[] state, DataType dataType) {
        StringBuffer sb = new StringBuffer();
        sb.append("\"");
        if (dataType instanceof GeneralDataType) {
            boolean first = true;
            for (int i : state) {
                if (!first) {
                    sb.append(" ");
                } else {
                    first = false;
                }

                sb.append(dataType.getCode(i));
            }
        } else {
            throw new RuntimeException("Only GeneralDataType currently accepted.");
            /*for (int i : state) {
                if (dataType instanceof Codons) {
                    sb.append(dataType.getTriplet(i));
                } else {
                    sb.append(dataType.getChar(i));
                }
            }*/
        }
        sb.append("\"");
        return sb.toString();
    }

    public int[] getStatesForNode(Tree tree, NodeRef node) {
        if (tree != treeModel) {
            throw new RuntimeException("Can only reconstruct states on treeModel given to constructor");
        }

        if (!likelihoodKnown) {
            calculateLogLikelihood();
            likelihoodKnown = true;
        }

        if (!areStatesRedrawn) {
            redrawAncestralStates();
        }
        return reconstructedStates[node.getNumber()];
    }

    //TreeTraitProvider interface
    public TreeTrait[] getTreeTraits() {
        return treeTraits.getTreeTraits();
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

        /*public void update(double distance, NodeRef node, IntervalType intervalType) {
            this.intervalLength = distance;
            this.node = node;
            this.intervalType = intervalType;
            this.leftChild = null;
            this.rightChild = null;
            this.incremented = false;
            for (int i = 0; i < startLineageProbs.length; i++) {
                this.startLineageProbs[i] = 0.0;
                this.startLineageProbsSquared[i] = 0.0;
                this.endLineageProbs[i] = 0.0;
                this.endLineageProbsSquared[i] = 0.0;
            }
            this.startLineageProbs[this.patternIndex] = 1.0;
            this.startLineageProbsSquared[this.patternIndex] = 1.0;
            this.endLineageProbs[this.patternIndex] = 1.0;
            this.endLineageProbsSquared[this.patternIndex] = 1.0;
        }*/

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

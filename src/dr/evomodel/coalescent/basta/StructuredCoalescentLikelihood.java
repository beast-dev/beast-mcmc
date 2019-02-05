/*
 * StructuredCoalescentLikelihood.java
 *
 * Copyright (c) 2002-2019 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.TaxonList;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.coalescent.AbstractCoalescentLikelihood;
import dr.evomodel.substmodel.GeneralSubstitutionModel;
import dr.evomodel.tree.TreeChangedEvent;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.*;

import java.util.*;

/**
 * @author Guy Baele
 *
 * Implementation of BASTA: Bayesian structured coalescent approximation.
 * Original paper: Nicola De Maio, Chieh-Hsi Wu, Kathleen O'Reilly and Daniel Wilson
 * "New routes to phylogeography: a Bayesian structured coalescent approximation".
 * PLOS Genetics 11, e1005421; doi: 10.1371/journal.pgen.1005421
 */
public class StructuredCoalescentLikelihood extends AbstractCoalescentLikelihood implements Citable {

    private static final boolean DEBUG = false;
    private static final boolean MATRIX_DEBUG = false;

    //private static final boolean USE_BEAGLE = false;
    private static final boolean ASSOC_MULTIPLICATION = true;

    public StructuredCoalescentLikelihood(Tree tree, BranchRateModel branchRateModel, Parameter popSizes, PatternList patternList,
                                          GeneralSubstitutionModel generalSubstitutionModel, int subIntervals,
                                          TaxonList includeSubtree, List<TaxonList> excludeSubtrees) throws TreeUtils.MissingTaxonException {

        super(StructuredCoalescentLikelihoodParser.STRUCTURED_COALESCENT, tree, includeSubtree, excludeSubtrees);

        this.treeModel = (TreeModel)tree;
        this.patternList = patternList;

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
        this.subIntervals = subIntervals;

        this.nodeProbDist = new HashMap<NodeRef,ProbDist>();
        this.activeLineageList = new ArrayList<ProbDist>();
        this.tempLineageList = new ArrayList<ProbDist>();

        this.maxCoalescentIntervals = treeModel.getTaxonCount() * 2 - 2;
        //System.out.println("maxCoalescentIntervals = " + maxCoalescentIntervals);
        this.currentCoalescentInterval = 0;
        this.migrationMatrices = new double[maxCoalescentIntervals][this.demes*this.demes];
        this.storedMigrationMatrices = new double[maxCoalescentIntervals][this.demes*this.demes];

    }

    // **************************************************************
    // Likelihood IMPLEMENTATION
    // **************************************************************

    /**
     * Calculates the log likelihood of this set of coalescent intervals,
     * given a demographic model.
     */
    public double calculateLogLikelihood() {

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
        return logLikelihood;
    }

    public double getLogLikelihood() {
        if (!likelihoodKnown) {
            logLikelihood = calculateLogLikelihood();
            likelihoodKnown = true;
        }
        return logLikelihood;
    }

    //based on the traverseTree method in OldAbstractCoalescentLikelihood
    private double traverseTree(Tree tree, NodeRef root, PatternList patternList) {
        double MULTIFURCATION_LIMIT = 1e-9;

        ArrayList<ComparableDouble> times = new ArrayList<ComparableDouble>();
        ArrayList<Integer> childs = new ArrayList<Integer>();
        ArrayList<NodeRef> nodes = new ArrayList<NodeRef>();
        collectAllTimes(tree, root, nodes, times, childs);
        int[] indices = new int[times.size()];

        HeapSort.sort(times, indices);

        double lnL = 0.0;

        // start is the time of the first tip
        double start = times.get(indices[0]).doubleValue();

        int i = 0;
        int j = 0;
        //System.out.println("** (j=" + j + ") \n");
        while (i < times.size()) {

            int lineagesRemoved = 0;
            int lineagesAdded = 0;

            final double finish = times.get(indices[i]).doubleValue();
            //int j = i;
            double next = finish;

            if (DEBUG) {
                System.out.println("start = " + start + " ; finish = " + finish);
            }

            double intervalLength = finish - start;

            while (Math.abs(next - finish) < MULTIFURCATION_LIMIT) {
                final int children = childs.get(indices[i]);
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
                        if (treeModel.isExternal(nodes.get(indices[j]))) {
                            ProbDist newProbDist = new ProbDist(demes, 0.0, nodes.get(indices[j]));
                            newProbDist.setIntervalType(IntervalType.SAMPLE);
                            newProbDist.setStartLineageProb(patternList.getPattern(0)[patternList.getTaxonIndex(treeModel.getNodeTaxon(newProbDist.node).getId())], 1.0);
                            newProbDist.computeEndLineageDensities(0.0, null);
                            if (DEBUG) {
                                System.out.println(newProbDist);
                            }
                            tempLineageList.add(newProbDist);
                            //activeLineageList.add(newProbDist);
                            //System.out.println("currently active lineages: " + activeLineageList.size());
                            //coalescentIntervalList.get(intervalCount).add(newProbDist);
                            //add to HashMap
                            nodeProbDist.put(nodes.get(indices[j]),newProbDist);
                        } else {
                            if (treeModel.getChildCount(nodes.get(indices[j])) > 2) {
                                throw new RuntimeException("Structured coalescent currently only allows strictly bifurcating trees.");
                            }
                            NodeRef leftChild = treeModel.getChild(nodes.get(indices[j]), 0);
                            NodeRef rightChild = treeModel.getChild(nodes.get(indices[j]), 1);
                            ProbDist newProbDist = new ProbDist(demes, intervalLength, nodes.get(indices[j]), leftChild, rightChild);
                            newProbDist.setIntervalType(IntervalType.COALESCENT);
                            tempLineageList.add(newProbDist);
                            //activeLineageList.add(newProbDist);
                            //System.out.println("currently active lineages: " + activeLineageList.size());
                            //coalescentIntervalList.get(intervalCount).add(newProbDist);
                            //add to HashMap
                            //nodeProbDist.put(nodes.get(indices[j]), newProbDist);
                        }
                        //System.out.println("** (j=" + j + ") " + treeModel.getNodeHeight(nodes.get(indices[j])));
                        j++;
                        //System.out.println("** (j=" + j + ") ");
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
                        if (treeModel.isExternal(nodes.get(indices[j]))) {
                            ProbDist newProbDist = new ProbDist(demes, 0.0, nodes.get(indices[j]));
                            newProbDist.setIntervalType(IntervalType.SAMPLE);
                            newProbDist.setStartLineageProb(patternList.getPattern(0)[patternList.getTaxonIndex(treeModel.getNodeTaxon(newProbDist.node).getId())], 1.0);
                            newProbDist.computeEndLineageDensities(0.0, null);
                            if (DEBUG) {
                                System.out.println(newProbDist);
                            }
                            tempLineageList.add(newProbDist);
                            //activeLineageList.add(newProbDist);
                            //System.out.println("currently active lineages: " + activeLineageList.size());
                            //coalescentIntervalList.get(intervalCount).add(newProbDist);
                            //add to HashMap
                            nodeProbDist.put(nodes.get(indices[j]),newProbDist);
                        } else {
                            throw new RuntimeException("First interval cannot be a coalescent event.");
                            /*if (treeModel.getChildCount(nodes.get(indices[j])) > 2) {
                                throw new RuntimeException("Structured coalescent currently only allows strictly bifurcating trees.");
                            }
                            NodeRef leftChild = treeModel.getChild(nodes.get(indices[j]), 0);
                            NodeRef rightChild = treeModel.getChild(nodes.get(indices[j]), 1);
                            ProbDist newProbDist = new ProbDist(demes, finish-start, nodes.get(indices[j]), leftChild, rightChild);*/
                            //coalescentIntervalList.get(intervalCount).add(newProbDist);
                            //add to HashMap
                            //nodeProbDist.put(nodes.get(indices[j]), newProbDist);
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
                    if (treeModel.isExternal(nodes.get(indices[j]))) {
                        ProbDist newProbDist = new ProbDist(demes, intervalLength, nodes.get(indices[j]));
                        newProbDist.setIntervalType(IntervalType.SAMPLE);
                        newProbDist.setStartLineageProb(patternList.getPattern(0)[patternList.getTaxonIndex(treeModel.getNodeTaxon(newProbDist.node).getId())], 1.0);
                        newProbDist.computeEndLineageDensities(0.0, null);
                        tempLineageList.add(newProbDist);
                        //activeLineageList.add(newProbDist);
                        //System.out.println("currently active lineages: " + activeLineageList.size());
                        //coalescentIntervalList.get(intervalCount).add(newProbDist);
                        //add to HashMap
                        nodeProbDist.put(nodes.get(indices[j]),newProbDist);
                    } else {
                        if (treeModel.getChildCount(nodes.get(indices[j])) > 2) {
                            throw new RuntimeException("Structured coalescent currently only allows strictly bifurcating trees.");
                        }
                        NodeRef leftChild = treeModel.getChild(nodes.get(indices[j]), 0);
                        NodeRef rightChild = treeModel.getChild(nodes.get(indices[j]), 1);
                        ProbDist leftProbDist = nodeProbDist.get(leftChild);
                        ProbDist rightProbDist = nodeProbDist.get(rightChild);
                        tempLineageList.remove(leftProbDist);
                        tempLineageList.remove(rightProbDist);
                        if (DEBUG) {
                            System.out.println("currently active lineages: " + activeLineageList.size());
                        }
                        ProbDist newProbDist = new ProbDist(demes, intervalLength, nodes.get(indices[j]), leftChild, rightChild);
                        newProbDist.setIntervalType(IntervalType.COALESCENT);
                        lnL += newProbDist.computeCoalescedLineage(leftProbDist, rightProbDist);
                        if (!treeModel.isRoot(nodes.get(indices[j]))) {
                            tempLineageList.add(newProbDist);
                            //activeLineageList.add(newProbDist);
                            //System.out.println("currently active lineages: " + activeLineageList.size());
                            //coalescentIntervalList.get(intervalCount).add(newProbDist);
                            //add to HashMap
                            nodeProbDist.put(nodes.get(indices[j]), newProbDist);
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
                    System.out.println("Computing expected lineage counts");
                }
                computeExpectedLineageCounts();
            }

            //and then compute the log likelihood
            if (!(finish == 0.0)) {
                lnL += computeLogLikelihood(intervalLength);
                if (DEBUG) {
                    System.out.println("Computing (log) likelihood contribution; logP = " + lnL);
                }
            }

            //update the list of active lineages for next iteration
            activeLineageList.clear();
            for (ProbDist pd : tempLineageList) {
                activeLineageList.add(pd);
            }
            //tempLineageList.clear();
            if (DEBUG) {
                System.out.println("currently active lineages: " + activeLineageList.size() + "\n");
            }

        }

        if (DEBUG) {
            System.out.println("Structured coalescent lnL = " + lnL);
        }

        this.finalCoalescentInterval = this.currentCoalescentInterval;
        //System.out.println("finalCoalescentInterval = " + finalCoalescentInterval);
        this.currentCoalescentInterval = 0;
        this.matricesKnown = true;

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
            System.out.println("active lineage list length = " + activeLineageList.size());
        }

        for (int i = 0; i < demes; i++) {
            double startProbs = 0.0;
            double endProbs = 0.0;
            for (ProbDist pd : activeLineageList) {
                startProbs += pd.getStartLineageProb(i)*pd.getStartLineageProb(i);
                endProbs += pd.getEndLineageProb(i)*pd.getEndLineageProb(i);
            }
            if (DEBUG) {
                System.out.println("  startProbs = " + startProbs + " ; endProbs = " + endProbs);
            }
            intervalOne += (1.0/popSizes.getParameterValue(i)) * (startExpected[i]*startExpected[i] - startProbs);
            intervalTwo += (1.0/popSizes.getParameterValue(i)) * (endExpected[i]*endExpected[i] - endProbs);
        }
        intervalOne *= -intervalLength/4.0;
        intervalTwo *= -intervalLength/4.0;
        if (DEBUG) {
            System.out.println("interval 1 (log) likelihood = " + intervalOne);
            System.out.println("interval 2 (log) likelihood = " + intervalTwo);
            System.out.println("total (log) likelihood = " + (intervalOne + intervalTwo));
        }

        return intervalOne + intervalTwo;
    }

    /**
     * Iterate over all the currently active lineages and compute the expected lineage counts.
     */
    private void computeExpectedLineageCounts() {
        for (int i = 0; i < this.demes; i++) {
            this.startExpected[i] = 0.0;
            this.endExpected[i] = 0.0;
            for (ProbDist pd : this.activeLineageList) {
                this.startExpected[i] += pd.getStartLineageProb(i);
                this.endExpected[i] += pd.getEndLineageProb(i);
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

            if (MATRIX_DEBUG) {
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

        for (ProbDist pd : activeLineageList) {
            pd.incrementIntervalLength(increment, migrationMatrices[this.currentCoalescentInterval]);
            if (DEBUG) {
                System.out.println("  " + pd);
            }
        }
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

    // **************************************************************
    // ModelListener IMPLEMENTATION
    // **************************************************************

    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (DEBUG) {
            System.out.println("handleModelChangedEvent: " + model.getModelName() + ", " + object + ", " + index);
        }
        if (model == treeModel) {
            //for all the nodes that are older than the event, set updateProbDist (still to implement) to true
            //then trigger a recalculation that makes use of an adjusted traverseTree method (that checks whether
            //or not the ProbDist needs to be updated
            likelihoodKnown = false;
            matricesKnown = false;
        } else if (model == branchRateModel) {
            likelihoodKnown = false;
            matricesKnown = false;
        } else if (model == generalSubstitutionModel) {
            likelihoodKnown = false;
            matricesKnown = false;
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
        //a change in one of the popSizes does not affect matrix exponentiation
        matricesKnown = true;
    }

    protected void storeState() {
        //System.out.println("STORESTATE, matricesKnown=" + this.matricesKnown);
        //super.storeState();
        for (int i = 0; i < this.finalCoalescentInterval; i++) {
            System.arraycopy(this.migrationMatrices[i],0, this.storedMigrationMatrices[i], 0, demes*demes);
        }
        /*double[] tmp;
        for (int i = 0; i < this.finalCoalescentInterval; i++) {
            tmp = this.storedMigrationMatrices[i];
            this.storedMigrationMatrices[i] = this.migrationMatrices[i];
            this.migrationMatrices[i] = tmp;
        }*/
        /*if (DEBUG) {
            System.out.println("STORESTATE, matricesKnown=" + this.matricesKnown);
        }*/
        storedLikelihoodKnown = likelihoodKnown;
        storedLogLikelihood = logLikelihood;
    }

    protected void restoreState() {
        //System.out.println("RESTORESTATE, matricesKnown=" + this.matricesKnown);
        //super.restoreState();
        double[] tmp;
        for (int i = 0; i < this.finalCoalescentInterval; i++) {
            tmp = this.migrationMatrices[i];
            this.migrationMatrices[i] = this.storedMigrationMatrices[i];
            this.storedMigrationMatrices[i] = tmp;
        }
        /*if (DEBUG) {
            System.out.println("RESTORESTATE, matricesKnown=" + this.matricesKnown);
        }*/
        likelihoodKnown = storedLikelihoodKnown;
        logLikelihood = storedLogLikelihood;
    }

    public void makeDirty() {
        likelihoodKnown = false;
        matricesKnown = false;
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

        //keep track of whether the interval length has already been incremented for matrix exponentiation
        private boolean incremented = false;

        //coalescent or sample
        private IntervalType intervalType;

        //child nodes only server a purpose when dealing with a coalescent event
        private NodeRef leftChild = null;
        private NodeRef rightChild = null;

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

        //compute the probability distribution of lineages among demes for a coalescent event
        public double computeCoalescedLineage(ProbDist leftProbDist, ProbDist rightProbDist) {
            double sum = 0.0;
            double[] sumComponents = new double[demes];
            for (int i = 0; i < demes; i++) {
                sumComponents[i] = (leftProbDist.getEndLineageProb(i) * rightProbDist.getEndLineageProb(i))/popSizes.getParameterValue(i);
                sum += sumComponents[i];
            }

            if (DEBUG) {
                System.out.println("coalescent end lineage prob 0 = " + leftProbDist.getEndLineageProb(0));
                System.out.println("coalescent end lineage prob 1 = " + leftProbDist.getEndLineageProb(1));
            }
            for (int i = 0; i < demes; i++) {
                this.setStartLineageProb(i, sumComponents[i]/sum);
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
            if (lineageLength == 0.0) {
                for (int k = 0; k < demes; k++) {
                    this.setEndLineageProb(k, this.getStartLineageProb(k));
                }
            } else {
                //TODO this should be possible in parallel for each lineage within the same coalescent interval
                for (int k = 0; k < demes; k++) {
                    double value = 0.0;
                    for (int l = 0; l < demes; l++) {
                        value += this.startLineageProbs[l] * migrationMatrix[l*demes+k];
                    }
                    this.setEndLineageProb(k, value);
                }
            }

        }

        public void incrementIntervalLength(double increment, double[] migrationMatrix) {
            this.intervalLength = this.intervalLength + increment;
            if (incremented) {
                if (ASSOC_MULTIPLICATION) {
                for (int i = 0; i < demes; i++) {
                    startLineageProbs[i] = endLineageProbs[i];
                }
                computeEndLineageDensities(increment, migrationMatrix);
                } else {
                    throw new RuntimeException("Only incremental matrix exponentiation allowed for performance reasons.");
                }
            } else {
                computeEndLineageDensities(this.intervalLength, migrationMatrix);
            }
            this.incremented = true;
        }

        public IntervalType getIntervalType() {
            return this.intervalType;
        }

        public void setIntervalType(IntervalType intervalType) {
            this.intervalType = intervalType;
        }

        public NodeRef getLeftChild() {
            return this.leftChild;
        }

        public NodeRef getRightChild() {
            return this.rightChild;
        }

        public double getStartLineageProb(int index) {
            return this.startLineageProbs[index];
        }

        public double getEndLineageProb(int index) {
            return this.endLineageProbs[index];
        }

        public void setStartLineageProb(int index, double value) {
            this.startLineageProbs[index] = value;
        }

        private void setEndLineageProb(int index, double value) {
            this.endLineageProbs[index] = value;
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

    //the tree model
    private TreeModel treeModel;

    //the branch rate model
    private BranchRateModel branchRateModel;

    //the population sizes
    private Parameter popSizes;

    //the discrete trait data
    private PatternList patternList;

    //expected starting lineage counts
    private double[] startExpected;

    //expected ending lineage counts
    private double[] endExpected;

    //map with a probability distribution of lineages for each node
    private HashMap<NodeRef,ProbDist> nodeProbDist;

    //list of currently active lineages
    private ArrayList<ProbDist> activeLineageList;

    //temporary list of lineagess
    private ArrayList<ProbDist> tempLineageList;

    //the migration model
    private GeneralSubstitutionModel generalSubstitutionModel;

    //number of demes for the structured coalescent model
    private int demes;

    //number of subintervals across each branch of the tree
    private int subIntervals;

    //variables that allow to use storeState and restoreState
    private int maxCoalescentIntervals;
    private int currentCoalescentInterval;
    private double[][] migrationMatrices;
    private int finalCoalescentInterval;
    private double[][] storedMigrationMatrices;

    private boolean matricesKnown;

}

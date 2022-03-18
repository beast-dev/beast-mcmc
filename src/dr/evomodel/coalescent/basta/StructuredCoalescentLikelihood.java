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

import java.util.*;

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

        this.maxCoalescentIntervals = treeModel.getTaxonCount() * 2 - 2;
        this.currentCoalescentInterval = 0;
        this.migrationMatrices = new double[maxCoalescentIntervals][this.demes*this.demes];
        this.storedMigrationMatrices = new double[maxCoalescentIntervals][this.demes*this.demes];

        this.intervalStartProbs = new double[intervals.getIntervalCount()][];
        this.intervalEndProbs = new double[intervals.getIntervalCount()][];
        this.intervalStartSquareProbs = new double[intervals.getIntervalCount()][];
        this.intervalEndSquareProbs = new double[intervals.getIntervalCount()][];
        this.coalescentLeftProbs = new double[intervals.getIntervalCount()][];
        this.coalescentRightProbs = new double[intervals.getIntervalCount()][];
        for (int i = 0; i < intervals.getIntervalCount(); i++) {
            this.intervalStartProbs[i] = new double[demes];
            this.intervalEndProbs[i] = new double[demes];
            this.intervalStartSquareProbs[i] = new double[demes];
            this.intervalEndSquareProbs[i] = new double[demes];
            this.coalescentLeftProbs[i] = new double[nodeCount*demes];
            this.coalescentRightProbs[i] = new double[nodeCount*demes];
        }
        //keep track of currently active lineage probabilities and use to compute / populate the lists above
        //this.activeLineages = new ArrayList<double[]>();
        this.activeLineages = new double[treeModel.getNodeCount()*demes];
        this.activeNodeNumbers = new ArrayList<Integer>();

        this.addedLineages = new boolean[nodeCount];
        for (int i = 0; i < addedLineages.length; i++) {
            addedLineages[i] = false;
        }
        this.addedLength = 0;

        this.times = new ArrayList<ComparableDouble>();
        this.children = new ArrayList<Integer>();
        this.nodes = new ArrayList<NodeRef>();

        this.storedTimes = new ArrayList<ComparableDouble>();
        this.storedChildren = new ArrayList<Integer>();
        this.storedNodes = new ArrayList<NodeRef>();

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
     *
     * @param startingInterval
     */
    public double calculateLogLikelihood(int startingInterval) {

        System.out.println("\n>calculateLogLikelihood");

        double logL = 0.0;

        //TODO only iterate over a subset of the intervals
        int intervalCount = intervals.getIntervalCount();
        for (int i = 0; i < intervalCount; i++) {
            double intervalLength = intervals.getInterval(i);

            //TODO this could be funky for zero-length intervals
            if (intervalLength != 0.0) {
                //get all the information from the various lists
                double[] lineageStartCount = intervalStartProbs[i];
                double[] lineageStartCountSquare = intervalStartSquareProbs[i];
                double[] lineageEndCount = intervalEndProbs[i];
                double[] lineageEndCountSquare = intervalEndSquareProbs[i];

                double halfLength = intervalLength/2.0;
                if (halfLength != 0.0) {
                    double densityOne = 0.0;
                    double densityTwo = 0.0;
                    //TODO I believe there is a calculation problem here
                    for (int j = 0; j < demes; j++) {
                        densityOne += (((lineageStartCount[j]*lineageStartCount[j])-lineageStartCountSquare[j])/(2.0*popSizes.getParameterValue(j)));
                        densityTwo += (((lineageEndCount[j]*lineageEndCount[j])-lineageEndCountSquare[j])/(2.0*popSizes.getParameterValue(j)));

                        System.out.println("lineageStartCount[j] = " + lineageStartCount[j]);
                        System.out.println("lineageStartCountSquare[j] = " + lineageStartCountSquare[j]);
                        System.out.println("lineageEndCount[j] = " + lineageEndCount[j]);
                        System.out.println("lineageEndCountSquare[j] = " + lineageEndCountSquare[j]);
                    }
                    logL += -halfLength*densityOne;
                    System.out.println("logL = " + logL);
                    logL += -halfLength*densityTwo;
                    System.out.println("logL = " + logL);
                }

                if (intervals.getIntervalType(i) == IntervalType.COALESCENT) {
                    System.out.println("coalescent contribution");
                    double contribution = 0.0;
                    for (int j = 0; j < demes; j++) {
                        contribution += (coalescentLeftProbs[i][j]*coalescentRightProbs[i][j])/ popSizes.getParameterValue(j);
                        System.out.println("coalescentLeftProbs[i][j] = " + coalescentLeftProbs[i][j]);
                        System.out.println("coalescentRightProbs[i][j] = " + coalescentRightProbs[i][j]);
                    }
                    logL += Math.log(contribution);
                    System.out.println("logL = " + logL);
                } else {
                    //do nothing
                }

            } else {
                //do nothing
            }

            System.out.println("logL(interval " + i + ") = " + logL);

        }

        System.out.println("final logL = " + logL);
        return logL;

    }

    /**
     * Compute all required probability distributions for calculating the overall structured coalescent density.
     * This methods computes equations 11 and 12 from the BASTA manuscript and populates the following arrays
     * in order to enable a straight forward structured coalescent density computation in calculateLogLikelihood:
     * intervalStartProbs, intervalStartSquareProbs, intervalEndProbs and intervalEndSquareProbs, and for the
     * contributions at the coalescent events: coalescentLeftProbs and coalescentRightProbs
     */
    private void computeProbabilityDistributions() {

        System.out.println("\n>computeProbabilityDistributions");

        int intervalCount = intervals.getIntervalCount();
        System.out.println("number of intervals: " + intervalCount);

        //first sampling event is not considered so take this into account
        NodeRef node = intervals.getSamplingNode(-1);
        System.out.println("sampling node: "+ treeModel.getNodeTaxon(node).getId());
        /*double[] lineageCount = new double[demes];
        for (int i = 0; i < demes; i++) {
            lineageCount[i] = 0.0;
        }*/
        //set start probabilities of first sampling interval

        //lineageCount[patternList.getPattern(0)[patternList.getTaxonIndex(treeModel.getNodeTaxon(node).getId())]] = 1.0;
        this.activeLineages[node.getNumber()*demes+patternList.getPattern(0)[patternList.getTaxonIndex(treeModel.getNodeTaxon(node).getId())]] = 1.0;
        this.activeNodeNumbers.add(node.getNumber());

        //this.intervalStartProbs[0] = Arrays.copyOf(lineageCount, demes);
        //this.intervalStartSquareProbs[0] = Arrays.copyOf(lineageCount, demes);
        this.intervalStartProbs[0] = Arrays.copyOfRange(this.activeLineages, node.getNumber()*demes, node.getNumber()*demes+demes);
        this.intervalStartSquareProbs[0] = Arrays.copyOfRange(this.activeLineages, node.getNumber()*demes, node.getNumber()*demes+demes);
        //add to list of active lineages
        //this.activeLineages.add(lineageCount);

        //print content of active lineages
        printActiveLineages();

        //each interval is defined by its ending event, i.e. does the interval end in a sampling or coalescent event
        for (int i = 0; i < intervalCount; i++) {

            System.out.println("interval type: " + intervals.getIntervalType(i));

            if (intervals.getIntervalType(i) == IntervalType.COALESCENT) {

                //compute end probabilities of sampling interval
                double intervalLength = intervals.getInterval(i);

                //matrix exponentiation to compute end interval probabilities; equation 11
                incrementActiveLineages(this.activeLineages, intervalLength, i);

                printActiveLineages();

                //get coalescent node from interval
                node = intervals.getCoalescentNode(i);
                //get child nodes
                NodeRef leftChild = treeModel.getChild(node, 0);
                NodeRef rightChild = treeModel.getChild(node, 1);

                //access probability densities from both child nodes to compute equation 12
                //get end lineage densities for the 2 child nodes from the previous interval
                double[] temp = new double[demes];
                double sum = 0.0;
                for (int k = 0; k < demes; k++) {
                    temp[k] = (activeLineages[leftChild.getNumber()*demes+k]*activeLineages[rightChild.getNumber()*demes+k])/popSizes.getParameterValue(k);
                    sum += temp[k];
                }
                //store the resulting coalescent probability density
                for (int k = 0; k < demes; k++) {
                    this.activeLineages[node.getNumber()*demes+k] = temp[k]/sum;
                }

                //TODO add information to coalescentLeftProbs and coalescentRightProbs really necessary
                //TODO merge into one of the previous loops over k?
                //this code mostly to keep the calculateLogLikelihood function as clean as possible
                for (int k = 0; k < demes; k++) {
                    this.coalescentLeftProbs[i][k] = activeLineages[leftChild.getNumber()*demes+k];
                    this.coalescentRightProbs[i][k] = activeLineages[rightChild.getNumber()*demes+k];
                }

                //remove 2 nodes from active lineage list and add a new one
                this.activeNodeNumbers.remove((Integer)leftChild.getNumber());
                this.activeNodeNumbers.remove((Integer)rightChild.getNumber());
                this.activeNodeNumbers.add(node.getNumber());

                printActiveLineages();

            } else if (intervals.getIntervalType(i) == IntervalType.SAMPLE) {
                //check for zero-length interval
                if (intervals.getInterval(i) == 0.0) {
                    //multiple samples at same sampling time
                    System.out.println("zero-length interval");

                    node = intervals.getSamplingNode(i+1);
                    System.out.println("sampling node: "+ treeModel.getNodeTaxon(node).getId());
                    this.activeLineages[node.getNumber()*demes+patternList.getPattern(0)[patternList.getTaxonIndex(treeModel.getNodeTaxon(node).getId())]] = 1.0;
                    this.activeNodeNumbers.add(node.getNumber());

                    printActiveLineages();

                } else {
                    //one sample at a time

                    //compute end probabilities of sampling interval
                    double intervalLength = intervals.getInterval(i);

                    //first compute Sd for the first interval half
                    //TODO check these computations
                    for (int k = 0; k < demes; k++) {
                        this.intervalStartProbs[i][k] = 0.0;
                        this.intervalStartSquareProbs[i][k] = 0.0;
                        for (int l = 0; l < activeNodeNumbers.size(); l++) {
                            this.intervalStartProbs[i][k] += activeLineages[activeNodeNumbers.get(l)*demes+l];
                            this.intervalStartSquareProbs[i][k] += activeLineages[activeNodeNumbers.get(l)*demes+l]*activeLineages[activeNodeNumbers.get(l)*demes+l];
                        }
                    }

                    //matrix exponentiation to compute end interval probabilities; equation 11
                    incrementActiveLineages(this.activeLineages, intervalLength, i);

                    /*double[] lineageCountSquare = new double[demes];
                    for (int j = 0; j < demes; j++) {
                        lineageCountSquare[j] = lineageCount[j]*lineageCount[j];
                    }
                    //computing Sd (see supplementary)
                    double[] sumActiveLineages = new double[demes];
                    for (int j = 0; j < demes; j++) {
                        sumActiveLineages[j] = 0.0;
                    }
                    for (double[] lineage : activeLineages) {
                        for (int j = 0; j < demes; j++) {
                            sumActiveLineages[j] += lineage[j];
                        }
                    }*/

                    //compute Sd for the second interval half
                    //TODO check these computations
                    for (int k = 0; k < demes; k++) {
                        this.intervalEndProbs[i][k] = 0.0;
                        this.intervalEndSquareProbs[i][k] = 0.0;
                        for (int l = 0; l < activeNodeNumbers.size(); l++) {
                            this.intervalEndProbs[i][k] += activeLineages[activeNodeNumbers.get(l)*demes+l];
                            this.intervalEndSquareProbs[i][k] += activeLineages[activeNodeNumbers.get(l)*demes+l]*activeLineages[activeNodeNumbers.get(l)*demes+l];
                        }
                    }

                    printIntervalContributions(i);

                    //this.intervalEndProbs[i] = Arrays.copyOf(sumActiveLineages, demes);
                    //this.intervalEndSquareProbs[i] = Arrays.copyOf(lineageCountSquare, demes);

                    //get the node number of the sampling node
                    node = intervals.getSamplingNode(i);
                    System.out.println("sampling node: "+ treeModel.getNodeTaxon(node).getId());

                    //add to the active lineages for next interval iteration
                    /*double[] newLineageCount = new double[demes];
                    for (int k = 0; k < demes; k++) {
                        newLineageCount[k] = 0.0;
                    }*/

                    //set start probabilities of first sampling interval
                    this.activeLineages[node.getNumber()*demes+patternList.getPattern(0)[patternList.getTaxonIndex(treeModel.getNodeTaxon(node).getId())]] = 1.0;
                    this.activeNodeNumbers.add(node.getNumber());

                    //add to list of active lineages
                    //this.activeLineages.add(newLineageCount);

                    printActiveLineages();

                }
            }
        }

    }

    private void printIntervalContributions(int interval) {
        System.out.print("  starting lineage count: ");
        for (int k = 0; k < demes; k++) {
            System.out.print(this.intervalStartProbs[interval][k] + " ");
        }
        System.out.println();
        System.out.print("  starting lineage count squared: ");
        for (int k = 0; k < demes; k++) {
            System.out.print(this.intervalStartSquareProbs[interval][k] + " ");
        }
        System.out.println();
        System.out.print("  ending lineage count: ");
        for (int k = 0; k < demes; k++) {
            System.out.print(this.intervalEndProbs[interval][k] + " ");
        }
        System.out.println();
        System.out.print("  ending lineage count squared: ");
        for (int k = 0; k < demes; k++) {
            System.out.print(this.intervalEndSquareProbs[interval][k] + " ");
        }
        System.out.println();
    }

    private void printActiveLineages() {
        System.out.println("printActiveLineages");
        //print content of active lineages
        for (int number : this.activeNodeNumbers) {
            System.out.print("  active lineage: ");
            for (int k = 0; k < demes; k++) {
                System.out.print(this.activeLineages[number*demes+k] + " ");
            }
            System.out.println();
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
     * Compute all the probability densities (equation 11) but bypass using an ArrayList and work directly on a double
     * array by using an offset equal to nodeNumber times number of demes.
     * The results have to be copied after calling this method to be used after tree changed events.
     *
     * @param activeLineages padded array of probability densities for the current interval to be updated
     *                     (comes out of an arraylist of double arrays)
     * @param intervalLength lenght of the current coalescent interval over which to update using equation 11
     * @param interval denotes the index of the interval to get the correct matrix exponential
     */
    private void incrementActiveLineages(double[] activeLineages, double intervalLength, int interval) {

        System.out.println("incrementActiveLineages with intervalLength: " + intervalLength);

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

        //compute all dot products / probability densities (i.e. equation 11) for all active nodes
        for (int active : activeNodeNumbers) {
            double[] temp = new double[demes];
            for (int k = 0; k < demes; k++) {
                temp[k] = USE_TRANSPOSE ?
                        rdot(demes, activeLineages, active*demes, 1, migrationMatrices[interval], k * demes, 1) :
                        rdot(demes, activeLineages, active*demes, 1, migrationMatrices[interval], k, demes);
            }
            //copy back from temp
            for (int k = 0; k < demes; k++) {
                activeLineages[active*demes+k] = temp[k];
            }
        }
    }

    /**
     * Compute all the probability densities (equation 11)
     *
     * @param intervalLength
     * @param interval denotes the index of the interval to get the correct matrix exponential
     */
    private void incrementActiveLineages(ArrayList<double[]> lineageCount, double intervalLength, int interval) {

        System.out.println("incrementActiveLineages with intervalLength: " + intervalLength);

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

        for (double[] linCount : lineageCount) {
            //temporary array
            double[] temp = new double[linCount.length];
            for (int k = 0; k < demes; k++) {
                    /*double value = 0.0;
                    for (int l = 0; l < demes; l++) {
                        value += this.startLineageProbs[l] * migrationMatrix[l*demes+k];
                    }
                    this.endLineageProbs[k] = value;*/
                temp[k] = USE_TRANSPOSE ?
                        rdot(demes, linCount, 0, 1, migrationMatrices[interval], k * demes, 1) :
                        rdot(demes, linCount, 0, 1, migrationMatrices[interval], k, demes);
            }
            for (int k = 0; k < demes; k++) {
                linCount[k] = temp[k];
            }
        }

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
            areStatesRedrawn = false;
        } else if (model == branchRateModel) {
            matricesKnown = false;
            //the following to accommodate events stemming from the upDownOperator
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
        likelihoodKnown = false;
        matricesKnown = false;
        areStatesRedrawn = false;
    }

    public TreeModel getTreeModel() {
        return this.treeModel;
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

    private boolean[] addedLineages;
    private int addedLength;

    private BigFastTreeIntervals intervals;

    //probability densities at the start and end of each coalescent interval
    //first index is the number of the coalescent interval
    //second index is the deme number
    private double[][] intervalStartProbs;
    private double[][] intervalEndProbs;
    private double[][] intervalStartSquareProbs;
    private double[][] intervalEndSquareProbs;
    //keep track of information to compute structured coalescent density at coalescent events
    private double[][] coalescentLeftProbs;
    private double[][] coalescentRightProbs;

    //TODO exchange for one very large array instead of a list of arrays to store the probability densities in each interval?
    //TODO work with offsets to only use one array; offset: interval number, offset 2: node number
    //offset is the node number; only need one double array to keep track of active probability densities
    //private ArrayList<double[]> activeLineages;
    private double[] activeLineages;
    //also keep track of the active node numbers per interval
    private ArrayList<Integer> activeNodeNumbers;

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

}

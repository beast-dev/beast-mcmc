/*
 * StructuredCoalescentLikelihood.java
 *
 * Copyright (c) 2002-2022 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.evolution.datatype.DataType;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.evolution.util.TaxonList;
import dr.evolution.util.Units;
import dr.evomodel.bigfasttree.BestSignalsFromBigFastTreeIntervals;
import dr.evomodel.bigfasttree.IntervalChangedEvent;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.DefaultBranchRateModel;
import dr.evomodel.substmodel.GeneralSubstitutionModel;
import dr.evomodel.tree.TreeChangedEvent;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.AbstractModelLikelihood;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * @author Guy Baele
 * @author Marc A. Suchard
 *
 * Implementation of BASTA: Bayesian structured coalescent approximation.
 * Original paper: Nicola De Maio, Chieh-Hsi Wu, Kathleen O'Reilly and Daniel Wilson
 * "New routes to phylogeography: a Bayesian structured coalescent approximation".
 * PLOS Genetics 11, e1005421; doi: 10.1371/journal.pgen.1005421
 */
public class FasterStructuredCoalescentLikelihood extends AbstractModelLikelihood implements Units, Citable {

    private static final boolean DEBUG = false;
    private static final boolean MATRIX_DEBUG = false;
    private static final boolean UPDATE_DEBUG = false;

    //private static final boolean USE_BEAGLE = false;
    private static final boolean ASSOC_MULTIPLICATION = true;
    private static final boolean USE_TRANSPOSE = false;

    private static final boolean MINIMUM_EVALUATION = false;

    public FasterStructuredCoalescentLikelihood(Tree tree, BranchRateModel branchRateModel, Parameter popSizes, PatternList patternList,
                                                DataType dataType, String tag, GeneralSubstitutionModel generalSubstitutionModel, int subIntervals,
                                                TaxonList includeSubtree, List<TaxonList> excludeSubtrees, boolean useMAP) throws TreeUtils.MissingTaxonException {

        super(StructuredCoalescentLikelihoodParser.STRUCTURED_COALESCENT);

        this.treeModel = (TreeModel)tree;
        this.patternList = patternList;
        this.dataType = dataType;
        this.useMAP = useMAP;

        if (tree instanceof TreeModel) {
            //System.out.println("initial tree = " + (TreeModel) tree);
            //this.intervals = new BigFastTreeIntervals((TreeModel) tree);
            this.intervals = new BestSignalsFromBigFastTreeIntervals((TreeModel) tree);
            addModel(intervals);
            //addModel((TreeModel) tree);
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

        int nodeCount = treeModel.getNodeCount();
        this.intervalCount = intervals.getIntervalCount();

        this.matricesKnown = new boolean[intervalCount];
        this.storedMatricesKnown = new boolean[intervalCount];

        this.migrationMatrices = new double[intervalCount][this.demes*this.demes];
        this.storedMigrationMatrices = new double[intervalCount][this.demes*this.demes];

        this.intervalStartProbs = new double[intervalCount][];
        this.intervalEndProbs = new double[intervalCount][];
        this.intervalStartSquareProbs = new double[intervalCount][];
        this.intervalEndSquareProbs = new double[intervalCount][];
        this.coalescentLeftProbs = new double[intervalCount][];
        this.coalescentRightProbs = new double[intervalCount][];
        for (int i = 0; i < intervalCount; i++) {
            this.intervalStartProbs[i] = new double[demes];
            this.intervalEndProbs[i] = new double[demes];
            this.intervalStartSquareProbs[i] = new double[demes];
            this.intervalEndSquareProbs[i] = new double[demes];
            this.coalescentLeftProbs[i] = new double[nodeCount*demes];
            this.coalescentRightProbs[i] = new double[nodeCount*demes];
        }
        //keep track of currently active lineage probabilities and use to compute / populate the lists above
        //this.activeLineages = new ArrayList<double[]>();
        //TODO make this dependent on the interval to store all probability densities (?)
        this.activeLineages = new double[nodeCount*demes];
        this.activeNodeNumbers = new ArrayList<>(nodeCount);

        this.likelihoodKnown = false;

        this.temp = new double[demes];
    }

    final double[] temp;

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
        }
        return logLikelihood;
    }

    /**
     * Calculates the log structured coalescent density across all coalescent intervals.
     */
    public double calculateLogLikelihood() {

        //TODO write if-else clause that determines whether a full likelihood recalculation is required
        //TODO if not, write a method that only recomputes part of the probability distributions
        computeProbabilityDistributions(0);

        return calculateLogLikelihood(0);

    }

    /**
     * Calculates the log structured coalescent density from a given interval number up to the root.
     *
     * @param startingInterval the interval from which to start the calculations up to the root
     */
    public double calculateLogLikelihood(int startingInterval) {

        //System.out.println("\n>calculateLogLikelihood");

        double logL = 0.0;

        //TODO only iterate over a subset of the intervals
        for (int i = 0; i < intervalCount; i++) {

            //System.out.println("interval: " + i + " (" + intervals.getIntervalType(i) + ")");

            double intervalLength = intervals.getInterval(i);

            if (intervalLength != 0.0) {
                //get all the information from the various lists
                double[] lineageStartCount = intervalStartProbs[i];
                double[] lineageStartCountSquare = intervalStartSquareProbs[i];
                double[] lineageEndCount = intervalEndProbs[i];
                double[] lineageEndCountSquare = intervalEndSquareProbs[i];

                double halfLength = intervalLength/2.0;
                //System.out.println("half interval length = " + halfLength);
                if (halfLength != 0.0) {
                    double densityOne = 0.0;
                    double densityTwo = 0.0;
                    for (int j = 0; j < demes; j++) {
                        densityOne += (((lineageStartCount[j]*lineageStartCount[j])-lineageStartCountSquare[j])/(2.0*popSizes.getParameterValue(j)));
                        densityOne += (((lineageEndCount[j]*lineageEndCount[j])-lineageEndCountSquare[j])/(2.0*popSizes.getParameterValue(j)));

                        //System.out.println("lineageStartCount[" + j + "] = " + lineageStartCount[j]);
                        //System.out.println("lineageStartCountSquare[" + j + "] = " + lineageStartCountSquare[j]);
                        //System.out.println("lineageEndCount[" + j + "] = " + lineageEndCount[j]);
                        //System.out.println("lineageEndCountSquare[" + j +"] = " + lineageEndCountSquare[j]);
                    }
                    logL += -halfLength*densityOne;
                    //System.out.println("logL = " + logL);
                    logL += -halfLength*densityTwo;
                    //System.out.println("logL = " + logL);
                }

                if (intervals.getIntervalType(i) == IntervalType.COALESCENT) {
                    //System.out.println("coalescent contribution");
                    double contribution = 0.0;
//                    for (int j = 0; j < demes; j++) {
//                        contribution += (coalescentLeftProbs[i][j]*coalescentRightProbs[i][j])/ popSizes.getParameterValue(j);
//                        //System.out.println("coalescentLeftProbs[i][j] = " + coalescentLeftProbs[i][j]);
//                        //System.out.println("coalescentRightProbs[i][j] = " + coalescentRightProbs[i][j]);
//                    }
                    contribution = coalescentLeftProbs[i][0]; // Use cached value
                    logL += Math.log(contribution);
                    //System.out.println("logL = " + logL);
                } else {
                    //do nothing
                }

            } else {
                //do nothing
            }

            //System.out.println("logL(interval " + i + ") = " + logL);

        }

        //System.out.println("final logL = " + logL);
        return logL;

    }

    private void newHardWork(double[] probs, double[] squareProbs) {
//        for (int k = 0; k < demes; k++) {
//            probs[k] = 0.0;
//            squareProbs[k] = 0.0;
//            for (int l : activeNodeNumbers) {
//                probs[k] += activeLineages[l*demes+k];
//                squareProbs[k] += activeLineages[l*demes+k]*activeLineages[l*demes+k];
//            }
//        }
        for (int k = 0; k < demes; k++) {
            probs[k] = 0.0;
            squareProbs[k] = 0.0;
        }
        for (int l : activeNodeNumbers) {
            for (int k = 0; k < demes; k++) {
                probs[k] += activeLineages[l*demes+k];
                squareProbs[k] += activeLineages[l*demes+k]*activeLineages[l*demes+k];
            }
        }
    }

    private void handleCoalescense(NodeRef node, int offset, int i) {
                //compute end probabilities of sampling interval
                double intervalLength = intervals.getInterval(i);

                //printActiveLineages();

                //first compute Sd for the first interval half
//                for (int k = 0; k < demes; k++) {
//                    this.intervalStartProbs[i][k] = 0.0;
//                    this.intervalStartSquareProbs[i][k] = 0.0;
//                    for (int l = 0; l < activeNodeNumbers.size(); l++) {
//                        this.intervalStartProbs[i][k] += activeLineages[activeNodeNumbers.get(l)*demes+k];
//                        this.intervalStartSquareProbs[i][k] += activeLineages[activeNodeNumbers.get(l)*demes+k]*activeLineages[activeNodeNumbers.get(l)*demes+k];
//                    }
//                }
                newHardWork(this.intervalStartProbs[i], this.intervalStartSquareProbs[i]);
//                hardWork(i);

                //matrix exponentiation to compute end interval probabilities; equation 11
                incrementActiveLineages(this.activeLineages, intervalLength, i);

                //printActiveLineages();

                //get coalescent node from interval
                node = intervals.getCoalescentNode(i);

                offset = node.getNumber()*demes;
                //get child nodes
                NodeRef leftChild = treeModel.getChild(node, 0);
                NodeRef rightChild = treeModel.getChild(node, 1);

                final int leftOffset = leftChild.getNumber() * demes;
                final int rightOffset = rightChild.getNumber() * demes;

                //access probability densities from both child nodes to compute equation 12
                //get end lineage densities for the 2 child nodes from the previous interval
//                double[] temp = new double[demes];
                double sum = 0.0;
                for (int k = 0; k < demes; k++) {
                    temp[k] = (activeLineages[leftOffset+k]*activeLineages[rightOffset+k])/popSizes.getParameterValue(k);
                    sum += temp[k];
                }
                this.coalescentLeftProbs[i][0] = sum; // Cache this value directly, instead of recomputing later

                //store the resulting coalescent probability density
                for (int k = 0; k < demes; k++) {
                    this.activeLineages[offset+k] = temp[k] / sum;
                }

                //compute Sd for the second interval half
//                for (int k = 0; k < demes; k++) {
//                    this.intervalEndProbs[i][k] = 0.0;
//                    this.intervalEndSquareProbs[i][k] = 0.0;
//                    for (int l = 0; l < activeNodeNumbers.size(); l++) {
//                        this.intervalEndProbs[i][k] += activeLineages[activeNodeNumbers.get(l)*demes+k];
//                        this.intervalEndSquareProbs[i][k] += activeLineages[activeNodeNumbers.get(l)*demes+k]*activeLineages[activeNodeNumbers.get(l)*demes+k];
//                    }
//                }
                    newHardWork(this.intervalEndProbs[i], this.intervalEndSquareProbs[i]);

                //TODO merge into one of the previous loops over k?
                //this code mostly to keep the calculateLogLikelihood function as clean as possible
//                for (int k = 0; k < demes; k++) {
//                    this.coalescentLeftProbs[i][k] = activeLineages[leftChild.getNumber()*demes+k];
//                    this.coalescentRightProbs[i][k] = activeLineages[rightChild.getNumber()*demes+k];
//                }

                doShit(node, leftChild, rightChild);
//
//                //remove 2 nodes from active lineage list and add a new one
//                this.activeNodeNumbers.remove((Integer)leftChild.getNumber());
//                this.activeNodeNumbers.remove((Integer)rightChild.getNumber());
//                this.activeNodeNumbers.add(node.getNumber());

                //printActiveLineages();
    }

    private void doShit(NodeRef node, NodeRef leftChild, NodeRef rightChild) {
            //remove 2 nodes from active lineage list and add a new one
            this.activeNodeNumbers.remove((Integer)leftChild.getNumber());
            this.activeNodeNumbers.remove((Integer)rightChild.getNumber());
            this.activeNodeNumbers.add(node.getNumber());
    }


    /**
     * Compute all required probability distributions for calculating the overall structured coalescent density.
     * This methods computes equations 11 and 12 from the BASTA manuscript and populates the following arrays
     * in order to enable a straight forward structured coalescent density computation in calculateLogLikelihood:
     * intervalStartProbs, intervalStartSquareProbs, intervalEndProbs and intervalEndSquareProbs, and for the
     * contributions at the coalescent events: coalescentLeftProbs and coalescentRightProbs
     */
    private void computeProbabilityDistributions(int startingInterval) {

        //System.out.println("\n>computeProbabilityDistributions");

        //TODO still clear the activeNodeNumbers if startingInterval != 0 ?
        this.activeNodeNumbers.clear();

        //first sampling event is not considered so take this into account
        NodeRef node = intervals.getSamplingNode(-1);
        //System.out.println("sampling node: "+ treeModel.getNodeTaxon(node).getId());

        int offset = node.getNumber() * demes;

        //set start probabilities of first sampling interval
        for (int k = 0; k < demes; k++) {
            this.activeLineages[offset+k] = 0.0;
        }
        this.activeLineages[offset+patternList.getPattern(0)[patternList.getTaxonIndex(treeModel.getNodeTaxon(node).getId())]] = 1.0;
        this.activeNodeNumbers.add(node.getNumber());

        //this.intervalStartProbs[0] = Arrays.copyOf(lineageCount, demes);
        //this.intervalStartSquareProbs[0] = Arrays.copyOf(lineageCount, demes);
        this.intervalStartProbs[0] = Arrays.copyOfRange(this.activeLineages, offset, offset+demes);
        this.intervalStartSquareProbs[0] = Arrays.copyOfRange(this.activeLineages, offset, offset+demes);
        //add to list of active lineages
        //this.activeLineages.add(lineageCount);

        //print content of active lineages
        //printActiveLineages();

        //each interval is defined by its ending event, i.e. does the interval end in a sampling or coalescent event
        for (int i = 0; i < intervalCount; i++) {

            //System.out.println("interval type: " + intervals.getIntervalType(i));

            if (intervals.getIntervalType(i) == IntervalType.COALESCENT) {

                handleCoalescense(node, offset, i);

            } else if (intervals.getIntervalType(i) == IntervalType.SAMPLE) {

                handleSampling(node, offset, i);

            }
        }

    }

    private void handleSampling(NodeRef node, int offset, int i) {
        //check for zero-length interval
        if (intervals.getInterval(i) == 0.0) {
            //multiple samples at same sampling time
            //System.out.println("zero-length interval");

            node = intervals.getSamplingNode(i);
            //System.out.println("sampling node: "+ treeModel.getNodeTaxon(node).getId());
            for (int k = 0; k < demes; k++) {
                this.activeLineages[node.getNumber()*demes+k] = 0.0;
            }
            this.activeLineages[node.getNumber()*demes+patternList.getPattern(0)[patternList.getTaxonIndex(treeModel.getNodeTaxon(node).getId())]] = 1.0;
            this.activeNodeNumbers.add(node.getNumber());

            //TODO initiate caching for samples with identical sampling time and location here?

            //printActiveLineages();

        } else {
            //one sample at a time

            double intervalLength = intervals.getInterval(i);

            //printActiveLineages();

            //first compute Sd for the first interval half
//            for (int k = 0; k < demes; k++) {
//                this.intervalStartProbs[i][k] = 0.0;
//                this.intervalStartSquareProbs[i][k] = 0.0;
//                for (int l : activeNodeNumbers) {
//                    this.intervalStartProbs[i][k] += activeLineages[l*demes+k];
//                    this.intervalStartSquareProbs[i][k] += activeLineages[l*demes+k]*activeLineages[l*demes+k];
//                }
//            }
            newHardWork(intervalStartProbs[i], intervalStartSquareProbs[i]);

            //matrix exponentiation to compute end interval probabilities; equation 11
            incrementActiveLineages(this.activeLineages, intervalLength, i);

            //compute Sd for the second interval half
//            for (int k = 0; k < demes; k++) {
//                this.intervalEndProbs[i][k] = 0.0;
//                this.intervalEndSquareProbs[i][k] = 0.0;
//                for (int l : activeNodeNumbers) {
//                    this.intervalEndProbs[i][k] += activeLineages[l*demes+k];
//                    this.intervalEndSquareProbs[i][k] += activeLineages[l*demes+k]*activeLineages[l*demes+k];
//                }
//            }
            newHardWork(intervalEndProbs[i], intervalEndSquareProbs[i]);

            //printIntervalContributions(i);

            //get the node number of the sampling node
            node = intervals.getSamplingNode(i);
            //System.out.println("sampling node: "+ treeModel.getNodeTaxon(node).getId());

            //set start probabilities of first sampling interval
            for (int k = 0; k < demes; k++) {
                this.activeLineages[node.getNumber()*demes+k] = 0.0;
            }
            this.activeLineages[node.getNumber()*demes+patternList.getPattern(0)[patternList.getTaxonIndex(treeModel.getNodeTaxon(node).getId())]] = 1.0;
            this.activeNodeNumbers.add(node.getNumber());

            //printActiveLineages();

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

        //System.out.println("incrementActiveLineages with intervalLength: " + intervalLength);

        final double branchRate;
        synchronized (branchRateModel) {
            branchRate = branchRateModel.getBranchRate(treeModel, treeModel.getRoot());
        }

        if (!matricesKnown[interval]) {
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

            matricesKnown[interval] = true;
        }

//        double[] temp = new double[demes];

        //compute all dot products / probability densities (i.e. equation 11) for all active nodes
        for (int active : activeNodeNumbers) {
//            double[] temp = new double[demes];
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
     * @param intervalLength length of the coalescent interval
     * @param interval denotes the index of the interval to get the correct matrix exponential
     */
    //TODO remove this method
    private void incrementActiveLineages(ArrayList<double[]> lineageCount, double intervalLength, int interval) {

        //System.out.println("incrementActiveLineages with intervalLength: " + intervalLength);

        final double branchRate;
        synchronized (branchRateModel) {
            branchRate = branchRateModel.getBranchRate(treeModel, treeModel.getRoot());
        }

        if (!matricesKnown[interval]) {
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

            matricesKnown[interval] = true;
        }

        for (double[] linCount : lineageCount) {
            //temporary array
//            double[] temp = new double[linCount.length];
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
                System.out.println("TreeChangedEvent");
                if (((TreeChangedEvent) object).isNodeChanged()) {
                    System.out.println("current tree = " + treeModel);
                    System.out.println("isNodeChanged: " + ((TreeChangedEvent) object).getNode().getNumber());
                    System.out.println("root node number: " + treeModel.getRoot().getNumber());
                    //System.out.println(treeModel.getNodeHeight(((TreeChangedEvent) object).getNode());
                    //double changeHeight = treeModel.getNodeHeight(((TreeChangedEvent) object).getNode());
                } else if (((TreeChangedEvent) object).isHeightChanged()) {
                    System.out.println("isHeightChanged: " + ((TreeChangedEvent) object).getNode().getNumber());
                } else if (((TreeChangedEvent) object).isTreeChanged()) {
                    System.out.println("isTreeChanged");
                    // Full tree events result in a complete updating of the tree likelihood
                    // This event type is now used for EmpiricalTreeDistributions.
                    System.err.println("Full tree update event - these events currently aren't used\n" +
                            "so either this is in error or a new feature is using them so remove this message.");
                } else {
                    //System.out.println("else: ?");
                    //Other event types are ignored (probably trait changes).
                    System.err.println("Another tree event has occurred (possibly a trait change).");
                }
            } else if (object instanceof IntervalChangedEvent) {
                //these are the kinds of events we expect from BestSignalsFromBigFastTreeIntervals
                //System.out.println("IntervalChangedEvent");

                //for all the nodes that are older than the event, set needsUpdate to true
                //then trigger a recalculation that makes use of an adjusted traverseTree method (that checks whether
                //or not the ProbDist needs to be updated

            }

        } else if (model == branchRateModel) {
            for (int i = 0; i < intervalCount; i++) {
                matricesKnown[i] = false;
            }
            likelihoodKnown = false;
            areStatesRedrawn = false;
        } else if (model == generalSubstitutionModel) {
            for (int i = 0; i < intervalCount; i++) {
                matricesKnown[i] = false;
            }
            likelihoodKnown = false;
            areStatesRedrawn = false;
        } else if (model == intervals) {
            //TODO use an array of matricesKnown and set the ones between minimum and maxium changed height to false
            for (int i = 0; i < intervalCount; i++) {
                matricesKnown[i] = false;
            }
            likelihoodKnown = false;
            areStatesRedrawn = false;
        } else {
            throw new RuntimeException("Unknown handleModelChangedEvent source, exiting.");
        }

        fireModelChanged();

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
        for (int i = 0; i < intervalCount; i++) {
            matricesKnown[i] = true;
        }
    }

    protected void storeState() {
        //super.storeState();
        for (int i = 0; i < intervalCount; i++) {
            System.arraycopy(this.migrationMatrices[i],0, this.storedMigrationMatrices[i], 0, demes*demes);
            this.storedMatricesKnown[i] = this.matricesKnown[i];
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
        for (int i = 0; i < intervalCount; i++) {
            tmp = this.migrationMatrices[i];
            this.migrationMatrices[i] = this.storedMigrationMatrices[i];
            this.storedMigrationMatrices[i] = tmp;
            this.matricesKnown[i] = this.storedMatricesKnown[i];
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
        for (int i = 0; i < intervalCount; i++) {
            matricesKnown[i] = false;
        }
        likelihoodKnown = false;
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
        else if (incx == 1 && incy == 1){
            for (int c = 0, xi = dxIdx, yi = dyIdx; c < n; ++c, ++xi, ++yi) {
                s += dx[xi] * dy[yi];
            }
        } else {
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

    //private BigFastTreeIntervals intervals;
    private BestSignalsFromBigFastTreeIntervals intervals;

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
    //TODO one array as currently used isn't going to cut it to store all the probability densities
    private double[] activeLineages;
    //also keep track of the active node numbers per interval
    private ArrayList<Integer> activeNodeNumbers;

    //the migration model
    private GeneralSubstitutionModel generalSubstitutionModel;

    //number of demes for the structured coalescent model
    private int demes;

    //variables that allow to use storeState and restoreState
    private final int intervalCount;
    private double[][] migrationMatrices;
    private double[][] storedMigrationMatrices;

    private boolean[] matricesKnown;
    private boolean[] storedMatricesKnown;

}

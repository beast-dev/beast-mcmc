/*
 * StructuredColourSampler.java
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

package dr.evolution.colouring;

import dr.evolution.alignment.Alignment;
import dr.evolution.coalescent.structure.MetaPopulation;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.math.MathUtils;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.TreeMap;

/**
 * @author Gerton Lunter
 * @version $Id: StructuredColourSampler.java,v 1.11 2006/09/11 09:33:01 gerton Exp $
 *          <p/>
 *          Samples tree colourings from a proposal distribution based on
 *          Felsenstein mutation process; biased at the nodes for coalescent rates;
 *          and biased at edges using a mean-field approximation of the branch
 *          colourings.  This should improve performance for unequal population sizes
 */
public class StructuredColourSampler implements ColourSampler {

    static final int maxIterations = 1000;
    static final double tinyTime = 1.0e-6;          // to deal with discontinuous demograhpic models

    // debugging flags
    static final boolean debugMessages = false;
    static final boolean debugMeanColours = false;
    static final boolean debugNodePartials = false;
    static final boolean debugSampleLikelihoods = false;
    static final boolean debugRejectionSampler = false;
    static final boolean debugProposalProbabilityCalculator = false;
    double _totalIntegratedRate;
    static final DecimalFormat df = new DecimalFormat("###.####");

    // tuning parameters
    static final double propAffected = 0.0;        /* proportion of lineages deemed affected by conditioning on a single one */
    private boolean useNodeBias = false;
    private boolean useBranchBias = false;
    private boolean useSecondColourIteration = true;

    public StructuredColourSampler(Alignment tipColours, Tree tree, boolean nodeBias, boolean branchBias, boolean secondIteration) {

        if (tipColours.getSiteCount() != 1) {
            throw new IllegalArgumentException("Tip colour alignment must consist of a single column!");
        }

        nodeColours = new int[tree.getNodeCount()];
        colourCount = tipColours.getDataType().getStateCount();

        leafColourCounts = new int[colourCount];
        // initialize external node colours
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            NodeRef node = tree.getExternalNode(i);
            int colour = tipColours.getState(tipColours.getTaxonIndex(tree.getTaxonId(i)), 0);
            nodeColours[node.getNumber()] = colour;
            leafColourCounts[colour]++;
        }

        useNodeBias = nodeBias;
        useBranchBias = branchBias;
        useSecondColourIteration = secondIteration;

        initialize(tree);     // uses only tree constants (number of tips, etc)

    }

    public StructuredColourSampler(TaxonList[] tipColours, Tree tree, boolean nodeBias, boolean branchBias, boolean secondIteration) {

        nodeColours = new int[tree.getNodeCount()];
        colourCount = tipColours.length + 1;

        leafColourCounts = new int[colourCount];

        // initialize external node colours
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            NodeRef node = tree.getExternalNode(i);
            int colour = 0;
            for (int j = 0; j < tipColours.length; j++) {
                if (tipColours[j].getTaxonIndex(tree.getTaxonId(i)) != -1) {
                    colour = j + 1;
                }
            }
            nodeColours[node.getNumber()] = colour;
            leafColourCounts[colour]++;
        }

        useNodeBias = nodeBias;
        useBranchBias = branchBias;
        useSecondColourIteration = secondIteration;

        initialize(tree);     // uses only tree constants (number of tips etc)

    }

    public int[] getLeafColourCounts() {
        return leafColourCounts;
    }

    private void initialize(Tree tree) {

        nodePartials = new double[tree.getNodeCount()][colourCount];
        meanColourCounts = new double[tree.getNodeCount()][colourCount];
        nodeColoursEM = new int[tree.getNodeCount()][];
        nodePartialsEM = new double[tree.getNodeCount()][][];
        equilibriumColours = new double[colourCount];

    }

    private void computeIntervals(Tree tree, MetaPopulation mp) {

        // dumb implementation using a sorted map - slowish but simple
        TreeMap<Double, ArrayList<NodeRef>> intervals = new TreeMap<Double, ArrayList<NodeRef>>();
        int numnodes = tree.getNodeCount();
        for (int i = 0; i < numnodes; i++) {
            NodeRef node = tree.getNode(i);
            Double height = new Double(tree.getNodeHeight(node));
            if (intervals.containsKey(height)) {
                (intervals.get(height)).add(node);
            } else {
                ArrayList<NodeRef> list = new ArrayList<NodeRef>(1);
                list.add(node);
                intervals.put(height, list);
            }
        }
        node2Interval = new int[numnodes];
        interval2Height = new double[intervals.size()];
        avgN0 = new double[intervals.size()];
        avgN1 = new double[intervals.size()];
        Iterator<Double> iter = intervals.keySet().iterator();
        int interval = 0;
        while (iter.hasNext()) {
            Double height = iter.next();
            interval2Height[interval] = height.doubleValue();
            List<NodeRef> nodes = intervals.get(height);
            for (int i = 0; i < nodes.size(); i++) {
                node2Interval[(nodes.get(i)).getNumber()] = interval;
            }
            // now initialize the effective (i.e. harmonic average) N's
            if (interval > 0) {
                double prevtime = interval2Height[interval - 1];
                double curtime = height.doubleValue();
                avgN0[interval - 1] = (curtime - prevtime) / mp.getIntegral(prevtime, curtime, 0);
                avgN1[interval - 1] = (curtime - prevtime) / mp.getIntegral(prevtime, curtime, 1);
            }
            interval += 1;
        }
        numIntervals = interval;
    }


    /**
     * Main entry point. Colours the tree probabilistically with the given migration rates
     *
     * @param colourChangeMatrix the colour change rate parameters
     */
//    public TreeColouring sampleTreeColouring(Tree tree, ColourChangeMatrix colourChangeMatrix, double[] N) {
    public DefaultTreeColouring sampleTreeColouring(Tree tree, ColourChangeMatrix colourChangeMatrix, MetaPopulation mp) {

        //double[] N = mp.getPopulationSizes(0);

        populateEquilibriumColourArray(colourChangeMatrix);

        // Build array of node indices, arranged by height
        computeIntervals(tree, mp);

        DefaultTreeColouring colouring = new DefaultTreeColouring(2, tree);

        // Calculate root partials, and those of all other nodes
        logNodePartialsRescaling = 0.0;
        prune(tree, tree.getRoot(), colourChangeMatrix);

        // fill meanColourCounts array
        calculateMeanColourCounts(tree, colourChangeMatrix);

        // Prune again, but now using mean colour counts for each interval
        logNodePartialsRescaling = 0.0;
        double[] rootPartials = pruneEM(tree, tree.getRoot(), colourChangeMatrix, mp);

        // Re-calculate meanColourCounts, and do a final pruning step
        if (useSecondColourIteration) {
            calculateMeanColourCountsEM(tree, tree.getRoot(), colourChangeMatrix);
            logNodePartialsRescaling = 0.0;
            rootPartials = pruneEM(tree, tree.getRoot(), colourChangeMatrix, mp);
        }

        // Sampling is conditional on data; so normalize by the probability of the
        // data under the proposal distribution
        double normalization = 0.0;
        for (int i = 0; i < colourCount; i++) {
            normalization += equilibriumColours[i] * rootPartials[i];
        }
        double logNormalization = Math.log(normalization) + logNodePartialsRescaling;

        double logP = sampleEM(tree, tree.getRoot(), colourChangeMatrix, mp, colouring) - logNormalization;

        colouring.setLogProbabilityDensity(logP);

        //JFrame frame = new JFrame();
        //ColouredTreeComponent comp = new ColouredTreeComponent()

        if (debugProposalProbabilityCalculator) {
            double logP2 = getProposalProbability(colouring, tree, colourChangeMatrix, mp);
            System.out.println("Sampling proposal probability " + logP + ", re-calculated as " + logP2 + " (norm=" + logNormalization + ")");
        }

        return colouring;
    }

    /**
     * Returns proposal probability density of a tree colouring
     *
     * @param treeColouring
     * @param tree
     * @param colourChangeMatrix
     * @param mp
     * @return probability density
     */
    public double getProposalProbability(TreeColouring treeColouring, Tree tree, ColourChangeMatrix colourChangeMatrix, MetaPopulation mp) {

        //double[] N = mp.getPopulationSizes(0);

        populateEquilibriumColourArray(colourChangeMatrix);

        // Build array of node indices, arranged by height
        computeIntervals(tree, mp);

        // Calculate root partials, and those of all other nodes
        logNodePartialsRescaling = 0.0;
        prune(tree, tree.getRoot(), colourChangeMatrix);

        // fill meanColourCounts array
        calculateMeanColourCounts(tree, colourChangeMatrix);

        // Prune again, but now using mean colour counts for each interval
        logNodePartialsRescaling = 0.0;
        double[] rootPartials = pruneEM(tree, tree.getRoot(), colourChangeMatrix, mp);

        // Re-calculate meanColourCounts, and do a final pruning step
        if (useSecondColourIteration) {
            calculateMeanColourCountsEM(tree, tree.getRoot(), colourChangeMatrix);
            logNodePartialsRescaling = 0.0;
            rootPartials = pruneEM(tree, tree.getRoot(), colourChangeMatrix, mp);
        }

        // Sampling is conditional on data; so normalize by the probability of the
        // data under the proposal distribution
        double normalization = 0.0;
        for (int i = 0; i < colourCount; i++) {
            normalization += equilibriumColours[i] * rootPartials[i];
        }

        double logP = calculateEMProposal(tree, tree.getRoot(), colourChangeMatrix, mp, treeColouring);

        return logP - Math.log(normalization) - logNodePartialsRescaling;

    }


    /**
     * @param node
     */
    private int getColour(NodeRef node) {
        return nodeColours[node.getNumber()];
    }

    /**
     * @param node
     */
    private void setColour(NodeRef node, int colour) {
        if (colour >= 0 && colour < colourCount) {
            nodeColours[node.getNumber()] = colour;
        } else {
            throw new IllegalArgumentException("colour value " + colour + " + is outside of range of colours, [0, " + Integer.toString(colourCount - 1) + "]");
        }
    }


    /**
     * **************************************************************************************
     * <p/>
     * Probability- and sampling-related code follows
     */


    void populateEquilibriumColourArray(ColourChangeMatrix colourChangeMatrix) {

        // fill equilibrium colour array
        for (int clr = 0; clr < colourCount; clr++) {
            equilibriumColours[clr] = colourChangeMatrix.getEquilibrium(clr);
        }

    }


    /* Returns probability of colours at a node, based on pruning output */
    double[] getMeanColours(int nodeNum, ColourChangeMatrix colourChangeMatrix) {

        double c[] = new double[colourCount];
        double norm = 0.0;
        for (int i = 0; i < colourCount; i++) {
            c[i] = nodePartials[nodeNum][i] * equilibriumColours[i];
            norm += c[i];
        }
        for (int i = 0; i < colourCount; i++) {
            c[i] /= norm;
        }
        return c;
    }

    /* Returns probability of colours at an node-based intreval, based on pruning output */
    double[] getMeanColoursEM(int nodeNum, int relativeInterval, ColourChangeMatrix colourChangeMatrix) {

        double c[] = new double[colourCount];
        double norm = 0.0;
        for (int i = 0; i < colourCount; i++) {
            c[i] = nodePartialsEM[nodeNum][relativeInterval][i] * equilibriumColours[i];
            norm += c[i];
        }
        for (int i = 0; i < colourCount; i++) {
            c[i] /= norm;
        }
        return c;
    }


    // Recursively fills meanColourCounts with derived counts for intervals;
    // utility function for computeMeanColourCounts below
    void fillMeanColourCounts(Tree tree, NodeRef node, ColourChangeMatrix colourChangeMatrix) {

        if (!tree.isRoot(node)) {

            NodeRef parent = tree.getParent(node);
            int parentNum = parent.getNumber();
            int nodeNum = node.getNumber();

            double[] parentColours = getMeanColours(parentNum, colourChangeMatrix);
            double[] childColours = getMeanColours(nodeNum, colourChangeMatrix);

            for (int i = 0; i < colourCount; i++) {
                // calculate expected colour at parent and child end of branch
                double meanColour = (parentColours[i] + childColours[i]) / 2.0;
                meanColourCounts[node2Interval[nodeNum]][i] += meanColour;
                meanColourCounts[node2Interval[parentNum]][i] -= meanColour;
            }

        }

        if (!tree.isExternal(node)) {

            // Note: assuming binary tree!
            NodeRef leftChild = tree.getChild(node, 0);
            NodeRef rightChild = tree.getChild(node, 1);

            fillMeanColourCounts(tree, leftChild, colourChangeMatrix);
            fillMeanColourCounts(tree, rightChild, colourChangeMatrix);

        }

    }


    // Recursively fills meanColourCounts with derived counts for intervals;
    // utility function for computeMeanColourCountsEM below
    void calculateMeanColourCountsEM(Tree tree, NodeRef node, ColourChangeMatrix colourChangeMatrix) {

        if (tree.isRoot(node)) {

            // empty the meanColourCounts array
            for (int clr = 0; clr < colourCount; clr++) {
                for (int i = 0; i < meanColourCounts.length; i++) {
                    meanColourCounts[i][clr] = 0.0;
                }
            }

        } else {

            NodeRef parent = tree.getParent(node);
            int parentNum = parent.getNumber();
            int nodeNum = node.getNumber();

            int bottomInterval = node2Interval[nodeNum];
            int topInterval = node2Interval[parentNum];

            double[] belowColours = getMeanColoursEM(nodeNum, 0, colourChangeMatrix);
            double[] aboveColours;

            for (int relInterval = 0; relInterval < topInterval - bottomInterval; relInterval++) {

                int aboveInterval = relInterval + 1;

                if (aboveInterval + bottomInterval < topInterval) {
                    aboveColours = getMeanColoursEM(nodeNum, aboveInterval, colourChangeMatrix);
                } else {
                    if (!tree.isRoot(parent)) {
                        aboveColours = getMeanColoursEM(parentNum, 0, colourChangeMatrix);
                    } else {
                        aboveColours = belowColours;
                    }
                }

                for (int i = 0; i < colourCount; i++) {
                    double meanColour = (aboveColours[i] + belowColours[i]) / 2.0;
                    meanColourCounts[relInterval + bottomInterval][i] += meanColour;
                }

                belowColours = aboveColours;
            }

        }

        if (!tree.isExternal(node)) {

            // Note: assuming binary tree!
            NodeRef leftChild = tree.getChild(node, 0);
            NodeRef rightChild = tree.getChild(node, 1);

            calculateMeanColourCountsEM(tree, leftChild, colourChangeMatrix);
            calculateMeanColourCountsEM(tree, rightChild, colourChangeMatrix);

        }

        if (tree.isRoot(node) && debugMeanColours) {
            for (int i = 0; i < numIntervals; i++) {
                System.out.println("EM Interval " + i + " height=" + interval2Height[i] + " ");
                double sum = 0.0;
                for (int clr = 0; clr < colourCount; clr++) {
                    System.out.print("c=" + clr + " mean=" + meanColourCounts[i][clr] + " ");
                    sum += meanColourCounts[i][clr];
                }
                System.out.println(" sum=" + sum);
            }
        }
    }


    // Computes the meanColourCounts array for the intervals
    void calculateMeanColourCounts(Tree tree, ColourChangeMatrix colourChangeMatrix) {

        // empty the meanColourCounts array
        for (int clr = 0; clr < colourCount; clr++) {
            for (int i = 0; i < meanColourCounts.length; i++) {
                meanColourCounts[i][clr] = 0.0;
            }
        }

        // populate meanColourCounts array, so that sum of counts before and including an interval,
        // give average colour counts for that interval.
        fillMeanColourCounts(tree, tree.getRoot(), colourChangeMatrix);

        // Sum out, so that actual count gives average colour count on interval
        for (int clr = 0; clr < colourCount; clr++) {
            double sum = 0.0;
            for (int i = 0; i < meanColourCounts.length; i++) {
                sum += meanColourCounts[i][clr];
                meanColourCounts[i][clr] = sum;
            }
        }

        if (debugMeanColours) {
            for (int i = 0; i < numIntervals; i++) {
                System.out.println("Interval " + i + " height=" + interval2Height[i] + " ");
                double sum = 0.0;
                for (int clr = 0; clr < colourCount; clr++) {
                    System.out.print("c=" + clr + " mean=" + meanColourCounts[i][clr] + " ");
                    sum += meanColourCounts[i][clr];
                }
                System.out.println(" sum=" + sum);
            }
        }
    }


    /**
     * Calculate probability of data at descendants from node, given a color at the node ('partials'),
     * by a Felsenstein-like pruning algorithm.  (First step in the color sampling algorithm)
     * Side effect: updates nodePartials[] for this node and all its descendants.
     *
     * @param node
     * @return the partials of this node
     */
    private double[] prune(Tree tree, NodeRef node, ColourChangeMatrix mm) {

        double[] p = new double[colourCount];

        if (tree.isExternal(node)) {
            p[getColour(node)] = 1.0;
        } else {

            // Note: assuming binary tree!
            NodeRef leftChild = tree.getChild(node, 0);
            NodeRef rightChild = tree.getChild(node, 1);

            double[] left = prune(tree, leftChild, mm);
            double[] right = prune(tree, rightChild, mm);

            double nodeHeight = tree.getNodeHeight(node);

            double leftTime = nodeHeight - tree.getNodeHeight(tree.getChild(node, 0));
            double rightTime = nodeHeight - tree.getNodeHeight(tree.getChild(node, 1));

            double maxp = 0.0;
            for (int i = 0; i < colourCount; i++) {

                double leftSum = 0.0;
                double rightSum = 0.0;

                // looping over colours in left and right children
                for (int j = 0; j < colourCount; j++) {
                    // forwardTimeEvolution conditions on the parent state i, i.e. time
                    // runs in the natural direction (forward from parent to child)
                    leftSum += mm.forwardTimeEvolution(i, j, leftTime) * left[j];
                    rightSum += mm.forwardTimeEvolution(i, j, rightTime) * right[j];
                }
                p[i] = leftSum * rightSum;
                // This version does not condition on the formal variable, because
                // that tends to bias towards the small population, whereas no correction
                // for bias towards the large population on branches is made at this stage

                if (p[i] > maxp) {
                    maxp = p[i];
                }

            }

            // take care of underflow
            if (maxp < 1.0e-100) {
                for (int i = 0; i < colourCount; i++) {
                    p[i] *= 1.0e+100;
                }
                logNodePartialsRescaling -= Math.log(1.0e+100);
            }

        }

        nodePartials[node.getNumber()] = p;

        if (debugNodePartials) {
            prettyPrint("Node " + node.getNumber() + " prune=", p);
        }

        return p;
    }


    /* Calculates the result of multiplying the vector vec with the matrix Exp[{{-a,b},{c,-d}}] */
    static double[] matrixEvolve(double[] mx, int current) {

        double a = mx[0];
        double b = mx[1];
        double c = mx[2];
        double d = mx[3];

        /*
          double[] vec = new double[2];
          vec[current] = 1.0;
          double S = Math.sqrt( (a-d)*(a-d) + 4*b*c );
          double c00,c01,c10,c11;
          if (S < 1.0e-5 ) {
              c00 = Math.exp( -a );
              c11 = Math.exp( -d );
              c01 = b*c00;
              c10 = c*c11;
          } else {
              double T = Math.exp( -(a+d+S)/2.0 );
              double U = Math.exp( -(a+d-S)/2.0 );
              c00 = ( (d-a+S)*U - (d-a-S)*T ) / (2*S);
              c01 = ( b*(U-T) ) / S;
              c10 = ( c*(U-T) ) / S;
              c11 = ( (a-d+S)*U - (a-d-S)*T ) / (2*S);
          }
          double v0 = vec[0] * c00 + vec[1] * c10;
          vec[1] = vec[0] * c01 + vec[1] * c11;
          vec[0] = v0;
          return vec;
          */

        double S = Math.sqrt((a - d) * (a - d) + 4 * b * c);
        if (S < 1.0e-5) {
            if (current == 0) {
                a = Math.exp(-a);
                return new double[]{a, b * a};
            } else {
                d = Math.exp(-d);
                return new double[]{c * d, d};
            }
        } else {
            double T = Math.exp(-(a + d + S) / 2.0);
            double U = Math.exp(-(a + d - S) / 2.0);
            if (current == 0) {
                return new double[]{((d - a + S) * U - (d - a - S) * T) / (2 * S), (b * (U - T)) / S};
            } else {
                return new double[]{(c * (U - T)) / S, ((a - d + S) * U - (a - d - S) * T) / (2 * S)};
            }
        }
    }


    /* Calculates the result of multiplying the matrix Exp[{{-a,b},{c,-d}}] with vec */
    static void matrixPullBack(double[] mx, double[] vec) {

        double a = mx[0];
        double b = mx[1];
        double c = mx[2];
        double d = mx[3];

        double S = Math.sqrt((a - d) * (a - d) + 4 * b * c);
        double c00, c01, c10, c11;
        if (S < 1.0e-5) {
            c00 = Math.exp(-a);
            c11 = Math.exp(-d);
            c01 = b * c00;
            c10 = c * c11;
        } else {
            double T = Math.exp(-(a + d + S) / 2.0);
            double U = Math.exp(-(a + d - S) / 2.0);
            c00 = ((d - a + S) * U - (d - a - S) * T) / (2 * S);
            c01 = (b * (U - T)) / S;
            c10 = (c * (U - T)) / S;
            c11 = ((a - d + S) * U - (a - d - S) * T) / (2 * S);
        }
        double v0 = vec[0] * c00 + vec[1] * c01;
        vec[1] = vec[0] * c10 + vec[1] * c11;
        vec[0] = v0;
    }


    double[] calculateMatrixElts(int curinterval, NodeRef affectedNode, Tree tree, double time, double N0, double N1, ColourChangeMatrix mm) {

        double meancolour0 = meanColourCounts[curinterval][0];
        double meancolour1 = meanColourCounts[curinterval][1];
        double lineages = meancolour0 + meancolour1;

        // get coalescent rate for colour 0
        double coalescent0 = ((propAffected * (lineages - 1)) + (1.0 - propAffected) * (meancolour0 - 1)) / (2.0 * N0) * time;
        if (coalescent0 < 0.0) {
            coalescent0 = 0.0;
        }
        // same for colour 1
        double coalescent1 = ((propAffected * (lineages - 1)) + (1.0 - propAffected) * (meancolour1 - 1)) / (2.0 * N1) * time;
        if (coalescent1 < 0.0) {
            coalescent1 = 0.0;
        }
        if (!useBranchBias) {
            coalescent0 = 0.0;
            coalescent1 = 0.0;
        }
        // partial conditioning on survival; take away coalescences that are
        // common to both colours
        double commonCoalescentRate = Math.min(coalescent0, coalescent1);
        /* debug */
        //commonCoalescentRate = 0.0;
        coalescent0 -= commonCoalescentRate;
        coalescent1 -= commonCoalescentRate;
        // get forward rates
        double fwrate01 = mm.getForwardRate(0, 1) * time;
        double fwrate10 = mm.getForwardRate(1, 0) * time;
        // next, evolve 'with leakage'.  Matrix:
        // 0->1 : rate01;   diagonal 0 : -rate01-coalescent0
        // 1->0 : rate10;   diagonal 1 : -rate10-coalescent1
        return new double[]{fwrate01 + coalescent0, fwrate01, fwrate10, fwrate10 + coalescent1};
    }


    /* Calculates partial branch probabilities with mean-field coalescent correction */
    /* Inputs are the partials at the child node; results on branch are stored in nodePartialsEM[][][] */
    /* This version is for two colours only */
    double[] pruneBranchEM(ColourChangeMatrix mm, double[] inState, NodeRef parent, NodeRef child, Tree tree, MetaPopulation mp) {

        // Get interval above parent (highest) and above child (lowest)
        int topInterval = node2Interval[parent.getNumber()];
        int bottomInterval = node2Interval[child.getNumber()];

        // Make a result array and a current state
        double[][] partials = new double[topInterval - bottomInterval][2];
        double[] state = inState.clone();
        int curinterval = bottomInterval;

        while (curinterval != topInterval) {

            // Store current state
            partials[curinterval - bottomInterval][0] = state[0];
            partials[curinterval - bottomInterval][1] = state[1];
            // Calculate partials at top of current interval
            double time = interval2Height[curinterval + 1] - interval2Height[curinterval];
            //double effectiveN0 = time / mp.getIntegral( lowertime, uppertime, 0 );
            //double effectiveN1 = time / mp.getIntegral( lowertime, uppertime, 1 );
            double[] mxElts = calculateMatrixElts(curinterval, child, tree, time, avgN0[curinterval], avgN1[curinterval], mm);
            matrixPullBack(mxElts, state);
            // update for next interval
            curinterval += 1;
        }

        // Store branch partials (including child's, excluding parent's) at child
        nodePartialsEM[child.getNumber()] = partials;

        // Return parent's partials
        return state;
    }


    /**
     * Calculate probability of data at descendants from node, given a color at the node ('partials'),
     * by a Felsenstein-like pruning algorithm.  (First step in the color sampling algorithm)
     * Side effect: updates nodePartials[] and nodePartialsEM[][] for this node and all its descendants.
     * This version uses mean-field coalescence correction and formal variable correction
     */
    private double[] pruneEM(Tree tree, NodeRef node, ColourChangeMatrix mm, MetaPopulation mp) {

        double[] p = new double[colourCount];

        if (tree.isExternal(node)) {

            p[getColour(node)] = 1.0;

        } else {

            // Note: assuming binary tree!
            NodeRef leftChild = tree.getChild(node, 0);
            NodeRef rightChild = tree.getChild(node, 1);

            // Obtain partials at the child nodes
            double[] left = pruneEM(tree, leftChild, mm, mp);
            double[] right = pruneEM(tree, rightChild, mm, mp);

            // Calculate partials at parent node
            // (Side effect: calculates partials for each interval along branch in nodePartialsEM)
            double[] leftBranchPartials = pruneBranchEM(mm, left, node, leftChild, tree, mp);
            double[] rightBranchPartials = pruneBranchEM(mm, right, node, rightChild, tree, mp);

            double maxp = 0.0;
            for (int i = 0; i < colourCount; i++) {

                p[i] = leftBranchPartials[i] * rightBranchPartials[i];

                // Condition on the formal variable
                if (useNodeBias) {
                    p[i] *= mm.getEquilibrium(i) / mp.getDemographic(tree.getNodeHeight(node) - tinyTime, i);
                }

                if (p[i] > maxp) {
                    maxp = p[i];
                }

            }

            // take care of underflow
            if (maxp < 1.0e-100) {
                for (int i = 0; i < colourCount; i++) {
                    p[i] *= 1.0e+100;
                }
                logNodePartialsRescaling -= Math.log(1.0e+100);
            }

        }

        // Store results for the sampling step
        nodePartials[node.getNumber()] = p;

        if (debugNodePartials) {
            prettyPrint("Node " + node.getNumber() + " pruneEM=", p);
        }

        return p;
    }


    /**
     * Samples internal node colours (from root to tips) and events
     * Precondition: parent node has been sampled, descendants have not
     * Requires the results from Felsenstein Backwards pruning, in nodePartials(EM)[] (see pruneEM())
     * Side effect: updates nodeColours[]
     */
    private double sampleEM(Tree tree, NodeRef node, ColourChangeMatrix mm, MetaPopulation mp, DefaultTreeColouring colouring) {

        double[] forward;
        double[] posterior;
        int colour;
        DefaultBranchColouring history = null;
        double logLikelihood = 0.0;

        if (tree.isRoot(node)) {

            _totalIntegratedRate = 0.0; /* for debugging */

            forward = mm.getEquilibrium();
            double[] backward = nodePartials[node.getNumber()];
            posterior = new double[colourCount];
            double max = -1.0;
            double min = 1.0;
            for (int i = 0; i < forward.length; i++) {
                posterior[i] = forward[i] * backward[i];
                max = Math.max(max, posterior[i]);
                min = Math.min(min, posterior[i]);
            }
            if (debugMessages && min < 0.0) {
                System.out.println("Aargh, negative probabilities " + min + " " + max);
            }
            if (debugMessages && max < 1.0e-200) {
                System.out.println("Hmm, very small numbers indeed " + max);
            }
            // Sample a colour
            colour = MathUtils.randomChoicePDF(posterior);
            // Add likelihood of (unconditional) probability
            logLikelihood += Math.log(forward[colour]);

        } else {

            int nodeNum = node.getNumber();
            double[][] backward = nodePartialsEM[nodeNum];
            int bottomInterval = node2Interval[nodeNum];
            // Allocate array for interval node colours
            nodeColoursEM[nodeNum] = new int[backward.length];
            // Start from parent, and sample colours for interval nodes
            colour = getColour(tree.getParent(node));
            // Allocate event history
            history = new DefaultBranchColouring(colour, colour);
            // Allocate posterior probabilities
            posterior = new double[colourCount];

            // Go down the branch (from parent to child)
            for (int relinterval = backward.length - 1; relinterval >= 0; relinterval--) {

                int interval = relinterval + bottomInterval;
                double childHeight = interval2Height[interval];
                double uppertime = interval2Height[interval + 1];
                double time = uppertime - childHeight;
                double[] mxElts = calculateMatrixElts(interval, node, tree, time, avgN0[interval], avgN1[interval], mm);
                forward = matrixEvolve(mxElts, colour);
                // Calculate posterior for interval node
                for (int i = 0; i < colourCount; i++) {
                    posterior[i] = forward[i] * backward[relinterval][i];
                }
                // Sample a colour
                int childColour = MathUtils.randomChoicePDF(posterior);
                // Store colour
                nodeColoursEM[nodeNum][relinterval] = childColour;
                // Sample events conditional on parent and child colour
                //System.out.println("node="+node.getNumber()+" parent="+colour+" child="+childColour+" interval="+relinterval+" time="+time);
                logLikelihood += sampleConditionalBranchColouringEM(node, colour, childColour, time, childHeight, mxElts, history);

                // update colour
                colour = childColour;
            }

            // Store history for branch at child node
            colouring.setBranchColouring(node, history);
        }

        // store colour in node
        setColour(node, colour);

        // add contribution of formal variable - except if this is a leaf,
        // since there is no corresponding coalescent event
        if (!tree.isExternal(node) && useNodeBias) {
            double childTime = tree.getNodeHeight(node);
            logLikelihood += Math.log(mm.getEquilibrium(colour) / mp.getDemographic(childTime - tinyTime, colour));
        }

        // Recursively sample down the tree
        for (int i = 0; i < tree.getChildCount(node); i++) {
            NodeRef child = tree.getChild(node, i);
            logLikelihood += sampleEM(tree, child, mm, mp, colouring);
        }

        if (debugSampleLikelihoods && tree.isRoot(node)) {
            System.out.println("Sampled likelihood " + logLikelihood);
            System.out.println("Integrated exit rate " + _totalIntegratedRate);
            System.out.println("sampleEMProposal node=" + node.getNumber() + " logL=" + logLikelihood);
        }

        //System.out.println("sampleEM node="+node.getNumber()+"\tlogL="+logLikelihood);

        return logLikelihood;
    }

    /**
     * This samples events along an interval, and adds these events to a BranchColouring.
     * Also returns the (unconditional) log probability density for these events
     */
    private double sampleConditionalBranchColouringEM(
            NodeRef node,
            int parentColour,
            int childColour,
            double time,
            double childHeight,
            double[] matrixElements,
            DefaultBranchColouring initialBranchColouring) {

        DefaultBranchColouring history = new DefaultBranchColouring(parentColour, childColour);

        int iterations = 0;
        int currentColour;
        double currentHeight, dt;
        double logLikelihood;
        boolean reject, firstEvent;
        double _sumRate = 0.0;
        String _msg = "";

        // Reject until we get the child colour
        do {

            history.clear();
            currentColour = parentColour;
            currentHeight = time;
            logLikelihood = 0.0;
            reject = false;
            firstEvent = true;
            if (debugRejectionSampler) {
                _sumRate = 0.0;
                _msg = "sample (iter=" + iterations + ") from " + parentColour + " to " + childColour + " at " + df.format(childHeight) + " over " + df.format(time) + " rate=" + df.format(matrixElements[3 * parentColour] / time) + "\n";
            }

            // Sample events until we reach the child
            do {

                // Sample a waiting time
                double totalRate, migrationRate;
                if (currentColour == 0) {
                    totalRate = matrixElements[0] / time;
                    migrationRate = matrixElements[1] / time;
                } else {
                    totalRate = matrixElements[3] / time;
                    migrationRate = matrixElements[2] / time;
                }
                double U;

                do {
                    U = MathUtils.nextDouble();
                } while (U == 0.0);

                // Neat trick (Rasmus Nielsen):
                // If colours of parent and child differ, condition on at least 1 event
                if (firstEvent && (parentColour != childColour)) {

                    double minU = Math.exp(-totalRate * time);
                    U = minU + U * (1.0 - minU);

                }

                // Calculate the waiting time, and update currentHeight
                dt = -Math.log(U) / totalRate;
                currentHeight -= dt;

                if (currentHeight > 0) {
                    // Not yet reached the child.  Sample an event.
                    // Neat trick II: condition on first event being a migration event
                    // Also, don't bother calculating a random sample if coalescent events
                    //  have zero rate
                    if (debugRejectionSampler) {
                        _sumRate += totalRate * dt;
                        _msg = _msg + "rate=" + df.format(totalRate) + "\ttime=" + df.format(dt) + "\tintrate=" + df.format(totalRate * dt) + "\n";
                    }

                    if (firstEvent ||
                            (migrationRate == totalRate) ||
                            (MathUtils.nextDouble() < migrationRate / totalRate)) {
                        // migration event
                        //System.out.println("Event height="+(currentHeight + childHeight));
                        currentColour = 1 - currentColour;
                        // Add it to the list
                        history.addEvent(currentColour, currentHeight + childHeight);

                        // Update likelihood
                        logLikelihood += (-totalRate * dt) + Math.log(migrationRate);

                        // ... replaces this to keep in log space
                        // likelihood *= Math.exp( -totalRate * dt ) * migrationRate;
                    } else {
                        // coalescent event
                        reject = true;
                    }
                } else {
                    if (debugRejectionSampler) {
                        _sumRate += totalRate * (currentHeight + dt);
                        _msg = _msg + "rate=" + df.format(totalRate) + "\ttime=" + df.format(currentHeight + dt) + "\tintrate=" + df.format(totalRate * (currentHeight + dt)) + "\n";
                    }
                    // Update likelihood, no event until child
                    //System.out.println(" totalRate="+totalRate+" dt="+(currentHeight+dt)+" time="+time);

                    logLikelihood += -totalRate * (currentHeight + dt);
                    // ... replaces this to keep in log space
                    //likelihood *= Math.exp( -totalRate * (currentHeight + dt) );
                }

                firstEvent = false;

            } while ((!reject) && (currentHeight > 0));

            iterations += 1;
            if (currentColour != childColour) {
                reject = true;
            }

        } while (reject && (iterations < maxIterations));

        if (debugRejectionSampler) {
            _totalIntegratedRate += _sumRate;
            System.out.print(_msg);
        }

        if (reject) {
            // Extreme migration rates may cause difficulty for the rejection sampler
            // Print a warning and add a bogus event somewhat near where you'd want it.
            if (currentColour != childColour) {
                double previousEventHeight = currentHeight + dt;
                double finalEventHeight = 0.01 * previousEventHeight;
                history.addEvent(childColour, finalEventHeight + childHeight);

                if (debugMessages) {
                    System.out.println("dr.evolution.colouring.BranchColourSampler: failed to generate sample after " + maxIterations + " trials.");
                    System.out.println(": height=" + time);
                    System.out.println(": childColour=" + childColour);
                    System.out.println(": childHeight=" + childHeight);
                    System.out.println(": migration rate 0->1 = " + matrixElements[1]);
                    System.out.println(": migration rate 1->0 = " + matrixElements[2]);
                    System.out.println(": total rate 0 = " + matrixElements[0]);
                    System.out.println(": total rate 1 = " + matrixElements[3]);
                }
            }
        }

        // Add this interval's history to the previously sampled history
        initialBranchColouring.addHistory(history);

        //System.out.println("sampleConditionalBranchColouringEM child="+node.getNumber()+"\tchildH="+childHeight+"\tlogL="+Math.log(likelihood)+"\tevents="+history.getNumEvents()+"\tfinal clr="+childColour);

        return logLikelihood;
    }


    /**
     * Calculates proposal density for given colouring
     * Requires the results from Felsenstein Backwards pruning, in nodePartials(EM)[] (see pruneEM())
     */
    private double calculateEMProposal(Tree tree, NodeRef node, ColourChangeMatrix mm, MetaPopulation mp, TreeColouring colouring) {

        double[] forward;
        int colour;
        BranchColouring history = null;
        double logLikelihood = 0.0;

        if (tree.isRoot(node)) {

            forward = mm.getEquilibrium();
            // Get root colour
            colour = colouring.getNodeColour(node);
            // Add likelihood of (unconditional) probability
            logLikelihood += Math.log(forward[colour]);

        } else {

            int nodeNum = node.getNumber();
            double[][] backward = nodePartialsEM[nodeNum];
            int bottomInterval = node2Interval[nodeNum];
            // Get branch colouring
            history = colouring.getBranchColouring(node);

            // Go down the branch (from parent to child)
            for (int relinterval = backward.length - 1; relinterval >= 0; relinterval--) {

                int interval = relinterval + bottomInterval;
                double childHeight = interval2Height[interval];
                double uppertime = interval2Height[interval + 1];
                double time = uppertime - childHeight;
                double[] mxElts = calculateMatrixElts(interval, node, tree, time, avgN0[interval], avgN1[interval], mm);
                logLikelihood += calculateConditionalBranchColouringEM(node, time, childHeight, mxElts, history);
            }

            colour = colouring.getNodeColour(node);

        }

        if (!tree.isExternal(node) && useNodeBias) {
            double childTime = tree.getNodeHeight(node);
            logLikelihood += Math.log(mm.getEquilibrium(colour) / mp.getDemographic(childTime - tinyTime, colour));
        }

        // Recursively sample down the tree
        for (int i = 0; i < tree.getChildCount(node); i++) {
            NodeRef child = tree.getChild(node, i);
            logLikelihood += calculateEMProposal(tree, child, mm, mp, colouring);
        }

        //System.out.println("calculateEMProposal node="+node.getNumber()+"\tlogL="+logLikelihood);

        return logLikelihood;
    }

    /**
     * This samples events along an interval, and adds these events to a BranchColouring.
     * Also returns the (unconditional) log probability density for these events
     */
    private double calculateConditionalBranchColouringEM(
            NodeRef node,
            double time,
            double childHeight,
            double[] matrixElements,
            BranchColouring branchColouring) {

        double currentHeight = time + childHeight;
        int nextEvent = branchColouring.getNextForwardEvent(currentHeight);
        int currentColour = branchColouring.getForwardColourBelow(nextEvent - 1);
        double logLikelihood = 0.0;
        int numEvents = 0;

        while (currentHeight > childHeight) {

            double nextEventHeight;
            //System.out.println("nextEvent="+nextEvent+" total="+branchColouring.getNumEvents());
            if (nextEvent == branchColouring.getNumEvents() + 1) {
                nextEventHeight = childHeight - 1.0;   // past child
            } else {
                nextEventHeight = branchColouring.getForwardTime(nextEvent);
                //System.out.println("Event height="+nextEventHeight);
            }

            double dt = currentHeight - nextEventHeight;

            double totalRate, migrationRate;
            if (currentColour == 0) {
                totalRate = matrixElements[0] / time;
                migrationRate = matrixElements[1] / time;
            } else {
                totalRate = matrixElements[3] / time;
                migrationRate = matrixElements[2] / time;
            }

            if (nextEventHeight < childHeight) {
                // no event - include exit probability, and done.
                dt = currentHeight - childHeight;
                //System.out.println(" totalRate="+totalRate+" dt="+dt+" time="+time);

                logLikelihood += -totalRate * dt;
                // ... replaces this to keep in log space
                //likelihood *= Math.exp( -totalRate * dt );
            } else {
                // proper migration event

                logLikelihood += (-totalRate * dt) + Math.log(migrationRate);
                // ... replaces this to keep in log space
                //likelihood *= Math.exp( -totalRate * dt ) * migrationRate;

                currentColour = branchColouring.getForwardColourBelow(nextEvent);
                numEvents += 1;
            }

            currentHeight = nextEventHeight;
            nextEvent += 1;

        }

        //System.out.println("calculateConditionalBranchColouring child="+node.getNumber()+"\tchildH="+childHeight+"\tlogL="+Math.log(likelihood)+"\tevents="+numEvents+"\tfinal clr="+currentColour);

        return logLikelihood;
    }


    private void prettyPrint(String label, double[] vec) {
        System.out.print(label + "= (");
        for (double x : vec) {
            System.out.print(x + ", ");

        }
        System.out.println(")");

    }

    // Testing code
    static void testMatrix(double[] m, double[] r) {
        if (Math.abs(matrixEvolve(m, 0)[0] - r[0]) > 1.0e-6) {
            throw new Error("1");
        }
        if (Math.abs(matrixEvolve(m, 0)[1] - r[1]) > 1.0e-6) {
            throw new Error("2");
        }
        if (Math.abs(matrixEvolve(m, 1)[0] - r[2]) > 1.0e-6) {
            throw new Error("3");
        }
        if (Math.abs(matrixEvolve(m, 1)[1] - r[3]) > 1.0e-6) {
            throw new Error("4");
        }
        double vec0[] = {1, 0};
        double vec1[] = {0, 1};
        matrixPullBack(m, vec0);
        matrixPullBack(m, vec1);
        if (Math.abs(vec0[0] - r[0]) > 1.0e-6) {
            throw new Error("5");
        }
        if (Math.abs(vec0[1] - r[2]) > 1.0e-6) {
            throw new Error("7");
        }
        if (Math.abs(vec1[0] - r[1]) > 1.0e-6) {
            throw new Error("6");
        }
        if (Math.abs(vec1[1] - r[3]) > 1.0e-6) {
            throw new Error("8");
        }
    }

    public static void main(String[] args) {

        double m[] = {5.0, 3.0, 2.0, 3.0};
        double result[] = {0.0811818, 0.145616, 0.097077, 0.178259};
        testMatrix(m, result);
        System.out.println("First matrix OK");
        double m2[] = {1.0, 1.0, 0.0, 1.0};
        double result2[] = {0.367879, 0.367879, 0, 0.367879};
        testMatrix(m2, result2);
        System.out.println("Second matrix OK");
        double m3[] = {1.0, 0.0, 1.0, 1.0};
        double result3[] = {0.367879, 0.0, 0.367879, 0.367879};
        testMatrix(m3, result3);
        System.out.println("Third matrix OK");

    }

    private final int colourCount;

    private final int[] nodeColours;

    private final int[] leafColourCounts;

    private double[][] meanColourCounts; // [i][c] = approximate number of lineages coloured c in interval i

    private int[][] nodeColoursEM;       // [n][r] = colour of lineage above node n and interval r (relative)

    private double[][] nodePartials;     // [node][colour] = p(descendant data | node colour)

    private double[][][] nodePartialsEM;  // [node][rel.interval][colour] = p(desc | node clr)

    private double logNodePartialsRescaling;  // rescaling factor to prevent underflow

    private double[] equilibriumColours;

    private int[] node2Interval;       // linear index for interval corresponding to (child) node

    private double[] interval2Height;

    private double[] avgN0;
    private double[] avgN1;

    private int numIntervals;
}

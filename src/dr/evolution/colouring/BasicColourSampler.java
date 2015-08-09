/*
 * BasicColourSampler.java
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


/**
 * @author Alexei Drummond
 * @author Gerton Lunter
 * @author Andrew Rambaut
 *         <p/>
 *         This is the old version.  It samples like a substitution model, except that nodes are biased
 *         towards populations with low Ne, since coalescences are more likely to occur there.
 *         <p/>
 *         It seems to work less well than a straight substitution-like model (such as Greg Ewing uses).
 *         The reason is that although the bias is correct, along branches the bias works the other way.
 *         By incorporating bias at the nodes, the problem along branches gets worse, and this seems to
 *         affect the acceptance probabilities more than having the right bias at the nodes helps.  So,
 *         although the code still accepts the population sizes, it ignores it now.
 * @version $Id: BasicColourSampler.java,v 1.16 2006/09/11 09:33:01 gerton Exp $
 */
public class BasicColourSampler implements ColourSampler {

    static final int maxIterations = 1000;

    public BasicColourSampler(Alignment tipColours, Tree tree) {

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

        nodePartials = new double[tree.getNodeCount()][colourCount];
    }

    public BasicColourSampler(TaxonList[] tipColours, Tree tree) {

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

        nodePartials = new double[tree.getNodeCount()][colourCount];
    }

    public int[] getLeafColourCounts() {
        return leafColourCounts;
    }

    /**
     * Colours the tree probabilistically with the given migration rates
     *
     * @param colourChangeMatrix the colour change rate parameters
     */
    public DefaultTreeColouring sampleTreeColouring(Tree tree, ColourChangeMatrix colourChangeMatrix, MetaPopulation mp) {
        DefaultTreeColouring colouring = new DefaultTreeColouring(2, tree);

        double[] N = mp.getPopulationSizes(0);

        double[] rootPartials = prune(tree, tree.getRoot(), colourChangeMatrix, N);

        // Sampling is conditional on data; so normalize by the probability of the
        // data under the proposal distribution
        double normalization = 0.0;
        for (int i = 0; i < rootPartials.length; i++) {
            normalization += colourChangeMatrix.getEquilibrium(i) * rootPartials[i];
        }

        sampleInternalNodes(tree, tree.getRoot(), colourChangeMatrix);

        sampleBranchColourings(colouring, tree, tree.getRoot(), colourChangeMatrix);

        double logP = calculateLogProbabilityDensity(colouring, tree, tree.getRoot(), colourChangeMatrix, N) - Math.log(normalization);

        colouring.setLogProbabilityDensity(logP);

        return colouring;
    }

    /**
     * @param node
     */
    private final int getColour(NodeRef node) {
        return nodeColours[node.getNumber()];
    }

    /**
     * @param node
     */
    private final void setColour(NodeRef node, int colour) {
        if (colour >= 0 && colour < colourCount) {
            nodeColours[node.getNumber()] = colour;
        } else {
            throw new IllegalArgumentException("colour value " + colour + " + is outside of range of colours, [0, " + Integer.toString(colourCount - 1) + "]");
        }
    }

    /*****************************************************************************************
     *
     *   Probability- and sampling-related code follows
     *
     */


    /**
     * Calculate probability of data at descendants from node, given a color at the node ('partials'),
     * by a Felsenstein-like pruning algorithm.  (First step in the color sampling algorithm)
     * Side effect: updates nodePartials[] for this node and all its descendants.
     *
     * @param node
     * @return the partials of this node
     */
    private final double[] prune(Tree tree, NodeRef node, ColourChangeMatrix mm, double[] N) {

        double[] p = new double[colourCount];

        if (tree.isExternal(node)) {
            p[getColour(node)] = 1.0;
            return p;
        }

        // Note: assuming binary tree!
        NodeRef leftChild = tree.getChild(node, 0);
        NodeRef rightChild = tree.getChild(node, 1);

        double[] left = prune(tree, leftChild, mm, N);
        double[] right = prune(tree, rightChild, mm, N);

        double nodeHeight = tree.getNodeHeight(node);

        double leftTime = nodeHeight - tree.getNodeHeight(tree.getChild(node, 0));
        double rightTime = nodeHeight - tree.getNodeHeight(tree.getChild(node, 1));

        for (int i = 0; i < p.length; i++) {

            double leftSum = 0.0;
            double rightSum = 0.0;

            // looping over colours in left and right children
            for (int j = 0; j < left.length; j++) {
                // forwardTimeEvolution conditions on the parent state i, i.e. time
                // runs in the natural direction (forward from parent to child)
                leftSum += left[j] * mm.forwardTimeEvolution(i, j, leftTime);
                rightSum += right[j] * mm.forwardTimeEvolution(i, j, rightTime);
            }
            p[i] = leftSum * rightSum;

            // Condition on the formal variable
            // (Removed - this didn't work; without it the sampler is robust)
            /*
            if (N != null) {
                p[i] /= N[i];
            }
            */
        }

        nodePartials[node.getNumber()] = p;

        return p;
    }


    /**
     * Samples internal node colours (from root to tips)
     * Requires the results from Felsenstein Backwards pruning, in nodePartials[] (see prune())
     * Side effect: updates nodeColours[]
     */
    private final void sampleInternalNodes(Tree tree, NodeRef node, ColourChangeMatrix mm) {

        double[] backward = nodePartials[node.getNumber()];
        double[] forward;

        if (tree.isRoot(node)) {

            forward = mm.getEquilibrium();

        } else {

            NodeRef parent = tree.getParent(node);
            int parentColour = getColour(parent);
            double time = tree.getNodeHeight(parent) - tree.getNodeHeight(node);
            forward = new double[backward.length];
            for (int i = 0; i < backward.length; i++) {
                forward[i] = mm.forwardTimeEvolution(parentColour, i, time);
            }

        }

        // Calculate the (unnormalized) probability for each colour, given the data and the
        // parent colour
        for (int i = 0; i < backward.length; i++) {
            forward[i] *= backward[i];
        }

        // Sample a colour
        int colour = MathUtils.randomChoicePDF(forward);
        setColour(node, colour);

        // Recursively sample down the tree
        for (int i = 0; i < tree.getChildCount(node); i++) {
            NodeRef child = tree.getChild(node, i);
            if (!tree.isExternal(child)) {
                sampleInternalNodes(tree, child, mm);
            }
        }
    }

    /**
     * Samples the colours on a tree branch, between node and its parent, conditional on the colour at these nodes.
     *
     * @param node the node above which to sample changes
     */
    private void sampleBranchColourings(DefaultTreeColouring colouring, Tree tree, NodeRef node, ColourChangeMatrix mm) {

        if (!tree.isRoot(node)) {
            NodeRef parent = tree.getParent(node);

            int parentColour = getColour(parent);
            int childColour = getColour(node);
            double parentHeight = tree.getNodeHeight(parent);
            double childHeight = tree.getNodeHeight(node);

            // Sample migration events on this branch, as a list of ColourChange-s
            DefaultBranchColouring history = sampleConditionalBranchColouring(parentColour, parentHeight, childColour, childHeight, mm);

            // Assign these migrations to the branch (attached to the child)
            colouring.setBranchColouring(node, history);
        }

        for (int i = 0; i < tree.getChildCount(node); i++) {
            sampleBranchColourings(colouring, tree, tree.getChild(node, i), mm);
        }
    }

    private DefaultBranchColouring sampleConditionalBranchColouring(int parentColour, double parentHeight,
                                                                    int childColour, double childHeight, ColourChangeMatrix mm) {

        DefaultBranchColouring history = new DefaultBranchColouring(parentColour, childColour);
        int currentColour;
        double currentHeight;
        int iterationsLeft = maxIterations;
        double time;

        // Reject until we get the child colour
        do {

            history.clear();
            currentColour = parentColour;
            currentHeight = parentHeight;

            // Sample events until we reach the child
            do {

                // Sample a waiting time
                double totalRate = -mm.getForwardRate(currentColour, currentColour);
                double U;

                do {
                    U = MathUtils.nextDouble();
                } while (U == 0.0);

                // Neat trick (Rasmus Nielsen):
                // If colours of parent and child differ, sample conditioning on at least 1 event
                if ((parentColour != childColour) && (history.getNumEvents() == 0)) {

                    double minU = Math.exp(-totalRate * (parentHeight - childHeight));
                    U = minU + U * (1.0 - minU);

                }

                // Calculate the waiting time, and update currentHeight
                time = -Math.log(U) / totalRate;
                currentHeight -= time;

                if (currentHeight > childHeight) {
                    // Not yet reached the child.  "Sample" an event
                    currentColour = 1 - currentColour;
                    // Add it to the list
                    history.addEvent(currentColour, currentHeight);
                }

            } while (currentHeight > childHeight);

            iterationsLeft -= 1;

        } while ((currentColour != childColour) && (iterationsLeft > 0));

        if (currentColour != childColour) {
            // Extreme migration rates may cause difficulty for the rejection sampler
            // Print a warning and add a bogus event somewhat near where you'd want it.
            double previousEventHeight = currentHeight + time;
            double finalEventHeight = childHeight + 0.01 * (previousEventHeight - childHeight);
            history.addEvent(childColour, finalEventHeight);

            System.out.println("dr.evolution.colouring.BranchColourSampler: failed to generate sample after " + maxIterations + " trials.");
            System.out.println(": parentColour=" + parentColour);
            System.out.println(": parentHeight=" + parentHeight);
            System.out.println(": childColour=" + childColour);
            System.out.println(": childHeight=" + childHeight);
            System.out.println(": migration rate 0->1 = " + mm.getForwardRate(0, 1));
            System.out.println(": migration rate 1->0 = " + mm.getForwardRate(1, 0));
        }

        return history;

    }


    /**
     * Calculates log probability density of the proposal colouring of the tree on the branch leading to this node,
     * and everything descending from it.
     */
    private final double calculateLogProbabilityDensity(TreeColouring colouring, Tree tree, NodeRef node, ColourChangeMatrix mm, double[] N) {

        double p = 1.0;

        if (tree.isRoot(node)) {

            p = mm.getEquilibrium(colouring.getNodeColour(node));

        } else {

            NodeRef parent = tree.getParent(node);
            BranchColouring history = colouring.getBranchColouring(node);   // note - it is attached to the child node

            int fromColour = colouring.getNodeColour(parent);
            double fromHeight = tree.getNodeHeight(parent);

            // Loop over all events, forward in time (i.e. down the tree)
            for (int i = 1; i <= history.getNumEvents(); i++) {

                // get colour below this node.
                int toColour = history.getForwardColourBelow(i);
                // get new height
                double toHeight = history.getForwardTime(i);
                // factor in the exit probability
                p *= Math.exp(-(fromHeight - toHeight) * (-mm.getForwardRate(fromColour, fromColour)));
                // and the event itself
                p *= mm.getForwardRate(fromColour, toColour);

                fromHeight = toHeight;
                fromColour = toColour;

            }

            // Include the exit probability on the branch from the last migration event to the child.
            double toHeight = tree.getNodeHeight(node);
            p *= Math.exp(-(fromHeight - toHeight) * (-mm.getForwardRate(fromColour, fromColour)));

            // Include the contribution of the formal variable (if this is an internal node)
            // (Removed)
            /*
            if (!tree.isExternal(node) && N != null) {
                p /= N[fromColour];
            }
            */
        }

        double logP = Math.log(p);

        for (int i = 0; i < tree.getChildCount(node); i++) {
            logP += calculateLogProbabilityDensity(colouring, tree, tree.getChild(node, i), mm, N);
        }

        return logP;
    }


    //
    // Calculates the Lebesgue measure of the space of migration events for the number of migration events
    // on each branch as specified by TreeColouring. 
    //
    public static final double calculateLogNormalization(TreeColouring colouring, Tree tree, NodeRef node) {

        final double arbitraryScaleFactor = 1.0;

        double logn = 0.0;

        if (!tree.isRoot(node)) {

            double norm = 1.0;
            double t = tree.getNodeHeight(tree.getParent(node)) - tree.getNodeHeight(node);
            int events = colouring.getBranchColouring(node).getNumEvents();

            for (int i = 1; i <= events; i++) {

                norm *= t / i;

            }

            logn = arbitraryScaleFactor * Math.log(norm);

        }

        for (int i = 0; i < tree.getChildCount(node); i++) {
            logn += calculateLogNormalization(colouring, tree, tree.getChild(node, i));
        }

        return logn;

    }


    public double getProposalProbability(TreeColouring treeColouring, Tree tree, ColourChangeMatrix colourChangeMatrix, MetaPopulation mp) {

        throw new IllegalArgumentException("Not implemented for BasicColourSampler; you can only use <ColouredOperator>s");

    }


    private final int colourCount;

    private int[] nodeColours;

    private double[][] nodePartials;

    private final int[] leafColourCounts;
}

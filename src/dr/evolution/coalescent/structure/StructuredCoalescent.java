/*
 * StructuredCoalescent.java
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

package dr.evolution.coalescent.structure;

import dr.evolution.colouring.ColourChangeMatrix;
import dr.evolution.colouring.TreeColouring;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;

/**
 * @author Alexei Drummond
 * @author Gerton Lunter
 * @author Andrew Rambaut
 * @version $Id: StructuredCoalescent.java,v 1.17 2006/09/11 09:33:01 gerton Exp $
 */
public class StructuredCoalescent {

    static double tinyTime = 1.0e-6;        // make sure rounding errors don't influence population sizes (for Bayesian Skyline)

    public StructuredCoalescent() {

    }

    /**
     * @param intervals
     * @return the log likelihood
     */
//    public double calculateLogLikelihood(TreeColouring treeColouring, StructuredIntervalList intervals, ColourChangeMatrix mm, double[] N) {
    public double calculateLogLikelihood(TreeColouring treeColouring, StructuredIntervalList intervals, ColourChangeMatrix mm, MetaPopulation mp) {

        double logL = 0;
        Tree tree = treeColouring.getTree();

        double _totalIntegratedRate = 0.0;  /* debugging */

        // include the probability of the leaf colouring
        for (int i = 0; i < tree.getExternalNodeCount(); i++) {

            NodeRef leafNode = tree.getExternalNode(i);
            logL += Math.log(mm.getEquilibrium(treeColouring.getNodeColour(leafNode)));

        }

        double currentTime = 0.0;

        // walk up the tree
        for (int i = 0; i < intervals.getIntervalCount(); i++) {

            // integral of exit rate over time interval
            double intensity = 0.0;

            // get length of current interval
            double time = intervals.getInterval(i);

            // calculate the integrated exit rate
            for (int j = 0; j < intervals.getPopulationCount(); j++) {

                int lineages = intervals.getLineageCount(i, j);

                // coalescent contribution
                //exitRate += lineages * (lineages-1) / (2*N[j]);

                intensity += (lineages * (lineages - 1) / 2.0) * mp.getIntegral(currentTime, currentTime + time, j);

                // migration contribution
                intensity += (-mm.getBackwardRate(j, j)) * lineages * time;
            }

            _totalIntegratedRate += intensity;

            //System.out.print(" LL: interv="+i+"\trate="+df.format(exitRate)+"\tt="+df.format(time)+"\tintrate="+ df.format(time*exitRate) + "\tsurv="+df.format(survivalProbability));

            double pointDensity = 0.0;
            Event event = intervals.getEvent(i);

            // update current time
            currentTime += time;
            //currentTime = event.time;

            switch (event.getType()) {

                case COALESCENT: {

                    if (!(event.getAboveColour() == event.getBelowColour())) {
                        throw new Error("Coalescent event changes colour");
                    }
                    // old version:
                    //pointDensity = 1.0 / N[event.aboveColour];

                    // Make sure that rounding errors do not influence the population size (relevant for
                    // non-continuous population models, e.g. piecewise constant models for the Bayesian skyline)
                    pointDensity = 1.0 / mp.getDemographic(currentTime - tinyTime, event.getAboveColour());

                    //System.out.println(" Coalescent; density="+pointDensity);

                }
                break;

                case MIGRATION: {

                    if (!(event.getAboveColour() != event.getBelowColour())) {
                        throw new Error("Colour-change event fails to change colour");
                    }
                    pointDensity = mm.getBackwardRate(event.getBelowColour(), event.getAboveColour());

                    //System.out.println(" Migration; density="+pointDensity);

                }
                break;
                case SAMPLE: {

                    pointDensity = 1.0;
                    //System.out.println(" Sample");
                }
                break;
            }

            logL += Math.log(pointDensity) - intensity;
        }

        /*
        System.out.println("Full likelihood "+logL);
        System.out.println("Integrated exit rate "+_totalIntegratedRate);
        */

        return logL;
    }

    /**
     * Calculate the (2-population) structured coalescent likelihood for the given tree,
     * tip-colors, population sizes and migration rates.
     * @param tree
     * @param initial the probability of being in population 0 at the tips
     * @param N the population sizes of the two populations
     * @param m the migration rates between the two populations
     */
    /*public double calculateLogLikelihood(Tree tree, double[] initial, double[] N, double[] m) {

        double logL = 0.0;

        // get the coalescent and add sample intervals
        TreeIntervals intervals = new TreeIntervals(tree);

        // the partial probs at all nodes.
        double[] p = new double[tree.getNodeCount()];

        // initialize the partials at all external nodes.
        // index of initial array should correspond to the index in external node list of tree.
        for (int i = 0; i < initial.length; i++) {
            p[i] = initial[i];
        }

        // the coalescence rates among all pairs of lineages
        double[][] cl = new double[tree.getNodeCount()][tree.getNodeCount()];

        // loop through intervals, updating partials
        for (int i = 0; i < intervals.getIntervalCount(); i++) {

            double ti = intervals.getInterval(i);
            List lineages = intervals.getLineages(i);

            // calculate new partials after time ti for all lineages in the current interval
            for (int j = 0; j < lineages.size(); j++) {

                NodeRef lineage = (NodeRef)lineages.get(j);
                int li = lineage.getNumber();

                // update the partial for node j
                p[li] = timeEvolution(p[li], ti);
            }

            // if this interval is a coalescent then calculate coalescent contribution to likelihood
            if (intervals.getIntervalType(i) == IntervalType.COALESCENT) {

                // calculate coalescence rates for all pairs of current lineages
                for (int j = 0; j < lineages.size(); j++) {

                    NodeRef lineage1 = (NodeRef)lineages.get(j);
                    int i1 = lineage1.getNumber();

                    for (int k = j+1; k < lineages.size(); k++) {

                        NodeRef lineage2 = (NodeRef)lineages.get(k);
                        int i2 = lineage2.getNumber();

                        double p0 = p[i1]*p[i2]; // probability both lineages in population 0
                        double p1 = (1-p[i1])*(1-p[i2]); // probability both lineages in population 1
                        double psame = p0 + p1; // probability both lineages i1 and i2 in same population
                        double pdiff = 1.0 - psame;

                        cl[i1][i2] = cl[i2][i1] = p0/N[0] + p1/N[1];
                    }
                }

                // do something with cl to calculate coalescent density

                // get nodes involved in coalescent event
                NodeRef parent = intervals.getCoalescentNode(i);
                int pi = parent.getNumber();
                NodeRef left = tree.getChild(parent,0);
                NodeRef right = tree.getChild(parent,1);
                int i1 = left.getNumber();
                int i2 = right.getNumber();

                // calculate the partial prob for the beginning of the coalesced lineage
                double p0 = p[i1]*p[i2]; // probability both lineages in population 0
                double p1 = (1-p[i1])*(1-p[i2]); // probability both lineages in population 1
                double psame = p0 + p1; // probability both lineages in same population

                p[pi] = p0 / psame; // re-normalized probability coalescent event in population 0

                // add the likelihood contribution from co-location probability psame
                // need to account for different coalescent rates here as well!
                logL += Math.log(psame);

            } else {
                // interval ended in a sample event -- nothing left to do
            }

            // at this point all p values have been updated for the start of the next coalescent interval.
        }

        return logL;
    } */

    /*
    *//**
     * @return the probability of being in population 0 after time t, starting from prob p
     *//*
    private double timeEvolution(double p, double t) {
        throw new RuntimeException("Not implemented");
    }*/
}

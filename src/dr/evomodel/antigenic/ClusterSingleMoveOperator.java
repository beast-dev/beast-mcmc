/*
 * ClusterSingleMoveOperator.java
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

package dr.evomodel.antigenic;

import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;

/**
 * An operator to take a single element from one cluster and move it to a new or different cluster.
 * Based on Ed Baskerville's implementation in Spatial Guilds in the Serengeti Food Web Revealed by
 * a Bayesian Group Model.
 *
 * @author Trevor Bedford
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id$
 */
public class ClusterSingleMoveOperator extends SimpleMCMCOperator {
    public final static boolean DEBUG = false;

    public final static String CLUSTER_SINGLE_MOVE_OPERATOR = "clusterSingleMoveOperator";

    private final int N; // the number of items
    private int K; // the number of occupied clusters
    private final Parameter allocationParameter;

    public ClusterSingleMoveOperator(Parameter allocationParameter, double weight) {
        this.allocationParameter = allocationParameter;
        this.N = allocationParameter.getDimension();

        setWeight(weight);
    }


    /**
     * @return the parameter this operator acts on.
     */
    public Parameter getParameter() {
        return (Parameter) allocationParameter;
    }

    /**
     * @return the Variable this operator acts on.
     */
    public Variable getVariable() {
        return allocationParameter;
    }

    /**
     * change the parameter and return the hastings ratio.
     */
    public final double doOperation() {

        // get a copy of the allocations to work with...
        // allocations are: element X -> cluster Y
        int[] allocations = new int[allocationParameter.getDimension()];

        // construct cluster occupancy vector excluding the selected item and count
        // the unoccupied clusters.
        // occupancy is: cluster Y -> element count
        int[] occupancy = new int[N];

        int K = 0; // k = number of unoccupied clusters
        for (int i = 0; i < allocations.length; i++) {
            allocations[i] = (int) allocationParameter.getParameterValue(i);
            occupancy[allocations[i]] += 1;
            if (occupancy[allocations[i]] == 1) { // first item in cluster
                K++;
            }
        }

        // log hastings ratio
        double hastings = 0.0;

        // pick element to move
        int element = MathUtils.nextInt(N);
        int elementAssignment = allocations[element];
        int elementClusterSize = occupancy[elementAssignment]; // cluster size before move

        // pick second element as target, must be different from first element
        int target = MathUtils.nextInt(N);
        while (element == target) {
            target = MathUtils.nextInt(N);
        }
        int targetAssignment = allocations[target];
        int targetClusterSize = occupancy[targetAssignment]; // cluster size before move

        // if allocation of element differs from allocation of target
        // change element allocation to match target allocation
        if (elementAssignment != targetAssignment) {

            allocations[element] = targetAssignment;
            allocationParameter.setParameterValue(element, targetAssignment);

            if (elementClusterSize > 1) {
                // adjusting Hastings's ratio for differences in cluster size
                hastings = Math.log(elementClusterSize - 1) - Math.log(targetClusterSize);
            }

            if (DEBUG) {
                System.err.println("Move element " + element + " from cluster " + elementAssignment + " to cluster " + targetAssignment);
            }

        }

        // if allocation of element matches allocation of target
        // move element to a new cluster
        else {

            // find random unoccupied cluster
            int clusterIndex = 0;
            int targetIndex = MathUtils.nextInt(N-K);
            int newCluster = 0;
            for (int i = 0; i < N; i++) {
                if (occupancy[i] == 0){
                    if (clusterIndex == targetIndex) {
                        newCluster = i;
                    }
                    clusterIndex++;
                }
            }

            while (occupancy[newCluster] > 0) {
                newCluster++;
            }

            // move the element to this cluster
            allocations[element] = newCluster;
            allocationParameter.setParameterValue(element, newCluster);

            if (DEBUG) {
                System.err.println("Move element " + element + " from cluster " + elementAssignment + " to new cluster " + newCluster);
            }

        }

        // return log Hastings' ratio
        return hastings;
    }


    //MCMCOperator INTERFACE
    public final String getOperatorName() {
        return CLUSTER_SINGLE_MOVE_OPERATOR +"(" + allocationParameter.getId() + ")";
    }

    public final void optimize(double targetProb) {

        throw new RuntimeException("This operator cannot be optimized!");
    }

    public boolean isOptimizing() {
        return false;
    }

    public void setOptimizing(boolean opt) {
        throw new RuntimeException("This operator cannot be optimized!");
    }

    public double getMinimumAcceptanceLevel() {
        return 0.1;
    }

    public double getMaximumAcceptanceLevel() {
        return 0.4;
    }

    public double getMinimumGoodAcceptanceLevel() {
        return 0.20;
    }

    public double getMaximumGoodAcceptanceLevel() {
        return 0.30;
    }

    public String getPerformanceSuggestion() {
        if (Utils.getAcceptanceProbability(this) < getMinimumAcceptanceLevel()) {
            return "";
        } else if (Utils.getAcceptanceProbability(this) > getMaximumAcceptanceLevel()) {
            return "";
        } else {
            return "";
        }
    }

    public String toString() {
        return getOperatorName();
    }


    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return CLUSTER_SINGLE_MOVE_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

            Parameter allocationParameter = (Parameter) xo.getChild(Parameter.class);

            return new ClusterSingleMoveOperator(allocationParameter, weight);

        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "An operator that moves single elements between clusters.";
        }

        public Class getReturnType() {
            return ClusterSingleMoveOperator.class;
        }


        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                new ElementRule(Parameter.class)
        };
    };

    public int getStepCount() {
        return 1;
    }

}

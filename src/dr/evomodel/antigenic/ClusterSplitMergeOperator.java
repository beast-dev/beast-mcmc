/*
 * ClusterSplitMergeOperator.java
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

import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.inference.operators.MCMCOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;

/**
 * An operator to split or merge clusters.
 *
 * @author Andrew Rambaut
 * @author Marc Suchard
 * @version $Id: DirichletProcessGibbsOperator.java,v 1.16 2005/06/14 10:40:34 rambaut Exp $
 */
public class ClusterSplitMergeOperator extends SimpleMCMCOperator {
    public final static boolean DEBUG = false;

    public final static String CLUSTER_SPLIT_MERGE_OPERATOR = "clusterSplitMergeOperator";

    private final int N; // the number of items
    private int K; // the number of occupied clusters
    private final Parameter allocationParameter;
    private final MatrixParameter clusterLocations;

    public ClusterSplitMergeOperator(Parameter allocationParameter,
                                     MatrixParameter clusterLocations,
                                     double weight) {
        this.allocationParameter = allocationParameter;
        this.clusterLocations = clusterLocations;
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
        int[] allocations = new int[allocationParameter.getDimension()];

        // construct cluster occupancy vector excluding the selected item and count
        // the unoccupied clusters.
        int[] occupancy = new int[N];

        int[] occupiedIndices = new int[N];

        // for testing, set these to -1 to force out of bounds exception if
        // used but not set
        for (int i = 0; i < occupiedIndices.length; i++) {
            occupiedIndices[i] = -1;
        }

        int K = 0; // k = number of unoccupied clusters
        for (int i = 0; i < allocations.length; i++) {
            allocations[i] = (int) allocationParameter.getParameterValue(i);
            occupancy[allocations[i]] += 1;
            if (occupancy[allocations[i]] == 1) { // first item in cluster
                occupiedIndices[K] = allocations[i];
                K++;
            }
        }

        // Container for split/merge random variable (only 2 draws in 2D)
        int paramDim = clusterLocations.getParameter(0).getDimension();
        double[] splitDraw = new double[paramDim]; // Need to keep these for computing MHG ratio
        double scale = 1.0; // TODO make tunable
        double newClusterProb = 0.5;  // TODO Make tunable

        // always split when K = 1, always merge when K = N, otherwise 50:50
        boolean doSplit = K == 1 || (K != N && MathUtils.nextBoolean());

        if (doSplit) {
            // Split operation

            int cluster1;
            do {
                // pick an occupied cluster
                cluster1 = occupiedIndices[MathUtils.nextInt(K)];

                // For reversibility, merge step requires that both resulting clusters are occupied,
                // so we should resample until condition is true
            } while (occupancy[cluster1] == 0);

            // find the first unoccupied cluster
            int cluster2 = 0;
            while (occupancy[cluster2] > 0) {
                cluster2 ++;
            }

            int oldCount = occupancy[cluster1];
            do {
                occupancy[cluster1] = 0;
                occupancy[cluster2] = 0;
                for (int i = 0; i < allocations.length; i++) {
                    if (allocations[i] == cluster1 || allocations[i] == cluster2) {
                        boolean putInNewCluster = MathUtils.nextDouble() < newClusterProb;
                        if (putInNewCluster) {
                            allocations[i] = cluster2;
                            occupancy[cluster2]++;
                        } else {
                            allocations[i] = cluster1;
                            occupancy[cluster1]++;
                        }
                    }
                }
            } while (occupancy[cluster1] != 0 && occupancy[cluster2] != 0);

            K++;

            // set both clusters to a location based on the first cluster with some random jitter...
            Parameter param1 = clusterLocations.getParameter(cluster1);
            Parameter param2 = clusterLocations.getParameter(cluster2);

            double[] loc = param1.getParameterValues();
            for (int dim = 0; dim < param1.getDimension(); dim++) {
                splitDraw[dim] = MathUtils.nextGaussian();
                param1.setParameterValue(dim, loc[dim] + (splitDraw[dim] * scale));
                param2.setParameterValue(dim, loc[dim] - (splitDraw[dim] * scale)); // Move in opposite direction
            }

            if (DEBUG) {
                System.err.println("Split: " + (oldCount - occupancy[cluster1]) + " items from cluster " + cluster1 + " to create cluster " + cluster2);
            }
        } else {
            // Merge operation

            // pick 2 occupied clusters
            int cluster1 = occupiedIndices[MathUtils.nextInt(K)];
            int cluster2;

            do {
                cluster2 = occupiedIndices[MathUtils.nextInt(K)];

                // resample until cluster1 != cluster2 to maintain reversibility, because split assumes they are different
            } while (cluster1 == cluster2);

            for (int i = 0; i < allocations.length; i++) {
                if (allocations[i] == cluster2) {
                    allocations[i] = cluster1;

                    // keep occupancy up to date (remove if not need)
                    occupancy[cluster1] ++;
                    occupancy[cluster2] --;
                }
            }

            K--;

            // set the merged cluster to the mean location of the two original clusters
            Parameter loc1 = clusterLocations.getParameter(cluster1);
            Parameter loc2 = clusterLocations.getParameter(cluster2);
            for (int dim = 0; dim < loc1.getDimension(); dim++) {
                double average = (loc1.getParameterValue(dim) + loc2.getParameterValue(dim)) / 2.0;
                splitDraw[dim] = (loc1.getParameterValue(dim) - average) / scale; // Record that the reverse step would need to draw
                loc1.setParameterValue(dim, average);

                // Consider loc2 as the extra dimensions for dimension-matching
                // On second thought, maybe not a good idea
//                loc2.setParameterValue(dim, splitDraw[dim]);
            }
            if (DEBUG) {
                System.err.println("Merge: " + occupancy[cluster1] + "items into cluster " + cluster1 + " from " + cluster2);
            }

        }

        // set the final allocations (only for those that have changed)
        for (int i = 0; i < allocations.length; i++) {
            int k = (int) allocationParameter.getParameterValue(i);
            if (allocations[i] != k) {
                allocationParameter.setParameterValue(i, allocations[i]);
            }
        }

        // todo the Hastings ratio
        return 0.0;
    }


    //MCMCOperator INTERFACE
    public final String getOperatorName() {
        return CLUSTER_SPLIT_MERGE_OPERATOR +"(" + allocationParameter.getId() + ")";
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
        public final static String CHI = "chi";
        public final static String LIKELIHOOD = "likelihood";

        public String getParserName() {
            return CLUSTER_SPLIT_MERGE_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            double weight = xo.getDoubleAttribute(MCMCOperator.WEIGHT);

            Parameter allocationParameter = (Parameter) xo.getChild(Parameter.class);

            MatrixParameter locationsParameter = (MatrixParameter) xo.getElementFirstChild("locations");

            return new ClusterSplitMergeOperator(allocationParameter, locationsParameter, weight);

        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "An operator that splits and merges clusters.";
        }

        public Class getReturnType() {
            return ClusterSplitMergeOperator.class;
        }


        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private final XMLSyntaxRule[] rules = {
                AttributeRule.newDoubleRule(MCMCOperator.WEIGHT),
                new ElementRule(Parameter.class),
                new ElementRule("locations", new XMLSyntaxRule[] {
                        new ElementRule(MatrixParameter.class)
                })
        };
    };

    public int getStepCount() {
        return 1;
    }

}

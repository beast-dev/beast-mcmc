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

        int X = N; // X = number of unoccupied clusters
        for (int i = 0; i < allocations.length; i++) {
            allocations[i] = (int) allocationParameter.getParameterValue(i);
            occupancy[allocations[i]] += 1;
            if (occupancy[allocations[i]] == 1) { // first item in cluster
                X -= 1; // one fewer unoccupied
            }
        }
        K = N - X;

        // Container for split/merge random variable (only 2 draws in 2D)
        int paramDim = clusterLocations.getParameter(0).getDimension();
        double[] splitDraw = new double[paramDim]; // Need to keep these for computing MHG ratio
        double scale = 1.0; // TODO make tunable

        if (MathUtils.nextDouble() < 0.5) {
            // Split operation

            int cluster1;
            int cluster2;

            do {
                // pick an occupied cluster
                cluster1 = MathUtils.nextInt(K);
                cluster2 = K; // next available unoccupied cluster

                for (int i = 0; i < allocations.length; i++) {
                    if (allocations[i] == cluster1) {
                        if (MathUtils.nextDouble() < 0.5) {
                            allocations[i] = cluster2;
                            occupancy[cluster1] --;
                            occupancy[cluster2] ++;
                        }
                    }
                }

                // For reversibility, merge step requires that both resulting clusters are occupied,
                // so we should resample until condition is true
            } while (occupancy[cluster1] == 0);

            // set both clusters to a location based on the first cluster with some random jitter...
            Parameter param1 = clusterLocations.getParameter(cluster1);
            Parameter param2 = clusterLocations.getParameter(cluster2);

            double[] loc = param1.getParameterValues();
            for (int dim = 0; dim < param1.getDimension(); dim++) {
                splitDraw[dim] = MathUtils.nextGaussian();
                param1.setParameterValue(dim, loc[dim] + (splitDraw[dim] * scale));
                param2.setParameterValue(dim, loc[dim] - (splitDraw[dim] * scale)); // Move in opposite direction
            }

        } else {
            // Merge operation

            // pick 2 occupied clusters
            int cluster1 = MathUtils.nextInt(K);
            int cluster2;

            do {
                cluster2 = MathUtils.nextInt(K);

                // resample until cluster1 != cluster2 to maintain reversibility, because split assumes they are different
            } while (cluster1 == cluster2);

            if (cluster1 > cluster2) {
                // swap the cluster indices to keep the destination cluster lower
                int tmp = cluster1;
                cluster1 = cluster2;
                cluster2 = tmp;
            }

            for (int i = 0; i < allocations.length; i++) {
                if (allocations[i] == cluster2) {
                    allocations[i] = cluster1;
                    occupancy[cluster1] ++;
                    occupancy[cluster2] --;
                }
            }

            // set the merged cluster to the mean location of the two original clusters
            Parameter loc1 = clusterLocations.getParameter(cluster1);
            Parameter loc2 = clusterLocations.getParameter(cluster2);
            for (int dim = 0; dim < loc1.getDimension(); dim++) {
                double average = (loc1.getParameterValue(dim) + loc2.getParameterValue(dim)) / 2.0;
                splitDraw[dim] = (loc1.getParameterValue(dim) - average) / scale; // Record that the reverse step would need to draw
                loc1.setParameterValue(dim, average);
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

            MatrixParameter locationsParameter = (MatrixParameter) xo.getChild(MatrixParameter.class);

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
                new ElementRule(MatrixParameter.class)
        };
    };

    public int getStepCount() {
        return 1;
    }

}

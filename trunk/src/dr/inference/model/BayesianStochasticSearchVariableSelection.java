package dr.inference.model;

import cern.colt.bitvector.BitVector;

import cern.colt.bitvector.BitVector;
import dr.math.MathUtils;

/**
 * @author Marc Suchard
 */

public interface BayesianStochasticSearchVariableSelection {

    public Parameter getIndicators();

    public boolean validState();

    public class Utils {
        
    public static boolean connectedAndWellConditioned(double[] probability) {
            for(int i=0; i<probability.length; i++) {
                if(probability[i] == 0 || probability[i] > 1)
                    return false;
            }
            return true;
        }

        public static void randomize(Parameter indicators,int dim, boolean reversible) {
            do {
                for (int i = 0; i < indicators.getDimension(); i++)
                    indicators.setParameterValue(i,
                            (MathUtils.nextDouble() < 0.5) ? 0.0 : 1.0);
            } while (!(isStronglyConnected(indicators.getParameterValues(),
                    dim, reversible)));
        }

        /* Determines if the graph is strongly connected, such that there exists
        * a directed path from any vertex to any other vertex
        *
        */
        public static boolean isStronglyConnected(double[] indicatorValues, int dim, boolean reversible) {
            BitVector visited = new BitVector(dim);
            boolean connected = true;
            for (int i = 0; i < dim && connected; i++) {
                visited.clear();
                depthFirstSearch(i, visited, indicatorValues, dim, reversible);
                connected = visited.cardinality() == dim;
            }
            return connected;
        }

        private static boolean hasEdge(int i, int j, double[] indicatorValues,
                                      int dim, boolean reversible) {
            return i != j && indicatorValues[getEntry(i, j, dim, reversible)] == 1;
        }

        private static int getEntry(int i, int j, int dim, boolean reversible) {
            if (reversible && j > i)
                return getEntry(j,i,dim,false);            

            int entry = i * (dim - 1) + j;
            if (j > i)
                entry--;
            return entry;
        }

        private static void depthFirstSearch(int node, BitVector visited, double[] indicatorValues,
                                            int dim, boolean reversible) {
            visited.set(node);
            for (int v = 0; v < dim; v++) {
                if (hasEdge(node, v, indicatorValues, dim, reversible) && !visited.get(v))
                    depthFirstSearch(v, visited, indicatorValues, dim, reversible);
            }
        }
    }
}

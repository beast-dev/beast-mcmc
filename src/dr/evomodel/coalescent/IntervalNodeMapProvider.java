package dr.evomodel.coalescent;

import java.util.ArrayList;
import java.util.Arrays;

import dr.evolution.tree.Tree;
import dr.util.ComparableDouble;
import dr.util.HeapSort;

public interface IntervalNodeMapProvider {


    int[] getIntervalsForNode(int nodeNumber);
    int[] getNodeNumbersForInterval(int interval);
    double[] sortByNodeNumbers(double[] byIntervalOrder);

    interface IntervalNodeMapping extends IntervalNodeMapProvider {
        void addNode(int nodeNumbe);
        void setIntervalStartIndices(int intervalCount);
        void initializeMaps();
        public void storeMapping();
        public void restoreMapping();
    }
    class Default implements IntervalNodeMapping {
        private int[] nodeNumbersInIntervals;
        private int[] intervalStartIndices;
        private int[] intervalNumberOfNodes;

        private int[] storedNodeNumbersInIntervals;
        private int[] storedIntervalStartIndices;
        private int[] storedIntervalNumberOfNodes;

        private int nextIndex = 0;
        private int nIntervals;
        private Tree tree;
        private final int maxIndicesPerNode = 3;

        public Default (int maxIntervalCount, Tree tree) {
            nodeNumbersInIntervals = new int[maxIndicesPerNode * maxIntervalCount];
            storedNodeNumbersInIntervals = new int[maxIndicesPerNode * maxIntervalCount];

            intervalStartIndices = new int[maxIntervalCount];
            storedIntervalStartIndices = new int[maxIntervalCount];

            intervalNumberOfNodes = new int[maxIndicesPerNode * maxIntervalCount];
            storedIntervalNumberOfNodes = new int[maxIndicesPerNode * maxIntervalCount];
            this.tree = tree;
        }

        public void addNode(int nodeNumber) {
            nodeNumbersInIntervals[nextIndex] = nodeNumber;
            nextIndex++;
        }

        private void mapNodeInterval(int nodeNumber, int intervalNumber) {
            int index = 0;
            while(index < maxIndicesPerNode) {
                if (intervalNumberOfNodes[maxIndicesPerNode * nodeNumber + index] == -1) {
                    intervalNumberOfNodes[maxIndicesPerNode * nodeNumber + index] = intervalNumber;
                    break;
                } else {
                    index++;
                }
            }
            if (index == maxIndicesPerNode) {
                throw new RuntimeException("The node appears in more than" + maxIndicesPerNode + " intervals!");
            }
//                if (intervalNumberOfNodes[maxIndicesPerNode * nodeNumber] == -1 || intervalNumberOfNodes[maxIndicesPerNode * nodeNumber] == intervalNumber) {
//                    intervalNumberOfNodes[3 * nodeNumber] = intervalNumber;
//                } else if (intervalNumberOfNodes[2 * nodeNumber + 1] == -1) {
//                    intervalNumberOfNodes[3 * nodeNumber + 1] = intervalNumber;
//                } else {
//                    double[] testIntervals = new double[nIntervals];
//                    for (int i = 0; i < nIntervals - 1; i++) {
//                        testIntervals[i] = tree.getNodeHeight(tree.getNode(nodeNumbersInIntervals[intervalStartIndices[i + 1]]))
//                                - tree.getNodeHeight(tree.getNode(nodeNumbersInIntervals[intervalStartIndices[i]]));
//                    }
//                    throw new RuntimeException("The node appears in more than two intervals!");
//                }
        }

        public void setIntervalStartIndices(int intervalCount) {

            if (nodeNumbersInIntervals[nextIndex - 1] == nodeNumbersInIntervals[nextIndex - 2]) {
                nodeNumbersInIntervals[nextIndex - 1] = 0;
                nextIndex--;
            }

            int index = 1;
            mapNodeInterval(nodeNumbersInIntervals[0], 0);

            for (int i = 1; i < intervalCount; i++) {

                while(nodeNumbersInIntervals[index] != nodeNumbersInIntervals[index - 1]) {
                    mapNodeInterval(nodeNumbersInIntervals[index], i - 1);
                    index++;
                }

                intervalStartIndices[i] = index;
                mapNodeInterval(nodeNumbersInIntervals[index], i);
                index++;

            }

            while(index < nextIndex) {
                mapNodeInterval(nodeNumbersInIntervals[index], intervalCount - 1);
                index++;
            }

            nIntervals = intervalCount;
        }

        public void initializeMaps() {
            Arrays.fill(intervalNumberOfNodes, -1);
            Arrays.fill(intervalStartIndices, 0);
            Arrays.fill(nodeNumbersInIntervals, 0);
            nextIndex = 0;
        }

        @Override
        public int[] getIntervalsForNode(int nodeNumber) {
            int nonZeros = 0;
            while(intervalNumberOfNodes[maxIndicesPerNode * nodeNumber + nonZeros] != -1) {
                nonZeros++;
            }
            int[] result = new int[nonZeros];
            for(int i = 0; i < nonZeros; i++) {
                result[i] = intervalNumberOfNodes[maxIndicesPerNode * nodeNumber + i];
            }
            return result;
//                if(intervalNumberOfNodes[maxIndicesPerNode * nodeNumber + 1] == -1) {
//                    return new int[]{intervalNumberOfNodes[maxIndicesPerNode * nodeNumber]};
//                } else {
//                    return new int[]{intervalNumberOfNodes[maxIndicesPerNode * nodeNumber], intervalNumberOfNodes[maxIndicesPerNode * nodeNumber + 1]};
//                }
        }

        @Override
        public int[] getNodeNumbersForInterval(int interval) {
            assert(interval < nIntervals);

            final int startIndex = intervalStartIndices[interval];
            int endIndex;
            if (interval == nIntervals - 1) {
                endIndex = nextIndex - 1;
            } else {
                endIndex = intervalStartIndices[interval + 1] - 1;
            }

            int[] nodeNumbers = new int[endIndex - startIndex + 1];

            for (int i = 0; i < endIndex - startIndex + 1; i++) {
                nodeNumbers[i] = nodeNumbersInIntervals[startIndex + i];
            }
            return nodeNumbers;
        }

        @Override
        public double[] sortByNodeNumbers(double[] byIntervalOrder) {
            double[] sortedValues = new double[byIntervalOrder.length];
            int[] nodeIndices = new int[byIntervalOrder.length];
            ArrayList<ComparableDouble> mappedIntervals = new ArrayList<ComparableDouble>();
            for (int i = 0; i < nodeIndices.length; i++) {
                mappedIntervals.add(new ComparableDouble(getIntervalsForNode(i + tree.getExternalNodeCount())[0]));
            }
            HeapSort.sort(mappedIntervals, nodeIndices);
            for (int i = 0; i < nodeIndices.length; i++) {
                sortedValues[nodeIndices[i]] = byIntervalOrder[i];
            }
            return sortedValues;
        }

                    /**
             * Additional state information, outside of the sub-model is stored by this call.
             */
            public void storeMapping() {
                System.arraycopy(nodeNumbersInIntervals,0,storedNodeNumbersInIntervals,0,nodeNumbersInIntervals.length);
                System.arraycopy(intervalNumberOfNodes,0,storedIntervalNumberOfNodes,0,intervalNumberOfNodes.length);
                System.arraycopy(intervalStartIndices,0,storedIntervalStartIndices,0,intervalStartIndices.length);
            }

            /**
             * After this call the model is guaranteed to have returned its extra state information to
             * the values coinciding with the last storeState call.
             * Sub-models are handled automatically and do not need to be considered in this method.
             */
             public void restoreMapping() {
                int[] tmp = storedNodeNumbersInIntervals;
                storedNodeNumbersInIntervals = nodeNumbersInIntervals;
                nodeNumbersInIntervals = tmp;

                int[] tmp2 = storedIntervalNumberOfNodes;
                storedIntervalNumberOfNodes=intervalNumberOfNodes;
                intervalNumberOfNodes = tmp2;

                int[] tmp3= storedIntervalStartIndices;
                storedIntervalStartIndices = intervalStartIndices;
                intervalStartIndices =tmp3;
            }
    }

    class None implements IntervalNodeMapping {

        @Override
        public void addNode(int nodeNumber) {
            // Do nothing
        }

        @Override
        public void setIntervalStartIndices(int intervalCount) {
            // Do nothing
        }

        @Override
        public void initializeMaps() {
            // Do nothing
        }

        @Override
        public int[] getIntervalsForNode(int nodeNumber) {
            throw new RuntimeException("No intervalNodeMapping available. This function should not be called.");
        }

        @Override
        public int[] getNodeNumbersForInterval(int interval) {
            throw new RuntimeException("No intervalNodeMapping available. This function should not be called.");
        }

        @Override
        public double[] sortByNodeNumbers(double[] byIntervalOrder) {
            throw new RuntimeException("No intervalNodeMapping available. This function should not be called.");
        }
        public void storeMapping() {
            // Do nothing
        }
        public void restoreMapping() {
            // Do nothing
        }
    }
}





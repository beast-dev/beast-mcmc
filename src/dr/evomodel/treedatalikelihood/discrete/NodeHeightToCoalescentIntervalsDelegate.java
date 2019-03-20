package dr.evomodel.treedatalikelihood.discrete;

import dr.evomodel.coalescent.GMRFSkyrideLikelihood;
import dr.evomodel.coalescent.OldAbstractCoalescentLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.inference.model.Bounds;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

public class NodeHeightToCoalescentIntervalsDelegate extends AbstractNodeHeightTransformDelegate {

    private GMRFSkyrideLikelihood skyrideLikelihood;
    private Parameter coalescentIntervals;
    private OldAbstractCoalescentLikelihood.IntervalNodeMapping intervalNodeMapping;

    public NodeHeightToCoalescentIntervalsDelegate(TreeModel treeModel,
                                                   Parameter nodeHeights,
                                                   GMRFSkyrideLikelihood skyrideLikelihood) {

        super(treeModel, nodeHeights);

        this.skyrideLikelihood = skyrideLikelihood;
        this.intervalNodeMapping = skyrideLikelihood.getIntervalNodeMapping();
        this.coalescentIntervals = createProxyForCoalescentIntervals();
        this.coalescentIntervals.addBounds(new NodeHeightToCoalescentIntervalsDelegate.CoalescentIntervalBounds());
        addVariable(coalescentIntervals);

        this.proxyValuesKnown = false;
    }

    @Override
    double[] transform(double[] values, int from, int to) {
        setNodeHeights(values);
        skyrideLikelihood.setupCoalescentIntervals();
        return coalescentIntervals.getParameterValues();
    }

    @Override
    double[] inverse(double[] values, int from, int to) {
        if (values.length != coalescentIntervals.getDimension()) {
            throw new RuntimeException("Dimension mismatch!");
        }

        double currentHeight = 0.0;

        for (int i = 0; i < values.length; i++) {
            int[] nodeNumbers = intervalNodeMapping.getNodeNumbersForInterval(i);
            currentHeight += values[i];
            TreeModel.Node node = (TreeModel.Node) tree.getNode(nodeNumbers[nodeNumbers.length - 1]);
            node.heightParameter.setParameterValueQuietly(0, currentHeight);
        }
        tree.pushTreeChangedEvent();
        return nodeHeights.getParameterValues();
    }

    @Override
    String getReport() {
        return null;
    }

    @Override
    Parameter getParameter() {
        return coalescentIntervals;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {

    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        proxyValuesKnown = false;
    }

    private boolean proxyValuesKnown;

    private Parameter createProxyForCoalescentIntervals() {
        return new Parameter.Proxy("coalescentIntervals",
                skyrideLikelihood.getCoalescentIntervalDimension()) {

            private double[] proxy;
            private Bounds<Double> bounds = null;

            {
                proxy = new double[dim];
            }

            @Override
            public double getParameterValue(int dim) {
                updateCoalescentIntervals();
                return proxy[dim];
            }

            @Override
            public void setParameterValue(int dim, double value) { // This function is very expensive, avoid repeated calls
                setParameterValueQuietly(dim, value);
                updateAllNodeHeights();
                proxyValuesKnown = false;
            }

            @Override
            public void setParameterValueQuietly(int dim, double value) {
                proxy[dim] = value;
            }

            @Override
            public void setParameterValueNotifyChangedAll(int dim, double value) {
                setParameterValue(dim, value);
            }

            @Override
            public void addBounds(Bounds<Double> boundary) {
                if (bounds == null) {
                    bounds = boundary;
                } else {
                    throw new RuntimeException("Not yet implemented for multiple use of addBounds.");
                }
            }

            @Override
            public Bounds<Double> getBounds() {
                if (bounds == null) {
                    throw new NullPointerException(getParameterName() + " parameter: Bounds not set");
                }
                return bounds;
            }

            @Override
            public void fireParameterChangedEvent(int index, Parameter.ChangeType type) {
                updateAllNodeHeights();
            }

            private void updateCoalescentIntervals() {
                if (!proxyValuesKnown) {
                    System.arraycopy(skyrideLikelihood.getCoalescentIntervals(), 0,
                            proxy, 0, proxy.length);
                    ((NodeHeightToCoalescentIntervalsDelegate.CoalescentIntervalBounds) getBounds()).setupBounds();
                    proxyValuesKnown = true;
                }
            }

            private void updateAllNodeHeights() {
                updateCoalescentIntervals();
                inverse(proxy, 0, proxy.length);
            }
        };
    }

    private class CoalescentIntervalBounds implements Bounds<Double> {

        private final int dim;
        private double[] upperBounds;
        private double[] lowerBounds;


        public CoalescentIntervalBounds() {
            this.dim = coalescentIntervals.getDimension();
            this.upperBounds = new double[dim];
            this.lowerBounds = new double[dim];
            setupBounds();
        }

        public void setupBounds() {
            if (!proxyValuesKnown) {

                double maxIncrease = Double.POSITIVE_INFINITY;
                upperBounds[dim - 1] = Double.POSITIVE_INFINITY;
                for (int i = dim - 2; i > -1; i--) {
                    int[] nodeNumbers = intervalNodeMapping.getNodeNumbersForInterval(i + 1);
                    final double currentIncrease = tree.getNodeHeight(tree.getNode(nodeNumbers[1])) - tree.getNodeHeight(tree.getNode(nodeNumbers[0]));
                    if (maxIncrease > currentIncrease) {
                        maxIncrease = currentIncrease;
                    }
                    upperBounds[i] = skyrideLikelihood.getCoalescentInterval(i) + maxIncrease;
                }

                double maxDecrease = Double.POSITIVE_INFINITY;
                for (int i = dim - 1; i > -1; i--) {
                    int[] nodeNumbers = intervalNodeMapping.getNodeNumbersForInterval(i);
                    final double currentDecrease = tree.getNodeHeight(tree.getNode(nodeNumbers[nodeNumbers.length - 1]))
                            - tree.getNodeHeight(tree.getNode(nodeNumbers[nodeNumbers.length - 2]));
                    if (maxDecrease > currentDecrease) {
                        maxDecrease = currentDecrease;
                    }
                    lowerBounds[i] = skyrideLikelihood.getCoalescentInterval(i) - maxDecrease;
                }
            }
        }


        @Override
        public Double getUpperLimit(int dimension) {
            setupBounds();
            return upperBounds[dimension];
        }

        @Override
        public Double getLowerLimit(int dimension) {
            setupBounds();
            return lowerBounds[dimension];
        }

        @Override
        public int getBoundsDimension() {
            return dim;
        }
    }
}

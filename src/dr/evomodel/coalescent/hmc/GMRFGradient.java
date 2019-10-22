package dr.evomodel.coalescent.hmc;

import dr.evolution.coalescent.TreeIntervals;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.coalescent.GMRFMultilocusSkyrideLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.HessianWrtParameterProvider;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.util.ComparableDouble;
import dr.util.HeapSort;
import dr.xml.Reportable;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Mandev Gill
 */
public class GMRFGradient implements GradientWrtParameterProvider, HessianWrtParameterProvider, Reportable, Loggable {

    private final GMRFMultilocusSkyrideLikelihood skygridLikelihood;
    private final WrtParameter wrtParameter;
    private final Parameter parameter;

    public GMRFGradient(GMRFMultilocusSkyrideLikelihood skygridLikelihood,
                        WrtParameter wrtParameter) {
        this.skygridLikelihood = skygridLikelihood;
        this.wrtParameter = wrtParameter;
        parameter = wrtParameter.getParameter(skygridLikelihood);
    }

    @Override
    public Likelihood getLikelihood() {
        return skygridLikelihood;
    }

    @Override
    public Parameter getParameter() {
        return parameter;
    }

    @Override
    public int getDimension() {
        return parameter.getDimension();
    }

    @Override
    public double[] getGradientLogDensity() {
        return wrtParameter.getGradientLogDensity(skygridLikelihood);
    }

    @Override
    public double[] getDiagonalHessianLogDensity() {
        return wrtParameter.getDiagonalHessianLogDensity(skygridLikelihood);
    }

    @Override
    public double[][] getHessianLogDensity() {
        throw new RuntimeException("Not yet implemented");
    }

    @Override
    public String getReport() {
        String header = skygridLikelihood + "." + wrtParameter.name + "\n";

        header += GradientWrtParameterProvider.getReportAndCheckForError(this,
                wrtParameter.getParameterLowerBound(), Double.POSITIVE_INFINITY,
                tolerance) + " \n";

        if (wrtParameter != WrtParameter.NODE_HEIGHT){
            header += HessianWrtParameterProvider.getReportAndCheckForError(this, tolerance);
        }

        return header;
    }

    @Override
    public LogColumn[] getColumns() {
        return Loggable.getColumnsFromReport(this, "GMRFGradient report");
    }

    public enum WrtParameter {

        LOG_POPULATION_SIZES("logPopulationSizes") {
            @Override
            Parameter getParameter(GMRFMultilocusSkyrideLikelihood likelihood) {
                return likelihood.getPopSizeParameter();
            }

            @Override
            double[] getGradientLogDensity(GMRFMultilocusSkyrideLikelihood likelihood) {
                return likelihood.getGradientWrtLogPopulationSize();
            }

            @Override
            double[] getDiagonalHessianLogDensity(GMRFMultilocusSkyrideLikelihood likelihood) {
                return likelihood.getDiagonalHessianWrtLogPopulationSize();
            }

            @Override
            double getParameterLowerBound() { return Double.NEGATIVE_INFINITY; }

            @Override
            public void getWarning(GMRFMultilocusSkyrideLikelihood likelihood) {

            }
        },
        PRECISION("precision") {
            @Override
            Parameter getParameter(GMRFMultilocusSkyrideLikelihood likelihood) {
                return likelihood.getPrecisionParameter();
            }

            @Override
            double[] getGradientLogDensity(GMRFMultilocusSkyrideLikelihood likelihood) {
                return likelihood.getGradientWrtPrecision();
            }

            @Override
            double[] getDiagonalHessianLogDensity(GMRFMultilocusSkyrideLikelihood likelihood) {
                return likelihood.getDiagonalHessianWrtPrecision();
            }

            @Override
            double getParameterLowerBound() { return 0.0; }

            @Override
            public void getWarning(GMRFMultilocusSkyrideLikelihood likelihood) {

            }
        },
        REGRESSION_COEFFICIENTS("regressionCoefficients") {
            @Override
            Parameter getParameter(GMRFMultilocusSkyrideLikelihood likelihood) {
                List<Parameter> allBetas = likelihood.getBetaListParameter();

                if (allBetas.size() > 1) {
                    throw new RuntimeException("This is not the correct way of handling multidimensional parameters");
                }

                return allBetas.get(0);
            }

            @Override
            double[] getGradientLogDensity(GMRFMultilocusSkyrideLikelihood likelihood) {
                return likelihood.getGradientWrtRegressionCoefficients();
            }

            @Override
            double[] getDiagonalHessianLogDensity(GMRFMultilocusSkyrideLikelihood likelihood) {
                return likelihood.getDiagonalHessianWrtRegressionCoefficients();
            }

            @Override
            double getParameterLowerBound() { return Double.NEGATIVE_INFINITY; }

            @Override
            public void getWarning(GMRFMultilocusSkyrideLikelihood likelihood) {

            }
        },
        NODE_HEIGHT("nodeHeight") {

            Parameter parameter;

            @Override
            Parameter getParameter(GMRFMultilocusSkyrideLikelihood likelihood) {
                if (parameter == null) {
                    TreeModel treeModel = (TreeModel) likelihood.getTree(0);
                    parameter = treeModel.createNodeHeightsParameter(true, true, false);
                }
                return parameter;
            }

            @Override
            double[] getGradientLogDensity(GMRFMultilocusSkyrideLikelihood likelihood) {
                return getGradientWrtNodeHeights(likelihood);
            }

            @Override
            double[] getDiagonalHessianLogDensity(GMRFMultilocusSkyrideLikelihood likelihood) {
                return new double[likelihood.getTree(0).getInternalNodeCount()];
            }

            @Override
            double getParameterLowerBound() {
                return 0.0;
            }

            public void getWarning(GMRFMultilocusSkyrideLikelihood likelihood) {
                if (likelihood.nLoci() > 1) {
                    throw new RuntimeException("Not yet implemented for multiple loci.");
                }
            }

            private double[] getGradientWrtNodeHeights(GMRFMultilocusSkyrideLikelihood likelihood) {

                likelihood.getLogLikelihood();

                Tree tree = likelihood.getTree(0);

                double[] gradient = new double[tree.getInternalNodeCount()];
                double[] currentGamma = likelihood.getPopSizeParameter().getParameterValues();

                double ploidyFactor = 1 / likelihood.getPopulationFactor(0);

                final TreeIntervals intervals = likelihood.getTreeIntervals(0);

                int[] intervalIndices = new int[tree.getInternalNodeCount()];
                int[] gridIndices = new int[tree.getInternalNodeCount()];

                getGridIndexForInternalNodes(likelihood, 0, intervalIndices, gridIndices);

                for (int i = 0; i < tree.getInternalNodeCount(); i++) {
                    NodeRef node = tree.getNode(i + tree.getExternalNodeCount());

                    final int nodeIndex = getNodeHeightParameterIndex(node, tree);

                    final int numLineage = intervals.getLineageCount(intervalIndices[i]);

                    final double currentPopSize = Math.exp(-currentGamma[gridIndices[nodeIndex]]);

                    gradient[nodeIndex] += -currentPopSize * numLineage * (numLineage - 1);

                    if (!tree.isRoot(node)) {
                        final int nextNumLineage = intervals.getLineageCount(intervalIndices[i] + 1);
                        gradient[nodeIndex] -= -currentPopSize * nextNumLineage * (nextNumLineage - 1);
                    }
                }

                final double multiplier = 0.5 * ploidyFactor;
                for (int i = 0; i < gradient.length; i++) {
                    gradient[i] *= multiplier;
                }

                return gradient;
            }

            private int getNodeHeightParameterIndex(NodeRef node, Tree tree) {
                return node.getNumber() - tree.getExternalNodeCount();
            }

            private void getGridIndexForInternalNodes(GMRFMultilocusSkyrideLikelihood likelihood, int treeIndex,
                                                       int[] intervalIndices, int[] gridIndices) {
                Tree tree = likelihood.getTree(treeIndex);
                double[] sortedValues = new double[tree.getInternalNodeCount()];
                double[] nodeHeights = parameter.getParameterValues();
                int[] nodeIndices = new int[tree.getInternalNodeCount()];
                ArrayList<ComparableDouble> sortedInternalNodes = new ArrayList<ComparableDouble>();
                for (int i = 0; i < nodeIndices.length; i++) {
                    sortedInternalNodes.add(new ComparableDouble(tree.getNodeHeight(tree.getNode(tree.getExternalNodeCount() + i))));
                }
                HeapSort.sort(sortedInternalNodes, nodeIndices);
                for (int i = 0; i < nodeIndices.length; i++) {
                    sortedValues[i] = nodeHeights[nodeIndices[i]];
                }

                int gridIndex = 0;
                double[] gridPoints = likelihood.getGridPoints();
                int intervalIndex = 0;
                final TreeIntervals intervals = likelihood.getTreeIntervals(treeIndex);
                for (int i = 0; i < tree.getInternalNodeCount(); i++) {
                    while(gridIndex < gridPoints.length && gridPoints[gridIndex] < sortedValues[i]) {
                        gridIndex++;
                    }
                    gridIndices[nodeIndices[i]] = gridIndex;

                    while(intervalIndex < intervals.getIntervalCount() - 1 && intervals.getIntervalTime(intervalIndex) < sortedValues[i]) {
                        intervalIndex++;
                    }
                    intervalIndices[nodeIndices[i]] = intervalIndex;
                }
            }

        };

        WrtParameter(String name) {
            this.name = name;
        }

        abstract Parameter getParameter(GMRFMultilocusSkyrideLikelihood likelihood);

        abstract double[] getGradientLogDensity(GMRFMultilocusSkyrideLikelihood likelihood);

        abstract double[] getDiagonalHessianLogDensity(GMRFMultilocusSkyrideLikelihood likelihood);

        abstract double getParameterLowerBound();

        public abstract void getWarning(GMRFMultilocusSkyrideLikelihood likelihood);

        private final String name;

        public static WrtParameter factory(String match) {
            for (WrtParameter type : WrtParameter.values()) {
                if (match.equalsIgnoreCase(type.name)) {
                    return type;
                }
            }
            return null;
        }
    }

    private final static Double tolerance = 1E-4;
}

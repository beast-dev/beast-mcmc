package dr.evomodel.coalescent.hmc;

import dr.evolution.coalescent.IntervalType;
import dr.evolution.coalescent.TreeIntervals;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.coalescent.GMRFMultilocusSkyrideLikelihood;
import dr.evomodel.tree.TreeModel;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.HessianWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.xml.Reportable;

import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Mandev Gill
 */
public class GMRFGradient implements GradientWrtParameterProvider, HessianWrtParameterProvider, Reportable {

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

        header += HessianWrtParameterProvider.getReportAndCheckForError(this, tolerance);

        return header;
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

                int[] gridIndices = getGridIndexForInternalNodes(likelihood, 0);

                for (int i = 0; i < intervals.getIntervalCount(); i++) {
                    if (intervals.getIntervalType(i) == IntervalType.COALESCENT) {

                        final int nodeIndex = getNodeHeightParameterIndex(intervals.getCoalescentNode(i), tree);

                        final int numLineage = intervals.getLineageCount(i);

                        gradient[nodeIndex] += -Math.exp(-currentGamma[gridIndices[nodeIndex]]) * numLineage * (numLineage - 1);

                        if (!tree.isRoot(intervals.getCoalescentNode(i))) {
                            final int nextNumLineage = intervals.getLineageCount(i + 1);
                            gradient[nodeIndex] -= -Math.exp(-currentGamma[gridIndices[nodeIndex]]) * nextNumLineage * (nextNumLineage - 1);
                        }

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

            private int[] getGridIndexForInternalNodes(GMRFMultilocusSkyrideLikelihood likelihood, int treeIndex) {
                Tree tree = likelihood.getTree(treeIndex);
                TreeIntervals intervals = likelihood.getTreeIntervals(treeIndex);

                int[] indices = new int[tree.getInternalNodeCount()];

                int gridIndex = 0;
                double[] gridPoints = likelihood.getGridPoints();
                for (int i = 0; i < intervals.getIntervalCount(); i++) {
                    if (intervals.getIntervalType(i) == IntervalType.COALESCENT) {
                        while(gridPoints[gridIndex] < intervals.getInterval(i) && gridIndex < gridPoints.length - 1) {
                            gridIndex++;
                        }
                        indices[getNodeHeightParameterIndex(intervals.getCoalescentNode(i), tree)] = gridIndex;
                    }
                }

                return indices;
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

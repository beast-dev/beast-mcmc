package dr.evomodel.coalescent.hmc;

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
        },
        NODE_HEIGHT("nodeHeight") {

            Parameter parameter;

            @Override
            Parameter getParameter(GMRFMultilocusSkyrideLikelihood likelihood) {
                getWarning(likelihood);
                if (parameter == null) {
                    TreeModel treeModel = (TreeModel) likelihood.getTree(0);
                    parameter = treeModel.createNodeHeightsParameter(true, true, false);
                }
                return parameter;
            }

            @Override
            double[] getGradientLogDensity(GMRFMultilocusSkyrideLikelihood likelihood) {
                getWarning(likelihood);
                return likelihood.getGradientWrtNodeHeights();
            }

            @Override
            double[] getDiagonalHessianLogDensity(GMRFMultilocusSkyrideLikelihood likelihood) {
                getWarning(likelihood);
                return new double[likelihood.getTree(0).getInternalNodeCount()];
            }

            @Override
            double getParameterLowerBound() {
                return 0.0;
            }

            private void getWarning(GMRFMultilocusSkyrideLikelihood likelihood) {
                if (likelihood.nLoci() > 1) {
                    throw new RuntimeException("Not yet implemented for multiple loci.");
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

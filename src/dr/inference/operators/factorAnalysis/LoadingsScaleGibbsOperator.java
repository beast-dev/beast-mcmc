package dr.inference.operators.factorAnalysis;

import dr.inference.distribution.NormalStatisticsProvider;
import dr.inference.model.*;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.RejectionOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.inferencexml.operators.factorAnalysis.LoadingsOperatorParserUtilities;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.*;
import dr.xml.*;

public class LoadingsScaleGibbsOperator extends SimpleMCMCOperator implements GibbsOperator,
        RejectionOperator.RejectionProvider, VariableListener, Reportable {

    private final Parameter sccaleComponent;
    private final MatrixParameterInterface matrixComponent;

    private final FactorAnalysisStatisticsProvider statisticsProvider;
    private final int nTraits;
    private final int nFactors;

    private final NormalStatisticsProvider prior;

    private final double[] mean;
    private final double[][] variance;

    private final Parameter[] listeningParameters;
    private boolean needToUpdateStatistics = true;

    private double[][] cholesky;

    LoadingsScaleGibbsOperator(FactorAnalysisStatisticsProvider statisticsProvider,
                               NormalStatisticsProvider prior) {

        ScaledMatrixParameter loadings = (ScaledMatrixParameter) statisticsProvider.getAdaptor().getLoadings();
        this.sccaleComponent = loadings.getScaleParameter();
        this.matrixComponent = loadings.getMatrixParameter();

        this.statisticsProvider = statisticsProvider;
        this.nTraits = matrixComponent.getRowDimension();
        this.nFactors = matrixComponent.getColumnDimension();

        this.prior = prior;

        this.mean = new double[nFactors];
        this.variance = new double[nFactors][nFactors];

        this.listeningParameters = statisticsProvider.getAdaptor().getLoadingsDependentParameter();
        for (Parameter parameter : listeningParameters) {
            parameter.addParameterListener(this);
        }
    }

    @Override
    public String getOperatorName() {
        return LOADINGS_SCALE_OPERATOR;
    }

    @Override
    public double doOperation() {
        double draw[] = getProposedUpdate();

        for (int k = 0; k < nFactors; k++) {
            sccaleComponent.setParameterValueQuietly(k, draw[k]);
        }
        sccaleComponent.fireParameterChangedEvent();

        return 0;
    }

    @Override
    public double[] getProposedUpdate() {
        if (needToUpdateStatistics) {
            statisticsProvider.getAdaptor().drawFactors();
            updateMeanAndVariance();

            try {
                cholesky = new CholeskyDecomposition(variance).getL();
            } catch (IllegalDimension illegalDimension) {
                illegalDimension.printStackTrace();
            }

            needToUpdateStatistics = false;
        }

        double[] draw = MultivariateNormalDistribution.nextMultivariateNormalCholesky(mean, cholesky);
        return draw;
    }

    @Override
    public Parameter getParameter() {
        return sccaleComponent;
    }

    private void updateMeanAndVariance() {
        double[][] precision = new double[nFactors][nFactors];

        double[][] factorInnerProduct = new double[nFactors][nFactors];

        double[] meanBuffer = new double[nFactors];

        for (int j = 0; j < nTraits; j++) {

            double factorPrecision = statisticsProvider.getAdaptor().getColumnPrecision(j);

            statisticsProvider.getFactorInnerProduct(j, nFactors, factorInnerProduct);
            statisticsProvider.getFactorTraitProduct(j, nFactors, mean);

            for (int k1 = 0; k1 < nFactors; k1++) {

                double u1 = matrixComponent.getParameterValue(j, k1);

                precision[k1][k1] += factorInnerProduct[k1][k1] * u1 * u1 * factorPrecision;

                for (int k2 = k1 + 1; k2 < nFactors; k2++) {
                    double u2 = matrixComponent.getParameterValue(j, k2);
                    precision[k1][k2] += factorInnerProduct[k1][k2] * u1 * u2 * factorPrecision;
                    precision[k2][k1] = precision[k1][k2];
                }

                statisticsProvider.getFactorTraitProduct(j, nFactors, mean);
                meanBuffer[k1] += mean[k1] * u1 * factorPrecision;
            }
        }

        for (int i = 0; i < nFactors; i++) {
            double sd = prior.getNormalSD(i);
            double precI = 1.0 / (sd * sd);
            precision[i][i] += precI;
            meanBuffer[i] += precI * prior.getNormalMean(i);
        }

        Matrix varMat = (new SymmetricMatrix(precision)).inverse(); //TODO: don't invert then do Cholesky
        for (int i = 0; i < nFactors; i++) {
            for (int j = 0; j < nFactors; j++) {
                variance[i][j] = varMat.component(i, j);
            }
        }

        for (int k1 = 0; k1 < nFactors; k1++) {
            mean[k1] = 0;
            for (int k2 = 0; k2 < nFactors; k2++) {
                mean[k1] += variance[k1][k2] * meanBuffer[k2];
            }
        }
    }

    @Override
    public void variableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        for (Parameter parameter : listeningParameters) {
            if (variable == parameter) {
                needToUpdateStatistics = true;
                break;
            }
        }
    }

    @Override
    public String getReport() {
        updateMeanAndVariance();

        StringBuilder sb = new StringBuilder();
        sb.append(getOperatorName() + "Report:\n");
        sb.append("Scale mean:\n");
        sb.append(new Vector(mean));
        sb.append("\n\n");
        sb.append("Scale covariance:\n");
        sb.append(new Matrix(variance));
        sb.append("\n\n");

        return sb.toString();
    }


    private static final String LOADINGS_SCALE_OPERATOR = "loadingsScaleGibbsOperator";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            FactorAnalysisStatisticsProvider statisticsProvider =
                    LoadingsOperatorParserUtilities.parseAdaptorAndStatistics(xo);
            NormalStatisticsProvider prior = (NormalStatisticsProvider) xo.getChild(NormalStatisticsProvider.class);

            MatrixParameterInterface loadings = statisticsProvider.getAdaptor().getLoadings();
            if (!(loadings instanceof ScaledMatrixParameter)) {
                throw new XMLParseException("The loadings matrix is of class" + loadings.getClass() + ". It must be " +
                        "of class " + ScaledMatrixParameter.class);
            }
            return new LoadingsScaleGibbsOperator(statisticsProvider, prior);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return XMLSyntaxRule.Utils.concatenate(
                    LoadingsOperatorParserUtilities.statisticsProviderRules,
                    new XMLSyntaxRule[]{
                            new ElementRule(NormalStatisticsProvider.class)
                    }
            );
        }

        @Override
        public String getParserDescription() {
            return "Gibbs operator for scale component of loadings.";
        }

        @Override
        public Class getReturnType() {
            return LoadingsScaleGibbsOperator.class;
        }

        @Override
        public String getParserName() {
            return LOADINGS_SCALE_OPERATOR;
        }
    };

}

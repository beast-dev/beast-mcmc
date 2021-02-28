package dr.inference.operators.factorAnalysis;

import dr.inference.distribution.NormalStatisticsProvider;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.inference.model.ScaledMatrixParameter;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.inferencexml.operators.factorAnalysis.LoadingsOperatorParserUtilities;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.math.matrixAlgebra.CholeskyDecomposition;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.SymmetricMatrix;
import dr.xml.*;

public class LoadingsScaleGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {

    private final Parameter sccaleComponent;
    private final MatrixParameterInterface matrixComponent;

    private final FactorAnalysisStatisticsProvider statisticsProvider;
    private final int nTraits;
    private final int nFactors;

    private final NormalStatisticsProvider prior;

    LoadingsScaleGibbsOperator(FactorAnalysisStatisticsProvider statisticsProvider,
                               NormalStatisticsProvider prior) {

        ScaledMatrixParameter loadings = (ScaledMatrixParameter) statisticsProvider.getAdaptor().getLoadings();
        this.sccaleComponent = loadings.getScaleParameter();
        this.matrixComponent = loadings.getMatrixParameter();

        this.statisticsProvider = statisticsProvider;
        this.nTraits = matrixComponent.getRowDimension();
        this.nFactors = matrixComponent.getColumnDimension();

        this.prior = prior;
    }

    @Override
    public String getOperatorName() {
        return LOADINGS_SCALE_OPERATOR;
    }

    @Override
    public double doOperation() {

        double[][] precision = new double[nFactors][nFactors];

        double[][] factorInnerProduct = new double[nFactors][nFactors];

        double[] mean = new double[nFactors];
        double[] meanBuffer = new double[nFactors];

        for (int j = 0; j < nTraits; j++) {

            double factorPrecision = statisticsProvider.getAdaptor().getColumnPrecision(j);

            statisticsProvider.getFactorInnerProduct(j, nFactors, factorInnerProduct);
            statisticsProvider.getFactorTraitProduct(j, nFactors, meanBuffer);

            for (int k1 = 0; k1 < nFactors; k1++) {

                double u1 = matrixComponent.getParameterValue(j, k1);

                precision[k1][k1] += factorInnerProduct[k1][k1] * u1 * u1 * factorPrecision;

                for (int k2 = k1 + 1; k2 < nFactors; k2++) {
                    double u2 = matrixComponent.getParameterValue(j, k2);
                    precision[k1][k2] += factorInnerProduct[k1][k2] * u1 * u2 * factorPrecision;
                    precision[k2][k1] = precision[k1][k2];
                }

                statisticsProvider.getFactorTraitProduct(j, nFactors, meanBuffer);
                mean[k1] += meanBuffer[k1] * u1 * factorPrecision;
            }
        }

        for (int i = 0; i < nFactors; i++) {
            double sd = prior.getNormalSD(i);
            double precI = 1.0 / (sd * sd);
            precision[i][i] += precI;
            mean[i] += precI * prior.getNormalMean(i);
        }

        double[][] variance = (new SymmetricMatrix(precision)).inverse().toComponents();

        for (int k1 = 0; k1 < nFactors; k1++) {
            meanBuffer[k1] = 0;
            for (int k2 = 0; k2 < nFactors; k2++) {
                meanBuffer[k1] += variance[k1][k2] * mean[k2];
            }
        }

        double[][] cholesky = null;
        try {
            cholesky = new CholeskyDecomposition(variance).getL();
        } catch (IllegalDimension illegalDimension) {
            illegalDimension.printStackTrace();
        }

        double[] draw = MultivariateNormalDistribution.nextMultivariateNormalCholesky(meanBuffer, cholesky);
        for (int k = 0; k < nFactors; k++) {
            sccaleComponent.setParameterValueQuietly(k, draw[k]);
        }
        sccaleComponent.fireParameterChangedEvent();

        return 0;
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

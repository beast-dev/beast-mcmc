package dr.inference.operators.factorAnalysis;

import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.*;
import dr.inferencexml.operators.factorAnalysis.LoadingsOperatorParserUtilities;
import dr.xml.*;


public class SampledLoadingsGradient implements GradientWrtParameterProvider {

    private final MatrixParameterInterface loadings;

    private final FactorAnalysisStatisticsProvider statisticsProvider;
    private final FactorAnalysisOperatorAdaptor adaptor;

    private final double[][] scaledFactorTraitProducts;
    private final double[][][] precisions;
    private final int nFactors;
    private final int nTraits;
    private boolean statisticsKnown = false;

    private Likelihood likelihood;

    SampledLoadingsGradient(FactorAnalysisStatisticsProvider statisticsProvider) {

        this.statisticsProvider = statisticsProvider;
        this.adaptor = statisticsProvider.getAdaptor();
        this.loadings = adaptor.getLoadings();

        this.nFactors = adaptor.getNumberOfFactors();
        this.nTraits = adaptor.getNumberOfTraits();
        this.scaledFactorTraitProducts = new double[nTraits][nFactors];
        this.precisions = new double[nTraits][nFactors][nFactors];

        this.likelihood = new CompoundLikelihood(adaptor.getLikelihoods());
    }

    @Override
    public Likelihood getLikelihood() {
        return likelihood;
    }

    @Override
    public Parameter getParameter() {
        return loadings;
    }

    @Override
    public int getDimension() {
        return loadings.getDimension();
    }

    private void updateStatistics() {
        adaptor.drawFactors();
        for (int i = 0; i < nTraits; i++) {
            statisticsProvider.getScaledFactorInnerProduct(i, nFactors, precisions[i]);
            statisticsProvider.getScaledFactorTraitProduct(i, nFactors, scaledFactorTraitProducts[i]);
        }

        statisticsKnown = true;
    }


    @Override
    public double[] getGradientLogDensity() {
        updateStatistics();

        double[] gradient = new double[getDimension()];

        for (int i = 0; i < nTraits; i++) {
            for (int j = 0; j < nFactors; j++) {
                int index = j * nTraits + i;
                gradient[index] = scaledFactorTraitProducts[i][j];
                for (int k = 0; k < nFactors; k++) {
                    gradient[index] -= precisions[i][j][k] * loadings.getParameterValue(i, k);
                }
            }
        }
        return gradient;
    }


    private static final String SAMPLED_LOADINGS_GRADIENT = "sampledLoadingsGradient";


    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            FactorAnalysisStatisticsProvider statisticsProvider =
                    LoadingsOperatorParserUtilities.parseAdaptorAndStatistics(xo);
            return new SampledLoadingsGradient(statisticsProvider);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return LoadingsOperatorParserUtilities.statisticsProviderRules;
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return SampledLoadingsGradient.class;
        }

        @Override
        public String getParserName() {
            return SAMPLED_LOADINGS_GRADIENT;
        }
    };


}

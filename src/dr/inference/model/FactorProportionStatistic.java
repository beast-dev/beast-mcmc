package dr.inference.model;

import dr.inference.operators.factorAnalysis.FactorAnalysisOperatorAdaptor;
import dr.inference.operators.factorAnalysis.FactorAnalysisStatisticsProvider;
import dr.inferencexml.operators.factorAnalysis.LoadingsOperatorParserUtilities;
import dr.xml.*;

public class FactorProportionStatistic extends Statistic.Abstract implements Reportable {

    private final FactorAnalysisOperatorAdaptor adaptor;
    private final FactorAnalysisStatisticsProvider statisticsProvider;

    private final double[] relativeFactorProportions;
    private final double[] absoluteFactorProportions;
    private double relativeMarginalProportion;
    private double factorProportion;

    FactorProportionStatistic(FactorAnalysisOperatorAdaptor adaptor) {
        this.adaptor = adaptor;
        this.statisticsProvider = new FactorAnalysisStatisticsProvider(adaptor,
                FactorAnalysisStatisticsProvider.CacheProvider.NO_CACHE);

        int k = adaptor.getNumberOfFactors();
        this.relativeFactorProportions = new double[k];
        this.absoluteFactorProportions = new double[k];
    }


    @Override
    public int getDimension() {
        int k = adaptor.getNumberOfFactors();
        return 2 * k + 2;
    }

    @Override
    public double getStatisticValue(int dim) {
        int k = adaptor.getNumberOfFactors();

        if (dim == 0) {
            computeStatistics(); //TODO: add listeners instead
        }

        if (dim == 0) {
            return factorProportion;
        } else if (dim < (k + 1)) {
            return absoluteFactorProportions[dim - 1];
        } else if (dim < (2 * k + 1)) {
            return relativeFactorProportions[dim - (k + 1)];
        } else if (dim == getDimension()) {
            return relativeMarginalProportion;
        } else {
            throw new RuntimeException("Unknown dimension: " + dim);
        }
    }

    @Override
    public String getDimensionName(int dim) {
        int k = adaptor.getNumberOfFactors();
        final String dimName;
        if (dim == 0) {
            dimName = "factorProportion";
        } else if (dim < (k + 1)) {
            dimName = "absoluteProportion." + dim;
        } else if (dim < (2 * k + 1)) {
            dimName = "relativeProportion." + (dim - k);
        } else if (dim == getDimension()) {
            dimName = "relativeMarginalProportion";
        } else {
            throw new RuntimeException("Unknown dimension: " + dim);
        }

        return getStatisticName() + "." + dimName;
    }

    private void computeStatistics() {
        int k = adaptor.getNumberOfFactors();
        int n = adaptor.getNumberOfTaxa();
        int p = adaptor.getNumberOfTraits();

        double[][] factorInnerProduct = new double[k][k];
        statisticsProvider.getFactorInnerProduct(-1, k, factorInnerProduct, false);

        double[][] loadingsInnerProduct = new double[k][k];
        statisticsProvider.getLoadingsInnerProduct(loadingsInnerProduct);

        double[] factorMeans = new double[k];
        statisticsProvider.getFactorMeans(factorMeans);

        double[][] factorComponent = new double[k][k];
        double factorSum = 0;
        double marginalSum = 0;
        for (int f1 = 0; f1 < k; f1++) {
            for (int f2 = f1; f2 < k; f2++) {
                double value = loadingsInnerProduct[f1][f2] *
                        (factorInnerProduct[f1][f2] - n * factorMeans[f1] * factorMeans[f2]);
                factorComponent[f1][f2] = value;
                factorSum += value;
                if (f1 != f2) {
                    factorComponent[f2][f1] = value;
                    factorSum += value;
                } else {
                    marginalSum += value;
                }
            }
        }

        double errorSum = 0;
        for (int j = 0; j < p; j++) {
            errorSum += 1 / adaptor.getColumnPrecision(j);
        }

        errorSum *= (n - 1);

        double totalSum = factorSum + errorSum;

        factorProportion = factorSum / totalSum;
        relativeMarginalProportion = marginalSum / factorSum;
        for (int i = 0; i < k; i++) {
            relativeFactorProportions[i] = factorComponent[i][i] / factorSum;
            absoluteFactorProportions[i] = factorComponent[i][i] / totalSum;
        }
    }

    private static final String FACTOR_PROPORTION = "factorProportionStatistic";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            FactorAnalysisOperatorAdaptor adaptor =
                    LoadingsOperatorParserUtilities.parseFactorAnalsysisOperatorAdaptor(xo);
            return new FactorProportionStatistic(adaptor);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return LoadingsOperatorParserUtilities.adaptorRules;
        }

        @Override
        public String getParserDescription() {
            return "Absolute and relative contribution of each factor to the variance in the data";
        }

        @Override
        public Class getReturnType() {
            return FactorProportionStatistic.class;
        }

        @Override
        public String getParserName() {
            return FACTOR_PROPORTION;
        }
    };

    @Override
    public String getReport() {
        return null;
    }
}

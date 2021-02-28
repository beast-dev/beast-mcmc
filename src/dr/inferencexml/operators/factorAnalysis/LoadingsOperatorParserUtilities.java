package dr.inferencexml.operators.factorAnalysis;

import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.IntegratedFactorAnalysisLikelihood;
import dr.inference.model.LatentFactorModel;
import dr.inference.operators.factorAnalysis.FactorAnalysisOperatorAdaptor;
import dr.inference.operators.factorAnalysis.FactorAnalysisStatisticsProvider;
import dr.xml.*;


public class LoadingsOperatorParserUtilities {

    private final static String USE_CACHE = "cacheInnerProducts";

    static FactorAnalysisStatisticsProvider parseAdaptorAndStatistics(XMLObject xo) throws XMLParseException {


        final FactorAnalysisStatisticsProvider.CacheProvider cacheProvider;
        boolean useCache = xo.getAttribute(USE_CACHE, false);
        if (useCache) {
            cacheProvider = FactorAnalysisStatisticsProvider.CacheProvider.USE_CACHE;
        } else {
            cacheProvider = FactorAnalysisStatisticsProvider.CacheProvider.NO_CACHE;
        }

        FactorAnalysisOperatorAdaptor adaptor = parseFactorAnalsysisOperatorAdaptor(xo);

        return new FactorAnalysisStatisticsProvider(adaptor, cacheProvider);
    }


    static FactorAnalysisOperatorAdaptor parseFactorAnalsysisOperatorAdaptor(XMLObject xo) {
        LatentFactorModel factorModel = (LatentFactorModel) xo.getChild(LatentFactorModel.class);

        if (factorModel == null) {
            IntegratedFactorAnalysisLikelihood integratedModel =
                    (IntegratedFactorAnalysisLikelihood) xo.getChild(IntegratedFactorAnalysisLikelihood.class);
            TreeDataLikelihood treeLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);

            return new FactorAnalysisOperatorAdaptor.IntegratedFactors(integratedModel, treeLikelihood);

        } else {
            return new FactorAnalysisOperatorAdaptor.SampledFactors(factorModel);
        }
    }

    public static final XMLSyntaxRule[] adaptorRules = new XMLSyntaxRule[]{
            new XORRule(
                    new ElementRule(LatentFactorModel.class),
                    new AndRule(
                            new ElementRule(IntegratedFactorAnalysisLikelihood.class),
                            new ElementRule(TreeDataLikelihood.class)
                    )
            )
    };


    public static final XMLSyntaxRule[] statisticsProviderRules =
            XMLSyntaxRule.Utils.concatenate(adaptorRules,
                    new XMLSyntaxRule[]{AttributeRule.newBooleanRule(USE_CACHE, true)});
}

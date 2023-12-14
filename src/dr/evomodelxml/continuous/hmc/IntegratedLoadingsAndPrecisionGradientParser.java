package dr.evomodelxml.continuous.hmc;

import dr.evomodel.continuous.hmc.IntegratedLoadingsAndPrecisionGradient;
import dr.evomodel.continuous.hmc.IntegratedLoadingsGradient;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.evomodel.treedatalikelihood.continuous.IntegratedFactorAnalysisLikelihood;
import dr.inference.model.CompoundParameter;
import dr.util.TaskPool;
import dr.xml.ElementRule;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class IntegratedLoadingsAndPrecisionGradientParser extends IntegratedLoadingsGradientParser {

    public static final String PARSER_NAME = "integratedFactorAnalysisLoadingsAndPrecisionGradient";

    protected IntegratedLoadingsGradient factory(TreeDataLikelihood treeDataLikelihood,
                                                 ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                 IntegratedFactorAnalysisLikelihood factorAnalysisLikelihood,
                                                 ContinuousTraitPartialsProvider jointPartialsProvider,
                                                 TaskPool taskPool,
                                                 IntegratedLoadingsGradient.ThreadUseProvider threadUseProvider,
                                                 IntegratedLoadingsGradient.RemainderCompProvider remainderCompProvider,
                                                 CompoundParameter parameter)
            throws XMLParseException {

        return new IntegratedLoadingsAndPrecisionGradient(
                parameter,
                treeDataLikelihood,
                likelihoodDelegate,
                factorAnalysisLikelihood,
                jointPartialsProvider,
                taskPool,
                threadUseProvider,
                remainderCompProvider);

    }


    @Override
    public String getParserDescription() {
        return "Generates a gradient provider for the loadings matrix & precision when factors are integrated out";
    }

    @Override
    public Class getReturnType() {
        return IntegratedLoadingsAndPrecisionGradient.class;
    }

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        XMLSyntaxRule[] newRules = new XMLSyntaxRule[rules.length + 1];
        newRules[0] = new ElementRule(CompoundParameter.class);
        System.arraycopy(rules, 0, newRules, 1, rules.length);
        return newRules;
    }
}

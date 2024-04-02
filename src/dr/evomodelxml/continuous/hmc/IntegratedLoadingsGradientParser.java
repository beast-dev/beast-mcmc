package dr.evomodelxml.continuous.hmc;

import dr.evomodel.continuous.hmc.IntegratedLoadingsGradient;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.evomodel.treedatalikelihood.continuous.IntegratedFactorAnalysisLikelihood;
import dr.evomodel.treedatalikelihood.continuous.JointPartialsProvider;
import dr.inference.model.CompoundParameter;
import dr.util.TaskPool;
import dr.xml.*;

public class IntegratedLoadingsGradientParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "integratedFactorAnalysisLoadingsGradient";
    private static final String THREAD_TYPE = "threadType";
    private static final String PARALLEL = "parallel";
    private static final String SERIAL = "serial";
    private static final String REMAINDER_COMPUTATION = "remainderComputation";
    private static final String FULL = "full";
    private static final String SKIP = "skip";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        TreeDataLikelihood treeDataLikelihood = (TreeDataLikelihood)
                xo.getChild(TreeDataLikelihood.class);

        IntegratedFactorAnalysisLikelihood factorAnalysis = (IntegratedFactorAnalysisLikelihood)
                xo.getChild(IntegratedFactorAnalysisLikelihood.class);

        DataLikelihoodDelegate likelihoodDelegate = treeDataLikelihood.getDataLikelihoodDelegate();

        if (!(likelihoodDelegate instanceof ContinuousDataLikelihoodDelegate)) {
            throw new XMLParseException("TODO");
        }

        ContinuousDataLikelihoodDelegate continuousDataLikelihoodDelegate =
                (ContinuousDataLikelihoodDelegate) likelihoodDelegate;

        TaskPool taskPool = (TaskPool) xo.getChild(TaskPool.class);

        String threadType = xo.getAttribute(THREAD_TYPE, PARALLEL);
        IntegratedLoadingsGradient.ThreadUseProvider threadProvider;

        if (threadType.equalsIgnoreCase(PARALLEL)) {
            threadProvider = IntegratedLoadingsGradient.ThreadUseProvider.PARALLEL;
        } else if (threadType.equalsIgnoreCase(SERIAL)) {
            threadProvider = IntegratedLoadingsGradient.ThreadUseProvider.SERIAL;
        } else {
            throw new XMLParseException("The attribute " + THREAD_TYPE + " must have values \"" + PARALLEL +
                    "\" or \"" + SERIAL + "\".");
        }

        String remComp = xo.getAttribute(REMAINDER_COMPUTATION, SKIP);
        IntegratedLoadingsGradient.RemainderCompProvider remainderCompProvider;

        if (remComp.equalsIgnoreCase(SKIP)) {
            remainderCompProvider = IntegratedLoadingsGradient.RemainderCompProvider.SKIP;
        } else if (remComp.equalsIgnoreCase(FULL)) {
            remainderCompProvider = IntegratedLoadingsGradient.RemainderCompProvider.FULL;
        } else {
            throw new XMLParseException("The attribute " + REMAINDER_COMPUTATION + " must have values \"" + SKIP +
                    "\" or \"" + FULL + "\".");
        }

        if (taskPool != null && threadType != PARALLEL) {
            throw new XMLParseException("Cannot simultaneously provide " + TaskPoolParser.TASK_PARSER_NAME +
                    " and " + THREAD_TYPE + "=\"" + threadType + "\". Please either change to " + THREAD_TYPE +
                    "=\"" + PARALLEL + "\" or remove the " + TaskPoolParser.TASK_PARSER_NAME + " element.");
        }

        ContinuousTraitPartialsProvider partialsProvider = (JointPartialsProvider) xo.getChild(JointPartialsProvider.class);
        if (partialsProvider == null) partialsProvider = factorAnalysis;

        // for IntegratedLoadingsAndPrecisionGradient
        CompoundParameter parameter = (CompoundParameter) xo.getChild(CompoundParameter.class);
        if (parameter != null) {
            if (parameter.getParameterCount() != 2) {
                throw new XMLParseException("The parameter must have two elements, " +
                        "the first being the loadings matrix and the second being the precision matrix.");
            }
            if (parameter.getParameter(0) != factorAnalysis.getLoadings()) {
                throw new XMLParseException("The first element of the parameter must be the loadings matrix.");
            }
            if (parameter.getParameter(1) != factorAnalysis.getPrecision()) {
                throw new XMLParseException("The second element of the parameter must be the precision matrix.");
            }
        }


        // TODO Check dimensions, parameters, etc.

        return factory(
                treeDataLikelihood,
                continuousDataLikelihoodDelegate,
                factorAnalysis,
                partialsProvider,
                taskPool,
                threadProvider,
                remainderCompProvider,
                parameter);

    }

    protected IntegratedLoadingsGradient factory(TreeDataLikelihood treeDataLikelihood,
                                                 ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                 IntegratedFactorAnalysisLikelihood factorAnalysisLikelihood,
                                                 ContinuousTraitPartialsProvider jointPartialsProvider,
                                                 TaskPool taskPool,
                                                 IntegratedLoadingsGradient.ThreadUseProvider threadUseProvider,
                                                 IntegratedLoadingsGradient.RemainderCompProvider remainderCompProvider,
                                                 CompoundParameter parameter)
            throws XMLParseException {

        return new IntegratedLoadingsGradient(
                treeDataLikelihood,
                likelihoodDelegate,
                factorAnalysisLikelihood,
                jointPartialsProvider,
                taskPool,
                threadUseProvider,
                remainderCompProvider);

    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return "Generates a gradient provider for the loadings matrix when factors are integrated out";
    }

    @Override
    public Class getReturnType() {
        return IntegratedLoadingsGradient.class;
    }

    @Override
    public String getParserName() {
        return PARSER_NAME;
    }

    protected final XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(IntegratedFactorAnalysisLikelihood.class),
            new ElementRule(TreeDataLikelihood.class),
            new ElementRule(TaskPool.class, true),
            AttributeRule.newStringRule(THREAD_TYPE, true),
            AttributeRule.newStringRule(REMAINDER_COMPUTATION, true),
            new ElementRule(JointPartialsProvider.class, true)

    };
}

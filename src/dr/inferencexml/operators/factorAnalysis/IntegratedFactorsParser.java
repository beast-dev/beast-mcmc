package dr.inferencexml.operators.factorAnalysis;

import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.IntegratedFactorAnalysisLikelihood;
import dr.inference.operators.factorAnalysis.FactorAnalysisOperatorAdaptor;
import dr.xml.*;

/**
 * @author Gabriel Hassler
 */

public class IntegratedFactorsParser extends AbstractXMLObjectParser {

    private static String INTEGRATED_FACTORS = "integratedFactors";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        TreeDataLikelihood treeDataLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
        IntegratedFactorAnalysisLikelihood factorLikelihood =
                (IntegratedFactorAnalysisLikelihood) xo.getChild(IntegratedFactorAnalysisLikelihood.class);
        return new FactorAnalysisOperatorAdaptor.IntegratedFactors(factorLikelihood, treeDataLikelihood);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(TreeDataLikelihood.class),
                new ElementRule(IntegratedFactorAnalysisLikelihood.class)
        };
    }

    @Override
    public String getParserDescription() {
        return "Class used internally for drawing (latent) factors.";
    }

    @Override
    public Class getReturnType() {
        return FactorAnalysisOperatorAdaptor.IntegratedFactors.class;
    }

    @Override
    public String getParserName() {
        return INTEGRATED_FACTORS;
    }
}

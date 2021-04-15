package dr.evomodelxml.continuous.hmc;

import dr.evomodel.continuous.hmc.IntegratedLoadingsGradient;
import dr.evomodel.continuous.hmc.NormalizedIntegratedLoadingsGradient;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.IntegratedFactorAnalysisLikelihood;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.ScaledMatrixParameter;
import dr.util.TaskPool;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class NormalizedIntegratedLoadingsGradientParser extends IntegratedLoadingsGradientParser {

    private static final String NORMALIZED_GRADIENT = "normalizedIntegratedLoadingsGradient";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        return super.parseXMLObject(xo);
    }

    @Override
    protected NormalizedIntegratedLoadingsGradient factory(TreeDataLikelihood treeDataLikelihood,
                                                           ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                           IntegratedFactorAnalysisLikelihood factorAnalysisLikelihood,
                                                           TaskPool taskPool,
                                                           IntegratedLoadingsGradient.ThreadUseProvider threadUseProvider,
                                                           IntegratedLoadingsGradient.RemainderCompProvider remainderCompProvider)
            throws XMLParseException {

        MatrixParameterInterface parameter = factorAnalysisLikelihood.getLoadings();
        if (!(parameter instanceof ScaledMatrixParameter)) {
            throw new XMLParseException("The " + IntegratedFactorAnalysisLikelihood.INTEGRATED_FACTOR_Model +
                    " element must be constructed with a " + ScaledMatrixParameter.SCALED_MATRIX + ".");
        }

        return new NormalizedIntegratedLoadingsGradient(
                treeDataLikelihood,
                likelihoodDelegate,
                factorAnalysisLikelihood,
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
        return "Gradient for normalized component of loadings matrix.";
    }

    @Override
    public Class getReturnType() {
        return NormalizedIntegratedLoadingsGradient.class;
    }

    @Override
    public String getParserName() {
        return NORMALIZED_GRADIENT;
    }
}

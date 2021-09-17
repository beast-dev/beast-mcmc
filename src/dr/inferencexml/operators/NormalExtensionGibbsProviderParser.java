package dr.inferencexml.operators;

import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.preorder.ModelExtensionProvider;
import dr.inference.operators.repeatedMeasures.GammaGibbsProvider;
import dr.xml.*;

/**
 * @author Gabriel Hassler
 * @author Marc A. Suchard
 */

public class NormalExtensionGibbsProviderParser extends AbstractXMLObjectParser {

    private static final String NORMAL_EXTENSION = "normalExtension";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        ModelExtensionProvider.NormalExtensionProvider dataModel = (ModelExtensionProvider.NormalExtensionProvider)
                xo.getChild(ModelExtensionProvider.NormalExtensionProvider.class);

        TreeDataLikelihood likelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);

        return new GammaGibbsProvider.NormalExtensionGibbsProvider(dataModel, likelihood);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(ModelExtensionProvider.NormalExtensionProvider.class),
                new ElementRule(TreeDataLikelihood.class)
        };
    }

    @Override
    public String getParserDescription() {
        return "Provides sufficient statistics for normal precision with gamma prior for extended tree model.";
    }

    @Override
    public Class getReturnType() {
        return GammaGibbsProvider.NormalExtensionGibbsProvider.class;
    }

    @Override
    public String getParserName() {
        return NORMAL_EXTENSION;
    }
}

package dr.evomodelxml.continuous;

import dr.evomodel.continuous.FullyConjugateMultivariateTraitLikelihood;
import dr.evomodel.continuous.GaussianProcessFromTree;
import dr.xml.*;

/**
 * Created by max on 12/9/14.
 */
public class GaussianProcessFromTreeParser extends AbstractXMLObjectParser {
    public static final String GAUSSIAN_PROCESS_FROM_TREE="gaussianProcessFromTree";
    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(FullyConjugateMultivariateTraitLikelihood.class),
    };

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        FullyConjugateMultivariateTraitLikelihood traitModel=(FullyConjugateMultivariateTraitLikelihood)xo.getChild(FullyConjugateMultivariateTraitLikelihood.class);
        return new GaussianProcessFromTree(traitModel);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    @Override
    public String getParserDescription() {
        return "Returns a random draw of traits given a trait model and a prior";
    }

    @Override
    public Class getReturnType() {
        return GaussianProcessFromTree.class;
    }

    @Override
    public String getParserName() {
        return GAUSSIAN_PROCESS_FROM_TREE;
    }
}

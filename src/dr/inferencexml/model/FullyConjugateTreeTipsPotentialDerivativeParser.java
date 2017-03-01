package dr.inferencexml.model;

import dr.evomodel.continuous.FullyConjugateMultivariateTraitLikelihood;
import dr.inference.model.FullyConjugateTreeTipsPotentialDerivative;
import dr.xml.*;

/**
 * @author Max Tolkoff
 */
public class FullyConjugateTreeTipsPotentialDerivativeParser extends AbstractXMLObjectParser {
    public final static String FULLY_CONJUGATE_TREE_TIPS_POTENTIAL_DERIVATIVE = "fullyConjugateTreeTipsPotential";


    @Override
    public String getParserName() {
        return FULLY_CONJUGATE_TREE_TIPS_POTENTIAL_DERIVATIVE;
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        FullyConjugateMultivariateTraitLikelihood treeLikelihood = (FullyConjugateMultivariateTraitLikelihood) xo.getChild(FullyConjugateMultivariateTraitLikelihood.class);

        return new FullyConjugateTreeTipsPotentialDerivative(treeLikelihood);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(FullyConjugateMultivariateTraitLikelihood.class),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return FullyConjugateTreeTipsPotentialDerivative.class;
    }
}

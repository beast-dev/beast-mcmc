package dr.inferencexml.model;

import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.CrossValidationProvider;
import dr.inference.model.Parameter;
import dr.inference.model.TraitValidationProvider;
import dr.xml.*;

/**
 * @author Gabriel Hassler
 */

@Deprecated
public class TraitValidationParser extends AbstractXMLObjectParser {

    private static final String TRAIT_VALIDATION = "traitValidation";

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        TraitValidationProvider provider = TraitValidationProviderParser.parseTraitValidationProvider(xo);

        boolean logSum = xo.getAttribute(CrossValidatorParser.LOG_SUM, false);

        if (logSum) return new CrossValidationProvider.CrossValidatorSum(provider);
        return new CrossValidationProvider.CrossValidator(provider);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                new ElementRule(TreeDataLikelihood.class),
                AttributeRule.newStringRule(TreeTraitParserUtilities.TRAIT_NAME),
                AttributeRule.newStringRule(TraitValidationProviderParser.INFERRED_NAME),
                new ElementRule(TreeTraitParserUtilities.TRAIT_PARAMETER, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                }),
                new ElementRule(TraitValidationProviderParser.MASK, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                }, true),
                AttributeRule.newBooleanRule(CrossValidatorParser.LOG_SUM)
        };
    }

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return CrossValidationProvider.CrossValidator.class;
    }

    @Override
    public String getParserName() {
        return TRAIT_VALIDATION;
    }
}

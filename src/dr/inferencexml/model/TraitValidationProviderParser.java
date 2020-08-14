package dr.inferencexml.model;

import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.CrossValidationProvider;
import dr.inference.model.Parameter;
import dr.inference.model.TraitValidationProvider;
import dr.xml.*;

/**
 * @author Gabriel Hassler
 */

public class TraitValidationProviderParser extends AbstractXMLObjectParser {

    public final static String TRAIT_VALIDATION_PROVIDER = "traitValidationProvider";
    final static String MASK = "mask";
    final static String INFERRED_NAME = "inferredTrait";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        return parseTraitValidationProvider(xo);

    }

    public static TraitValidationProvider parseTraitValidationProvider(XMLObject xo) throws XMLParseException {
        String trueValuesName = xo.getStringAttribute(TreeTraitParserUtilities.TRAIT_NAME);
        String inferredValuesName = xo.getStringAttribute(INFERRED_NAME);

        TreeDataLikelihood treeLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);

        ContinuousDataLikelihoodDelegate delegate =
                (ContinuousDataLikelihoodDelegate) treeLikelihood.getDataLikelihoodDelegate();

        ContinuousTraitPartialsProvider dataModel = delegate.getDataModel();
        Tree treeModel = treeLikelihood.getTree();


        TreeTraitParserUtilities utilities = new TreeTraitParserUtilities();


        TreeTraitParserUtilities.TraitsAndMissingIndices returnValue =
                utilities.parseTraitsFromTaxonAttributes(xo, trueValuesName,
                        treeModel, true);

        Parameter trueParameter = returnValue.traitParameter;
        boolean[] trueMissing = returnValue.getMissingIndicators();
        Parameter missingParameter = null;
        if (xo.hasChildNamed(MASK)) {
            missingParameter = (Parameter) xo.getElementFirstChild(MASK);
        }


        String id = xo.getId();


        TraitValidationProvider provider = new TraitValidationProvider(trueParameter, dataModel, treeModel, id,
                missingParameter, treeLikelihood, inferredValuesName, trueMissing);

        return provider;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{

                new ElementRule(TreeDataLikelihood.class),
                AttributeRule.newStringRule(TreeTraitParserUtilities.TRAIT_NAME),
                AttributeRule.newStringRule(INFERRED_NAME),
                new ElementRule(TreeTraitParserUtilities.TRAIT_PARAMETER, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                }),
                new ElementRule(MASK, new XMLSyntaxRule[]{
                        new ElementRule(Parameter.class)
                }, true)
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
        return TRAIT_VALIDATION_PROVIDER;
    }


}
package dr.util;

import dr.evolution.tree.MutableTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitDataModel;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

public class TraitValidationProvider implements CrossValidationProvider {

    final Parameter trueTraits;
    final Parameter inferredTraits;

    TraitValidationProvider(Parameter trueTraits, Parameter inferredTraits) {

        this.trueTraits = trueTraits;
        this.inferredTraits = inferredTraits;

    }

    @Override
    public Parameter getTrueParameter() {
        return trueTraits;
    }

    @Override
    public Parameter getInferredParameter() {
        return inferredTraits;
    }

    @Override
    public int[] getRelevantDimensions() {
        return new int[0];
    }

    @Override
    public String getName(int dim) {
        return null;
    }

    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        final static String PARSER_NAME = "traitValidation";


        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
//            Parameter trueParameter = (Parameter) xo.getChild(TRUE_PARAMETER).getChild(Parameter.class);
            String trueValuesName = xo.getStringAttribute(TreeTraitParserUtilities.TRAIT_NAME);
            TreeDataLikelihood treeLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);

            ContinuousDataLikelihoodDelegate delegate =
                    (ContinuousDataLikelihoodDelegate) treeLikelihood.getDataLikelihoodDelegate();

            ContinuousTraitPartialsProvider dataModel = delegate.getDataModel();
            TreeModel treeModel = (TreeModel) treeLikelihood.getTree();


            TreeTraitParserUtilities utilities = new TreeTraitParserUtilities();


            TreeTraitParserUtilities.TraitsAndMissingIndices returnValue =
                    utilities.parseTraitsFromTaxonAttributes(xo, trueValuesName,
                            treeModel, true);

            Parameter trueParameter = returnValue.traitParameter;

            Parameter inferredParameter = dataModel.getParameter();

            TraitValidationProvider provider = new TraitValidationProvider(trueParameter, inferredParameter);

            CrossValidator crossValidator = new CrossValidator(provider);

            return crossValidator;

        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
//                    new ElementRule(TRUE_PARAMETER, new XMLSyntaxRule[]{
//                            new ElementRule(Parameter.class)
//                    }),
                    new ElementRule(TreeDataLikelihood.class),
                    AttributeRule.newStringRule(TreeTraitParserUtilities.TRAIT_NAME),
                    new ElementRule(TreeTraitParserUtilities.TRAIT_PARAMETER, new XMLSyntaxRule[]{
                            new ElementRule(Parameter.class)
                    })
            };
        }

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return CrossValidator.class;
        }

        @Override
        public String getParserName() {
            return PARSER_NAME;
        }


    };
}

package dr.inference.model;

import dr.evolution.tree.TreeTrait;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.xml.*;

import static dr.evomodel.treedatalikelihood.preorder.AbstractRealizedContinuousTraitDelegate.REALIZED_TIP_TRAIT;


import java.util.List;

public class TraitValidationProvider implements CrossValidationProvider {

    final Parameter trueTraits;
    final Parameter inferredTraits;
    final int[] missingInds;
    final String[] dimNames;
    final int dimTrait;
    final ContinuousTraitPartialsProvider dataModel;
    final TreeDataLikelihood treeLikehoood;
    final TreeTrait treeTrait;

    TraitValidationProvider(Parameter trueTraits,
                            ContinuousTraitPartialsProvider dataModel,
                            TreeModel treeModel,
                            String id,
                            Parameter missingParameter,
                            Boolean useTreeTraits,
                            TreeDataLikelihood treeLikelihood,
                            String inferredValuesName) {


        this.trueTraits = trueTraits;

        this.treeLikehoood = treeLikelihood;

        if (!useTreeTraits) {
            this.treeTrait = null;
            this.inferredTraits = dataModel.getParameter();

        } else {
            this.treeTrait = treeLikelihood.getTreeTrait(REALIZED_TIP_TRAIT + "." + inferredValuesName);
            this.inferredTraits = new Parameter.Default(inferredValuesName, trueTraits.getDimension());
        }

        this.dimTrait = dataModel.getTraitDimension();
        this.dataModel = dataModel;

        this.missingInds = setupMissingInds(missingParameter);

        int nMissing = missingInds.length;


        dimNames = new String[nMissing];

        setupDimNames(treeModel, id);

    }

    private int[] setupMissingInds(Parameter missingParameter) {
        int[] missingInds;
        int nMissing = 0;
        if (missingParameter == null) {
            List<Integer> missingList = dataModel.getMissingIndices();
            nMissing = missingList.size();

            missingInds = new int[nMissing];

            for (int i = 0; i < nMissing; i++) {
                missingInds[i] = missingList.get(i);
            }

        } else {

            for (int i = 0; i < missingParameter.getSize(); i++) {
                if (missingParameter.getParameterValue(i) == 1.0) {
                    nMissing += 1;
                }
            }

            missingInds = new int[nMissing];
            int counter = 0;

            for (int i = 0; i < missingParameter.getSize(); i++) {
                if (missingParameter.getParameterValue(i) == 1.0) {
                    missingInds[counter] = i;
                    counter += 1;
                }
            }

        }

        return missingInds;

    }

    private void setupDimNames(TreeModel treeModel, String id) {
        if (id == null) {
            id = PARSER.getParserName();
        }
        int dim = 0;
        for (int i : missingInds) {
            int taxonInd = i / dimTrait;
            int traitInd = i - taxonInd * dimTrait;
            String taxonName = treeModel.getTaxonId(taxonInd);
            dimNames[dim] = id + "." + taxonName + (traitInd + 1);
            dim += 1;
        }
    }

    @Override
    public Parameter getTrueParameter() {
        return trueTraits;
    }

    @Override
    public Parameter getInferredParameter() {
        if (this.treeTrait != null) {
            updateTraitsFromTree();
        }
        return inferredTraits;
    }

    private void updateTraitsFromTree() {
        double[] tipValues = (double[]) treeTrait.getTrait(treeLikehoood.getTree(), null);
        assert (tipValues.length == inferredTraits.getDimension());
        for (int i : missingInds) {
            inferredTraits.setParameterValueQuietly(i, tipValues[i]);
        }
        inferredTraits.fireParameterChangedEvent();
    }

    @Override
    public int[] getRelevantDimensions() {
        return missingInds;
    }

    @Override
    public String getName(int dim) {
        return dimNames[dim];
    }

    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        final static String PARSER_NAME = "traitValidation";
        final static String MASK = "mask";
        final static String TREE_TRAITS = "useTreeTraits";
        final static String INFERRED_NAME = "inferredTrait";


        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
//            Parameter trueParameter = (Parameter) xo.getChild(TRUE_PARAMETER).getChild(Parameter.class);
            String trueValuesName = xo.getStringAttribute(TreeTraitParserUtilities.TRAIT_NAME);
            TreeDataLikelihood treeLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);


            ContinuousDataLikelihoodDelegate delegate =
                    (ContinuousDataLikelihoodDelegate) treeLikelihood.getDataLikelihoodDelegate();


            ContinuousTraitPartialsProvider dataModel = delegate.getDataModel();
            TreeModel treeModel = (TreeModel) treeLikelihood.getTree();

            Boolean useTreeTraits = xo.getAttribute(TREE_TRAITS, false);
            String inferredValuesName = xo.getStringAttribute(INFERRED_NAME);

            if (useTreeTraits && inferredValuesName == null) {
                throw new XMLParseException("If " + TREE_TRAITS + "=\"true\", you must provide attribute " +
                        INFERRED_NAME + "=<traitName>.");
            }


            TreeTraitParserUtilities utilities = new TreeTraitParserUtilities();


            TreeTraitParserUtilities.TraitsAndMissingIndices returnValue =
                    utilities.parseTraitsFromTaxonAttributes(xo, trueValuesName,
                            treeModel, true);

            Parameter trueParameter = returnValue.traitParameter;
            Parameter missingParameter = null;
            if (xo.hasChildNamed(MASK)) {
                missingParameter = (Parameter) xo.getElementFirstChild(MASK);
            }

            if (missingParameter != null && useTreeTraits) {
                throw new XMLParseException(PARSER_NAME + " should not have both " + MASK + " element and " +
                        TREE_TRAITS + "=\"true\".");
            }

            String id = xo.getId();


            TraitValidationProvider provider = new TraitValidationProvider(trueParameter, dataModel, treeModel, id,
                    missingParameter, useTreeTraits, treeLikelihood, inferredValuesName);

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
                    AttributeRule.newBooleanRule(TREE_TRAITS),
                    AttributeRule.newStringRule(INFERRED_NAME, true),
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
            return CrossValidator.class;
        }

        @Override
        public String getParserName() {
            return PARSER_NAME;
        }


    };
}

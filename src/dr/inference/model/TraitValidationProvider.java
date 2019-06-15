package dr.inference.model;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.evomodel.treedatalikelihood.preorder.ContinuousExtensionDelegate;
import dr.evomodel.treedatalikelihood.preorder.ModelExtensionProvider;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.xml.*;

import java.util.List;

import static dr.evomodel.treedatalikelihood.preorder.AbstractRealizedContinuousTraitDelegate.REALIZED_TIP_TRAIT;

public class TraitValidationProvider implements CrossValidationProvider {

    private final Parameter trueTraits;
    private final int[] missingInds;
    private final String[] dimNames;
    private final String sumName;
    private final int dimTrait;
    private final ContinuousExtensionDelegate extensionDelegate;


    TraitValidationProvider(Parameter trueTraits,
                            ContinuousTraitPartialsProvider dataModel,
                            Tree treeModel,
                            String id,
                            Parameter missingParameter,
                            TreeDataLikelihood treeLikelihood,
                            String inferredValuesName,
                            List<Integer> trueMissingIndices) {


        this.trueTraits = trueTraits;
        this.dimTrait = dataModel.getTraitDimension();

        this.missingInds = setupMissingInds(dataModel, missingParameter, trueMissingIndices);
        int nMissing = missingInds.length;


        TreeTrait treeTrait = treeLikelihood.getTreeTrait(REALIZED_TIP_TRAIT + "." + inferredValuesName);

        if (dataModel instanceof ModelExtensionProvider) {
            this.extensionDelegate = ((ModelExtensionProvider) dataModel).getExtensionDelegate(
                    (ContinuousDataLikelihoodDelegate) treeLikelihood.getDataLikelihoodDelegate(),
                    treeTrait,
                    treeModel);
        } else { //Simply returns the tree traits
            this.extensionDelegate = new ContinuousExtensionDelegate(
                    (ContinuousDataLikelihoodDelegate) treeLikelihood.getDataLikelihoodDelegate(),
                    treeTrait,
                    treeModel);
        }

        this.dimNames = new String[nMissing];
        setupDimNames(treeModel, id);
        sumName = getSumName(id);

    }

    private int[] setupMissingInds(ContinuousTraitPartialsProvider dataModel, Parameter missingParameter,
                                   List<Integer> trueMissing) {
        int[] missingInds;
        int nMissing = 0;
        if (missingParameter == null) {
            List<Integer> missingList = dataModel.getMissingIndices();
            missingList.removeAll(trueMissing);
            nMissing = missingList.size();

            missingInds = new int[nMissing];

            for (int i = 0; i < nMissing; i++) {
                missingInds[i] = missingList.get(i);
            }

        } else {


            for (int i = 0; i < missingParameter.getSize(); i++) {
                if (missingParameter.getParameterValue(i) == 1.0 && !trueMissing.contains(i)) {
                    //TODO: search more efficiently through the `trueMissing` array
                    nMissing += 1;
                }
            }

            missingInds = new int[nMissing];
            int counter = 0;

            for (int i = 0; i < missingParameter.getSize(); i++) {
                if (missingParameter.getParameterValue(i) == 1.0 && !trueMissing.contains(i)) {
                    //TODO: (see above)
                    missingInds[counter] = i;
                    counter += 1;
                }
            }

        }

        return missingInds;

    }

    private String getId(String id) {
        if (id == null) {
            id = PARSER.getParserName();
        }
        return id;
    }

    private void setupDimNames(Tree treeModel, String id) {
        id = getId(id);
        int dim = 0;
        for (int i : missingInds) {
            int taxonInd = i / dimTrait;
            int traitInd = i - taxonInd * dimTrait;
            String taxonName = treeModel.getTaxonId(taxonInd);
            dimNames[dim] = id + "." + taxonName + (traitInd + 1);
            dim += 1;
        }
    }

    private String getSumName(String id) {
        id = getId(id);
        return id + ".TotalSum";
    }

    @Override
    public double[] getTrueValues() {
        return trueTraits.getParameterValues();
    }

    @Override
    public double[] getInferredValues() {
        return extensionDelegate.getExtendedValues();
    }

//    private void updateTraitsFromTree() {
//        double[] tipValues = (double[]) treeTrait.getTrait(treeLikehoood.getTree(), null);
//        assert (tipValues.length == inferredTraits.getDimension());
//        for (int i : missingInds) {
//            inferredTraits.setParameterValueQuietly(i, tipValues[i]);
//        }
//        inferredTraits.fireParameterChangedEvent();
//    }

    @Override
    public int[] getRelevantDimensions() {
        return missingInds;
    }

    @Override
    public String getName(int dim) {
        return dimNames[dim];
    }

    @Override
    public String getNameSum(int dim) {
        return sumName;
    }

    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        final static String PARSER_NAME = "traitValidation";
        final static String MASK = "mask";
        final static String INFERRED_NAME = "inferredTrait";
        final static String LOG_SUM = "logSum";


        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
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
            List<Integer> trueMissing = returnValue.missingIndices;
            Parameter missingParameter = null;
            if (xo.hasChildNamed(MASK)) {
                missingParameter = (Parameter) xo.getElementFirstChild(MASK);
            }


            String id = xo.getId();


            TraitValidationProvider provider = new TraitValidationProvider(trueParameter, dataModel, treeModel, id,
                    missingParameter, treeLikelihood, inferredValuesName, trueMissing);

            boolean logSum = xo.getAttribute(LOG_SUM, false);

            if (logSum) return new CrossValidatorSum(provider);
            return new CrossValidator(provider);

        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{

                    new ElementRule(TreeDataLikelihood.class),
                    AttributeRule.newStringRule(TreeTraitParserUtilities.TRAIT_NAME),
                    AttributeRule.newStringRule(INFERRED_NAME),
                    AttributeRule.newBooleanRule(LOG_SUM, true),
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

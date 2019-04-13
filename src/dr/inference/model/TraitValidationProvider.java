package dr.inference.model;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.tree.TreeTraitProvider;
import dr.evomodel.treedatalikelihood.ExtendedProcessSimulation;
import dr.evomodel.treedatalikelihood.ProcessSimulation;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.evomodel.treedatalikelihood.preorder.AbstractContinuousExtensionDelegate;
import dr.evomodel.treedatalikelihood.preorder.ModelExtensionProvider;
import dr.evomodel.treedatalikelihood.preorder.ProcessSimulationDelegate;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.xml.*;

import java.util.List;

import static dr.evomodel.treedatalikelihood.preorder.AbstractRealizedContinuousTraitDelegate.REALIZED_TIP_TRAIT;

public class TraitValidationProvider implements CrossValidationProvider {

    private final Parameter trueTraits;
    private final Parameter inferredTraits;
    private final int[] missingInds;
    private final String[] dimNames;
    private final String sumName;
    private final int dimTrait;
    private final ContinuousTraitPartialsProvider dataModel;
    private final TreeDataLikelihood treeLikehoood;
    private final TreeTrait treeTrait;

    TraitValidationProvider(Parameter trueTraits,
                            ContinuousTraitPartialsProvider dataModel,
                            Tree treeModel,
                            String id,
                            Parameter missingParameter,
                            ProcessSimulation processSimulation,
                            Boolean useTreeTraits,
                            TreeDataLikelihood treeLikelihood,
                            String inferredValuesName,
                            List<Integer> trueMissingIndices) {


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

        this.missingInds = setupMissingInds(missingParameter, trueMissingIndices);

        int nMissing = missingInds.length;


        dimNames = new String[nMissing];

        setupDimNames(treeModel, id);
        sumName = getSumName(id);

    }

    private int[] setupMissingInds(Parameter missingParameter, List<Integer> trueMissing) {
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

    @Override
    public String getNameSum(int dim) {
        return sumName;
    }

    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        final static String PARSER_NAME = "traitValidation";
        final static String MASK = "mask";
        final static String TREE_TRAITS = "useTreeTraits";
        final static String INFERRED_NAME = "inferredTrait";
        final static String LOG_SUM = "logSum";


        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
//            Parameter trueParameter = (Parameter) xo.getChild(TRUE_PARAMETER).getChild(Parameter.class);
            String trueValuesName = xo.getStringAttribute(TreeTraitParserUtilities.TRAIT_NAME);
            TreeDataLikelihood treeLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);


            ContinuousDataLikelihoodDelegate delegate =
                    (ContinuousDataLikelihoodDelegate) treeLikelihood.getDataLikelihoodDelegate();


            ContinuousTraitPartialsProvider dataModel = delegate.getDataModel();
            Tree treeModel = treeLikelihood.getTree();

            Boolean useTreeTraits = xo.getAttribute(TREE_TRAITS, false);

            String inferredValuesName = null;

            if (xo.hasAttribute(INFERRED_NAME)) {
                inferredValuesName = xo.getStringAttribute(INFERRED_NAME);
            }


            if (useTreeTraits && inferredValuesName == null) {
                throw new XMLParseException("If " + TREE_TRAITS + "=\"true\", you must provide attribute " +
                        INFERRED_NAME + "=<traitName>.");
            }


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

            if (missingParameter != null && useTreeTraits) {
                throw new XMLParseException(PARSER_NAME + " should not have both " + MASK + " element and " +
                        TREE_TRAITS + "=\"true\".");
            }
            inferredValuesName = dataModel.getParameter().getParameterName();

            ProcessSimulationDelegate.AbstractContinuousTraitDelegate simulationDelegate =
                    (ProcessSimulationDelegate.AbstractContinuousTraitDelegate) treeLikelihood.getTreeTrait(
                            REALIZED_TIP_TRAIT + "." + inferredValuesName);
            ProcessSimulation processSimulation;
            if (dataModel instanceof ModelExtensionProvider) {
                AbstractContinuousExtensionDelegate extensionDelegate =
                        ((ModelExtensionProvider) dataModel).getExtensionDelegate(simulationDelegate, inferredValuesName);
                processSimulation = new ExtendedProcessSimulation(treeLikelihood, simulationDelegate, extensionDelegate);
            } else {
                processSimulation = new ProcessSimulation(treeLikelihood, simulationDelegate);
            }


            String id = xo.getId();


            TraitValidationProvider provider = new TraitValidationProvider(trueParameter, dataModel, treeModel, id,
                    missingParameter, processSimulation, useTreeTraits, treeLikelihood, inferredValuesName, trueMissing);

            boolean logSum = xo.getAttribute(LOG_SUM, false);

            if (logSum) return new CrossValidatorSum(provider);
            return new CrossValidator(provider);

        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
//                    new ElementRule(TRUE_PARAMETER, new XMLSyntaxRule[]{
//                            new ElementRule(Parameter.class)
//                    }),
                    new ElementRule(TreeDataLikelihood.class),
                    AttributeRule.newStringRule(TreeTraitParserUtilities.TRAIT_NAME),
                    AttributeRule.newBooleanRule(TREE_TRAITS, true),
                    AttributeRule.newStringRule(INFERRED_NAME, true),
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

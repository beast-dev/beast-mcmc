package dr.evomodel.treelikelihood.utilities;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.util.Taxon;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.preorder.ContinuousExtensionDelegate;
import dr.evomodel.treedatalikelihood.preorder.ModelExtensionProvider;
import dr.evomodelxml.treedatalikelihood.TreeDataLikelihoodParser;
import dr.evomodelxml.treelikelihood.TraitLoggerParser;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.xml.*;

import static dr.evomodel.treedatalikelihood.preorder.AbstractRealizedContinuousTraitDelegate.REALIZED_TIP_TRAIT;
import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.TRAIT_NAME;


public class ModelExtensionTraitLogger implements Loggable, Reportable {

    private final ContinuousExtensionDelegate extensionDelegate;
    private final int traitDim;
    private double[] traits;


    ModelExtensionTraitLogger(ContinuousExtensionDelegate extensionDelegate, int traitDim) {

        this.extensionDelegate = extensionDelegate;
        this.traitDim = traitDim;

    }

    @Override
    public LogColumn[] getColumns() {

        Tree tree = extensionDelegate.getTree();
        String traitName = extensionDelegate.getTreeTrait().getTraitName();

        int n = tree.getExternalNodeCount();

        LogColumn[] columns = new LogColumn[n * traitDim];

        for (int i = 0; i < n; i++) {
            int finalI = i;
            NodeRef node = tree.getExternalNode(i);
            Taxon taxon = tree.getNodeTaxon(node);
            String taxonName = (taxon != null) ? taxon.getId() : null;

            for (int j = 0; j < traitDim; j++) {
                int dim = i * traitDim + j;

                int finalJ = j;
                columns[dim] = new LogColumn.Abstract(traitName + "." + taxonName + "." + (finalJ + 1)) {
                    @Override
                    protected String getFormattedValue() {
                        if (finalI == 0 && finalJ == 0) {
                            traits = extensionDelegate.getExtendedValues();
                        }
                        ;
                        return Double.toString(traits[finalI * traitDim + finalJ]);
                    }
                };
            }
        }
        ;
        return columns;
    }


    @Override
    public String getReport() {
        return null;
    }

    public static final String MODEL_EXTENSION_LOGGER = "modelExtensionTraitLogger";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            TreeDataLikelihood dataLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);

            ModelExtensionProvider extensionProvider = (ModelExtensionProvider)
                    xo.getChild(ModelExtensionProvider.class);

            if (extensionProvider == null) {
                throw new XMLParseException(TreeDataLikelihoodParser.TREE_DATA_LIKELIHOOD + " must have child of type "
                        + ContinuousExtensionDelegate.class + " to use " + MODEL_EXTENSION_LOGGER + ". If you are not " +
                        "using an extended model, use " + TraitLoggerParser.PARSER_NAME + " to log the traits.");
            }


            Tree tree = dataLikelihood.getTree();
            String traitName = xo.getStringAttribute(TRAIT_NAME);

            int traitDim = dataLikelihood.getDataLikelihoodDelegate().getTraitDim();

            TreeTrait treeTrait = dataLikelihood.getTreeTrait(REALIZED_TIP_TRAIT + "." + traitName);

            ContinuousExtensionDelegate extensionDelegate = extensionProvider.getExtensionDelegate(treeTrait, tree);

            return new ModelExtensionTraitLogger(extensionDelegate, traitDim);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(TreeDataLikelihood.class),
                AttributeRule.newStringRule(TRAIT_NAME)

        };

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return ModelExtensionTraitLogger.class;
        }

        @Override
        public String getParserName() {
            return MODEL_EXTENSION_LOGGER;
        }
    };
}

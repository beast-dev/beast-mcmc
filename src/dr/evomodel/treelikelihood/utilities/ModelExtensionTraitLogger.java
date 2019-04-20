package dr.evomodel.treelikelihood.utilities;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.util.Taxon;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.preorder.ContinuousExtensionDelegate;
import dr.evomodel.treedatalikelihood.preorder.ModelExtensionProvider;
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

        int[] inds = new int[traitDim];
        for (int i = 0; i < traitDim; i++) {
            inds[i] = i;
        }

        int n = 1000000;
        double[] mus = new double[inds.length];
        double[][] vars = new double[inds.length][inds.length];
        for (int i = 0; i < n; i++) {
            traits = extensionDelegate.getExtendedValues();
            for (int j = 0; j < inds.length; j++) {
                double val = traits[inds[j]];
                mus[j] += val;
                vars[j][j] += val * val;

                for (int k = (j + 1); k < inds.length; k++) {
                    double val2 = traits[inds[k]];
                    vars[j][k] += val * val2;
                }
            }
        }

        for (int i = 0; i < inds.length; i++) {
            mus[i] = mus[i] / n;
        }
        for (int i = 0; i < inds.length; i++) {
            vars[i][i] = vars[i][i] / n - mus[i] * mus[i];

            for (int j = (i + 1); j < inds.length; j++) {
                vars[i][j] = vars[i][j] / n - mus[i] * mus[j];
                vars[j][i] = vars[i][j];
            }
        }


        StringBuilder sb = new StringBuilder();
        sb.append(new dr.math.matrixAlgebra.Vector(mus));
        sb.append("\n\n");
        sb.append(new dr.math.matrixAlgebra.Matrix(vars));

        return sb.toString();
    }

    public static final String MODEL_EXTENSION_LOGGER = "modelExtensionTraitLogger";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            TreeDataLikelihood dataLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);

            ModelExtensionProvider extensionProvider = (ModelExtensionProvider)
                    xo.getChild(ModelExtensionProvider.class);


            Tree tree = dataLikelihood.getTree();
            String traitName = xo.getStringAttribute(TRAIT_NAME);

            int traitDim = dataLikelihood.getDataLikelihoodDelegate().getTraitDim();

            TreeTrait treeTrait = dataLikelihood.getTreeTrait(REALIZED_TIP_TRAIT + "." + traitName);

            ContinuousExtensionDelegate extensionDelegate = extensionProvider.getExtensionDelegate(
                    (ContinuousDataLikelihoodDelegate) dataLikelihood.getDataLikelihoodDelegate(), treeTrait, tree);

            return new ModelExtensionTraitLogger(extensionDelegate, traitDim);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(TreeDataLikelihood.class),
                new ElementRule(ModelExtensionProvider.class),
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

package dr.evomodel.treelikelihood.utilities;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evolution.util.Taxon;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.evomodel.treedatalikelihood.preorder.ContinuousExtensionDelegate;
import dr.evomodel.treedatalikelihood.preorder.ModelExtensionProvider;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.math.Polynomial;
import dr.xml.*;

import static dr.evomodel.treedatalikelihood.preorder.AbstractRealizedContinuousTraitDelegate.REALIZED_TIP_TRAIT;
import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.TRAIT_NAME;

/**
 * ModelExtensionTraitLogger - logs traits under extended models (i.e. repeated measures or factor models).
 *
 * @author Gabriel Hassler
 */

public class ModelExtensionTraitLogger implements Loggable, Reportable {

    private final ContinuousExtensionDelegate extensionDelegate;
    private final int traitDim;
    private final int[] logDims;
    private double[] traits;


    ModelExtensionTraitLogger(ContinuousExtensionDelegate extensionDelegate, int[] loggableDims, int traitDim) {

        this.extensionDelegate = extensionDelegate;
        this.traitDim = traitDim;
        this.logDims = loggableDims;

    }

    @Override
    public LogColumn[] getColumns() {

        Tree tree = extensionDelegate.getTree();
        String traitName = extensionDelegate.getTreeTrait().getTraitName();


        LogColumn[] columns = new LogColumn[logDims.length];

        for (int dim : logDims) {
            int taxaInd = dim / traitDim;
            int traitInd = dim - taxaInd * traitDim;

            NodeRef node = tree.getExternalNode(taxaInd);
            Taxon taxon = tree.getNodeTaxon(node);
            String taxonName = (taxon != null) ? taxon.getId() : ("taxon_" + taxaInd);

            columns[dim] = new LogColumn.Abstract(traitName + "." + taxonName + "." + (traitInd + 1)) {
                @Override
                protected String getFormattedValue() {
                    if (dim == 0) {
                        traits = extensionDelegate.getExtendedValues();
                    }
                    return Double.toString(traits[dim]);
                }
            };


        }

        return columns;
    }


    private enum IndicesProvider {
        ALL {
            @Override
            int[] getLoggableIndices(ModelExtensionProvider extensionProvider) {

                int dim = extensionProvider.getParameter().getDimension();
                int[] allDims = new int[dim];

                for (int i = 0; i < dim; i++) {
                    allDims[i] = i;
                }

                return allDims;
            }

            @Override
            String getName() {
                return "all";
            }
        },

        MISSING {
            @Override
            int[] getLoggableIndices(ModelExtensionProvider extensionProvider) {
                int[] missingInds = new int[extensionProvider.getMissingIndices().size()];

                for (int i = 0; i < missingInds.length; i++) {
                    missingInds[i] = extensionProvider.getMissingIndices().get(i);
                }

                return missingInds;
            }

            @Override
            String getName() {
                return "missing";
            }
        };

        abstract int[] getLoggableIndices(ModelExtensionProvider extensionProvider);

        abstract String getName();
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
    public static final String DIMENSIONS = "dimensions";

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

            String dims = xo.getAttribute(DIMENSIONS, IndicesProvider.ALL.getName());
            int[] logDims;

            if (dims.equalsIgnoreCase(IndicesProvider.ALL.getName())) {

                logDims = IndicesProvider.ALL.getLoggableIndices(extensionProvider);

            } else if (dims.equalsIgnoreCase(IndicesProvider.MISSING.getName())) {

                logDims = IndicesProvider.MISSING.getLoggableIndices(extensionProvider);

            } else {
                throw new XMLParseException("The attribte \"" + DIMENSIONS + "\" must have the value \"" +
                        IndicesProvider.ALL.getName() + "\" or \"" + IndicesProvider.MISSING.getName() + "\". " +
                        "You supplied the value \"" + dims + "\".");
            }

            return new ModelExtensionTraitLogger(extensionDelegate, logDims, traitDim);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(TreeDataLikelihood.class),
                new ElementRule(ModelExtensionProvider.class),
                AttributeRule.newStringRule(TRAIT_NAME),
                AttributeRule.newStringRule(DIMENSIONS, true)
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

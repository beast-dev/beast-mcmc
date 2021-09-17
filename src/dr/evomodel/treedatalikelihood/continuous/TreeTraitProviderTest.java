package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.preorder.ModelExtensionProvider;
import dr.math.matrixAlgebra.Matrix;
import dr.xml.*;

public class TreeTraitProviderTest implements Reportable {

    private final TreeDataLikelihood treeDataLikelihood;
    private final ModelExtensionProvider.NormalExtensionProvider modelExtensionProvider;
    private final Tree tree;
    private final TreeTrait tipTrait;
    private static final int REPS = 1000;


    public TreeTraitProviderTest(TreeDataLikelihood treeDataLikelihood,
                                 ModelExtensionProvider.NormalExtensionProvider modelExtension) {

        this.treeDataLikelihood = treeDataLikelihood;
        this.modelExtensionProvider = modelExtension;
        this.tree = treeDataLikelihood.getTree();
        this.tipTrait = treeDataLikelihood.getTreeTrait(modelExtension.getTipTraitName());

    }


    @Override
    public String getReport() {
        double[] meanTipTraits = (double[]) tipTrait.getTrait(tree, null);
        double[] meanTransformedTraits = modelExtensionProvider.transformTreeTraits(meanTipTraits);


        for (int i = 1; i < REPS; i++) { //start at 1 because they're already initialized with data
            treeDataLikelihood.fireModelChanged(); // force new sample
            double[] meanTipTraitsNew = (double[]) tipTrait.getTrait(tree, null);
            double[] meanTransformedTraitsNew = modelExtensionProvider.transformTreeTraits(meanTipTraitsNew);
            for (int j = 0; j < meanTipTraits.length; j++) {
                meanTipTraits[j] += meanTipTraitsNew[j];
            }
            for (int j = 0; j < meanTransformedTraits.length; j++) {
                meanTransformedTraits[j] += meanTransformedTraitsNew[j];
            }
        }

        for (int j = 0; j < meanTipTraits.length; j++) {
            meanTipTraits[j] /= REPS;
        }
        for (int j = 0; j < meanTransformedTraits.length; j++) {
            meanTransformedTraits[j] /= REPS;
        }

        int traitDim = modelExtensionProvider.getTraitDimension();
        int dataDim = modelExtensionProvider.getDataDimension();
        int taxonCount = tree.getTaxonCount();

        Matrix matrix = new Matrix(taxonCount, traitDim);
        Matrix matrixTransformed = new Matrix(taxonCount, dataDim);


        for (int taxon = 0; taxon < taxonCount; taxon++) {
            for (int trait = 0; trait < traitDim; trait++) {
                matrix.set(taxon, trait, meanTipTraits[taxon * traitDim + trait]);
            }
            for (int factor = 0; factor < dataDim; factor++) {
                matrixTransformed.set(taxon, factor, meanTransformedTraits[taxon * dataDim + factor]);
            }
        }

        StringBuilder sb = new StringBuilder("Normal extension gibbs report for trait " +
                modelExtensionProvider.getTipTraitName() + ":\n");
        sb.append("\ttaxon order:");
        for (int i = 0; i < taxonCount; i++) {
            sb.append(" " + tree.getTaxonId(i));
        }
        sb.append("\n");
        sb.append("\ttree trait values:\n");
        sb.append(matrix.toString(2));
        sb.append("\n");
        sb.append("\ttransformed trait values:\n");
        sb.append(matrixTransformed.toString(2));
        sb.append("\n\n");
        return sb.toString();
    }

    public static final AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {

        private static final String TRAIT_PROVIDER_TEST = "treeTraitReporter";


        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            TreeDataLikelihood dataLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
            ModelExtensionProvider.NormalExtensionProvider extensionProvider =
                    (ModelExtensionProvider.NormalExtensionProvider)
                            xo.getChild(ModelExtensionProvider.NormalExtensionProvider.class);
            return new TreeTraitProviderTest(dataLikelihood, extensionProvider);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(TreeDataLikelihood.class),
                    new ElementRule(ModelExtensionProvider.NormalExtensionProvider.class)
            };
        }

        @Override
        public String getParserDescription() {
            return "Calculates the average tree traits (and transformed tree traits)";
        }

        @Override
        public Class getReturnType() {
            return TreeTraitProviderTest.class;
        }

        @Override
        public String getParserName() {
            return TRAIT_PROVIDER_TEST;
        }
    };
}

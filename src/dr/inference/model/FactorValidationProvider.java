package dr.inference.model;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.IntegratedFactorAnalysisLikelihood;
import dr.xml.*;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import static dr.evomodel.treedatalikelihood.preorder.AbstractRealizedContinuousTraitDelegate.REALIZED_TIP_TRAIT;

public class FactorValidationProvider implements CrossValidationProvider {

    private final CompoundParameter data;
    private final MatrixParameterInterface loadings;
    private final int nFac;
    private final int nTrait;
    private final int nTaxa;
    private final TreeTrait treeTrait;
    private final Tree tree;
    private final String id;

    FactorValidationProvider(IntegratedFactorAnalysisLikelihood factorModel,
                             TreeDataLikelihood treeDataLikelihood,
                             String traitName,
                             String id) {
        this.data = factorModel.getParameter();
        this.loadings = factorModel.getLoadings();

        this.nFac = factorModel.getNumberOfFactors();
        this.nTrait = factorModel.getNumberOfTraits();
        this.nTaxa = factorModel.getNumberOfTaxa();
        this.treeTrait = treeDataLikelihood.getTreeTrait(REALIZED_TIP_TRAIT + "." + traitName);
        this.tree = treeDataLikelihood.getTree();
        if (id == null) {
            this.id = FACTOR_VALIDATION;
        } else {
            this.id = id;
        }

    }

    @Override
    public double[] getTrueValues() {
        double[] trueValues = new double[nTrait * nTaxa];

        for (int taxon = 0; taxon < nTaxa; taxon++) {
            int offset = nTrait * taxon;

            for (int trait = 0; trait < nTrait; trait++) {
                trueValues[trait + offset] = data.getParameter(taxon).getParameterValue(trait);
            }
        }


        return trueValues;
    }

    @Override
    public double[] getInferredValues() {
        double[] inferredValues = new double[nTrait * nTaxa];

        DenseMatrix64F loadingsMat = new DenseMatrix64F(nFac, nTrait);
        System.arraycopy(loadings.getParameterValues(), 0, loadingsMat.data, 0, loadings.getDimension());

        double[] factors = (double[]) treeTrait.getTrait(tree, null);

        DenseMatrix64F factorBuffer = new DenseMatrix64F(nFac, 1);
        DenseMatrix64F traitBuffer = new DenseMatrix64F(nTrait, 1);

        for (int i = 0; i < nTaxa; i++) {
            int factorOffset = nFac * i;
            int traitOffset = nTrait * i;

            System.arraycopy(factors, factorOffset, factorBuffer.data, 0, nFac);

            CommonOps.multTransA(loadingsMat, factorBuffer, traitBuffer);

            System.arraycopy(traitBuffer.data, 0, inferredValues, traitOffset, nTrait);


        }

        return inferredValues;
    }

    @Override
    public int[] getRelevantDimensions() {
        int dim = nTaxa * nTrait;
        int[] dims = new int[dim];
        for (int i = 0; i < dim; i++) {
            dims[i] = i;
        }
        return dims;
    }

    @Override
    public String getName(int dim) {

        String base = id;
        if (base == null) {
            base = treeTrait.getTraitName();
        }

        int taxon = dim / nTrait;
        int trait = dim - taxon * nTrait;

        return base + "." + tree.getTaxonId(taxon) + "." + trait;
    }

    @Override
    public String getNameSum(int dim) {
        return id + ".sum";
    }


    private static final String TRAIT = "traitName";
    private static final String FACTOR_VALIDATION = "factorValidation";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            IntegratedFactorAnalysisLikelihood factorModel = (IntegratedFactorAnalysisLikelihood)
                    xo.getChild(IntegratedFactorAnalysisLikelihood.class);

            TreeDataLikelihood treeLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);

            String traitName = xo.getStringAttribute(TRAIT);
            String id = null;
            if (xo.hasId()) {
                id = xo.getId();
            }

            return new FactorValidationProvider(factorModel, treeLikelihood, traitName, id);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(IntegratedFactorAnalysisLikelihood.class),
                    new ElementRule(TreeDataLikelihood.class),
                    AttributeRule.newStringRule(TRAIT)
            };
        }

        @Override
        public String getParserDescription() {
            return "Cross-validation between the latent factors and the observed data";
        }

        @Override
        public Class getReturnType() {
            return FactorValidationProvider.class;
        }

        @Override
        public String getParserName() {
            return FACTOR_VALIDATION;
        }
    };

}

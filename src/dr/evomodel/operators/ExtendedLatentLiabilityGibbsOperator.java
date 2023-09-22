package dr.evomodel.operators;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.continuous.OrderedLatentLiabilityLikelihood;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitPartialsProvider;
import dr.evomodel.treedatalikelihood.preorder.ContinuousExtensionDelegate;
import dr.evomodel.treedatalikelihood.preorder.ModelExtensionProvider;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.CompoundParameter;
import dr.inference.model.Parameter;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;
import dr.xml.*;
import org.ejml.data.DenseMatrix64F;

import java.util.ArrayList;
import java.util.Arrays;

import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.*;

public class ExtendedLatentLiabilityGibbsOperator extends SimpleMCMCOperator implements GibbsOperator, Reportable {

    private final ContinuousExtensionDelegate extensionDelegate;
    private final OrderedLatentLiabilityLikelihood latentLiabilityLikelihood;
    private final ModelExtensionProvider.NormalExtensionProvider dataModel;

    private static final int MAX_REJECT = 1000;

    ExtendedLatentLiabilityGibbsOperator(
            ModelExtensionProvider.NormalExtensionProvider dataModel,
            ContinuousExtensionDelegate extensionDelegate,
            OrderedLatentLiabilityLikelihood latentLiabilityLikelihood) {

        this.extensionDelegate = extensionDelegate;
        this.latentLiabilityLikelihood = latentLiabilityLikelihood;
        this.dataModel = dataModel;

        if (!dataModel.diagonalVariance() && dataModel.getDataDimension() > 1) {
            throw new RuntimeException(EXTENDED_LATENT_GIBBS +
                    " is only valid for extended models with diagonal variance.");
        }

    }

    @Override
    public String getOperatorName() {
        return EXTENDED_LATENT_GIBBS;
    }

    @Override
    public double doOperation() {
        int nTaxa = extensionDelegate.getTree().getExternalNodeCount();

        double[] means = extensionDelegate.getTransformedTraits();
        DenseMatrix64F variance = dataModel.getExtensionVariance();

        ArrayList<Integer> constrainedTraits = latentLiabilityLikelihood.getConstrainedTraits();

        double[] stdev = new double[constrainedTraits.size()];
        for (int i = 0; i < constrainedTraits.size(); i++) {
            int dim = constrainedTraits.get(i);
            stdev[i] = Math.sqrt(variance.get(dim, dim));
        }

        for (int taxon = 0; taxon < nTaxa; taxon++) {
            int traitDim = 0;
            int offset = dataModel.getDataDimension() * taxon;
            for (int trait : constrainedTraits) {
                double draw = Double.NaN;
                int reps = 0;
                int dim = offset + trait;
                double mean = means[dim];
                boolean isValid = false;
                while (!isValid && reps < MAX_REJECT) {
                    draw = MathUtils.nextGaussian() * stdev[traitDim] + mean;
                    isValid = latentLiabilityLikelihood.validTraitForTip(draw, taxon, trait);

                    reps++;
                }

                if (isValid) {
                    dataModel.getParameter().getParameter(taxon).setParameterValueQuietly(trait, draw);
                }

                traitDim++;
            }
        }

        dataModel.getParameter().fireParameterChangedEvent();
        return 0;
    }

    private static final String EXTENDED_LATENT_GIBBS = "extendedLatentLiabilityGibbsOperator";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            TreeDataLikelihood treeDataLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);

            ContinuousDataLikelihoodDelegate delegate =
                    (ContinuousDataLikelihoodDelegate) treeDataLikelihood.getDataLikelihoodDelegate();


            String traitName = xo.getAttribute(TRAIT_NAME, getTipTraitNameFromDataLikelihood(treeDataLikelihood));

            TreeTrait treeTrait = treeDataLikelihood.getTreeTrait(traitName);


            ContinuousTraitPartialsProvider superDataModel = delegate.getDataModel();
            ModelExtensionProvider.NormalExtensionProvider dataModel =
                    (ModelExtensionProvider.NormalExtensionProvider) superDataModel.getProviderForTrait(traitName);

            Tree tree = treeDataLikelihood.getTree();

            ContinuousExtensionDelegate extensionDelegate = dataModel.getExtensionDelegate(delegate, treeTrait, tree);

            OrderedLatentLiabilityLikelihood latentLiabilityLikelihood =
                    (OrderedLatentLiabilityLikelihood) xo.getChild(OrderedLatentLiabilityLikelihood.class);

            return new ExtendedLatentLiabilityGibbsOperator(dataModel, extensionDelegate, latentLiabilityLikelihood);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(TreeDataLikelihood.class),
                    new ElementRule(OrderedLatentLiabilityLikelihood.class),
                    AttributeRule.newStringRule(TreeTraitParserUtilities.TRAIT_NAME, true)
            };
        }

        @Override
        public String getParserDescription() {
            return "Gibbs sampler for the latent variables under the latent liability model when latent traits " +
                    "arise according to Gaussian extension of the process on the tree.";
        }

        @Override
        public Class getReturnType() {
            return ExtendedLatentLiabilityGibbsOperator.class;
        }

        @Override
        public String getParserName() {
            return EXTENDED_LATENT_GIBBS;
        }
    };

    @Override
    public String getReport() {
        CompoundParameter dataParameter = dataModel.getParameter();
        int nTaxa = dataParameter.getParameterCount();
        int nTraits = dataParameter.getParameter(0).getDimension();

        int reps = 20000;
        int m = dataParameter.getDimension();
        double[] mean = new double[m];
        double[] min = new double[m];
        double[] max = new double[m];

        Arrays.fill(min, Double.POSITIVE_INFINITY);
        Arrays.fill(max, Double.NEGATIVE_INFINITY);


        for (int i = 0; i < reps; i++) {
            doOperation();
            for (int j = 0; j < m; j++) {
                double value = dataParameter.getParameterValue(j);
                mean[j] += value;

                if (value < min[j]) {
                    min[j] = value;
                }
                if (value > max[j]) {
                    max[j] = value;
                }
            }
        }

        for (int i = 0; i < m; i++) {
            mean[i] /= reps;
        }

        StringBuilder sb = new StringBuilder();
        sb.append(EXTENDED_LATENT_GIBBS + " report:\n");

        for (int taxon = 0; taxon < nTaxa; taxon++) {
            Parameter taxonParameter = dataParameter.getParameter(taxon);
            String taxonID = taxonParameter.getId();
            int offset = taxon * nTraits;


            sb.append("\t" + taxonID + " mean:\t");
            for (int trait = 0; trait < nTraits; trait++) {
                sb.append(mean[offset + trait]);
                sb.append(" ");
            }
            sb.append("\n\t" + taxonID + " minimum:\t");
            for (int trait = 0; trait < nTraits; trait++) {
                sb.append(min[offset + trait]);
                sb.append(" ");
            }
            sb.append("\n\t" + taxonID + " maximum:\t");
            for (int trait = 0; trait < nTraits; trait++) {
                sb.append(max[offset + trait]);
                sb.append(" ");
            }
            sb.append("\n\n");


        }

        sb.append("\n\n");

        return sb.toString();
    }
}

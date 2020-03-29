package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.preorder.ContinuousExtensionDelegate;
import dr.inference.model.MatrixParameterInterface;
import dr.math.distributions.WishartSufficientStatistics;
import dr.math.interfaces.ConjugateWishartStatisticsProvider;
import dr.xml.*;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import static dr.evomodel.treedatalikelihood.preorder.AbstractRealizedContinuousTraitDelegate.REALIZED_TIP_TRAIT;

/**
 * @author Gabriel Hassler
 * @author Marc A. Suchard
 */
public class RepeatedMeasuresWishartStatistics implements ConjugateWishartStatisticsProvider {


    private final RepeatedMeasuresTraitDataModel traitModel;
    private final Tree tree;
    private final TreeTrait tipTrait;
    private final String traitName;
    private final ContinuousExtensionDelegate extensionDelegate;
    private final ContinuousDataLikelihoodDelegate likelihoodDelegate;
    private final double[] outerProduct;
    private final int dimTrait;
    private final int nTaxa;
    private final double[] buffer;
    private boolean forceResample;

    public RepeatedMeasuresWishartStatistics(RepeatedMeasuresTraitDataModel traitModel,
                                             TreeDataLikelihood treeLikelihood,
                                             boolean forceResample) {
        this.traitModel = traitModel;
        this.traitName = traitModel.getTraitName();
        this.tree = treeLikelihood.getTree();

        this.tipTrait = treeLikelihood.getTreeTrait(REALIZED_TIP_TRAIT + "." + traitName);
        this.likelihoodDelegate = (ContinuousDataLikelihoodDelegate) treeLikelihood.getDataLikelihoodDelegate();
        this.extensionDelegate = traitModel.getExtensionDelegate(likelihoodDelegate, tipTrait, tree);

        this.dimTrait = traitModel.getTraitDimension();
        this.nTaxa = tree.getExternalNodeCount();

        this.outerProduct = new double[dimTrait * dimTrait];
        this.buffer = new double[nTaxa * dimTrait];

        this.forceResample = forceResample;

    }


    @Override
    public MatrixParameterInterface getPrecisionParameter() {
        return traitModel.getExtensionPrecision();
    }

    @Override
    public WishartSufficientStatistics getWishartStatistics() {

        if (forceResample) {
            likelihoodDelegate.fireModelChanged();
        }
        double[] treeValues = (double[]) tipTrait.getTrait(tree, null);
        double[] dataValues = extensionDelegate.getExtendedValues(treeValues);

        DenseMatrix64F XminusY = DenseMatrix64F.wrap(nTaxa, dimTrait, buffer);
        DenseMatrix64F X = DenseMatrix64F.wrap(nTaxa, dimTrait, treeValues);
        DenseMatrix64F Y = DenseMatrix64F.wrap(nTaxa, dimTrait, dataValues);

        CommonOps.subtract(X, Y, XminusY);

        DenseMatrix64F outerProductMat = DenseMatrix64F.wrap(dimTrait, dimTrait, outerProduct);

        CommonOps.multTransA(XminusY, XminusY, outerProductMat);


        return new WishartSufficientStatistics(nTaxa, outerProduct);
    }

    public void setForceResample(Boolean b) {
        forceResample = b;
    }


    private static final boolean DEBUG = false;


    public static final String RM_WISHART_STATISTICS = "repeatedMeasuresWishartStatistics";
    private static final String FORCE_RESAMPLE = "forceResample";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            TreeDataLikelihood dataLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);

            RepeatedMeasuresTraitDataModel traitModel =
                    (RepeatedMeasuresTraitDataModel) xo.getChild(RepeatedMeasuresTraitDataModel.class);

            boolean forceResample = xo.getAttribute(FORCE_RESAMPLE, true);

            return new RepeatedMeasuresWishartStatistics(traitModel, dataLikelihood, forceResample);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(RepeatedMeasuresTraitDataModel.class),
                new ElementRule(TreeDataLikelihood.class),
                AttributeRule.newBooleanRule(FORCE_RESAMPLE, true)
        };

        @Override
        public String getParserDescription() {
            return null;
        }

        @Override
        public Class getReturnType() {
            return null;
        }

        @Override
        public String getParserName() {
            return RM_WISHART_STATISTICS;
        }
    };
}



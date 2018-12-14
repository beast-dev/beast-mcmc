package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.CompoundParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.math.distributions.WishartStatistics;
import dr.math.distributions.WishartSufficientStatistics;
import dr.math.interfaces.ConjugateWishartStatisticsProvider;
import dr.math.matrixAlgebra.WrappedVector;
import dr.math.matrixAlgebra.missingData.MissingOps;
import dr.xml.*;
import org.ejml.data.DenseMatrix64F;

import static dr.evomodel.treedatalikelihood.preorder.AbstractRealizedContinuousTraitDelegate.REALIZED_TIP_TRAIT;

/**
 * @author Gabriel Hassler
 * @author Marc A. Suchard
 */
public class RepeatedMeasuresWishartStatistics implements ConjugateWishartStatisticsProvider {


    private final RepeatedMeasuresTraitDataModel traitModel;
    private final WishartStatisticsWrapper treeWishartStatistics;
    private final TreeDataLikelihood treeLikelihood;
    private final TreeTrait tipTrait;
    private final String traitName;
    private final RepeatedMeasuresTraitSimulator traitSimulator;

    public RepeatedMeasuresWishartStatistics(RepeatedMeasuresTraitDataModel traitModel,
                                             TreeDataLikelihood treeLikelihood,
                                             WishartStatisticsWrapper treeWishartStatistics) {
        this.traitModel = traitModel;
        this.treeWishartStatistics = treeWishartStatistics;
        this.treeLikelihood = treeLikelihood;
        this.traitName = traitModel.getTraitName();
        this.traitSimulator = new RepeatedMeasuresTraitSimulator(traitModel, treeLikelihood);

        this.tipTrait = treeLikelihood.getTreeTrait(REALIZED_TIP_TRAIT + "." + traitName);
    }

    @Override
    public WishartSufficientStatistics getWishartStatistics() {

        treeWishartStatistics.simulateMissingTraits();
        return getRepeatedMeasuresStatistics();
    }

    @Override
    public MatrixParameterInterface getPrecisionParameter() {

        assert (traitModel.getSamplingPrecision() instanceof MatrixParameterInterface);

        return (MatrixParameterInterface) traitModel.getSamplingPrecision();
    }

    private WishartSufficientStatistics getRepeatedMeasuresStatistics() {


        double[] tipValues = (double[]) tipTrait.getTrait(treeLikelihood.getTree(), null);

        if (DEBUG) {
            System.err.println("tipValues: " + new WrappedVector.Raw(tipValues));
        }

        traitSimulator.simulateMissingData(tipValues);

        CompoundParameter traitParameter = traitModel.getParameter();


        int n = traitParameter.getParameterCount();
        int dim = traitModel.getTraitDimension();

        DenseMatrix64F XminusY = MissingOps.wrap(tipValues, 0, n, dim);


        for (int i = 0; i < n; i++) {

            for (int j = 0; j < dim; j++) {

                double traitValue = traitParameter.getParameterValue(j, i);
                double test_values = XminusY.get(i, j);

                XminusY.set(i, j, XminusY.get(i, j) - traitValue);
            }


        }

        DenseMatrix64F XminusYtXminusY = new DenseMatrix64F(dim, dim);

        org.ejml.ops.CommonOps.multTransA(XminusY, XminusY, XminusYtXminusY);

        double[] outerProduct = new double[dim * dim];

        MissingOps.unwrap(XminusYtXminusY, outerProduct, 0);

        return new WishartSufficientStatistics(n, outerProduct);
    }


    private static final boolean DEBUG = false;


    private static final String RM_WISHART_STATISTICS = "repeatedMeasuresWishartStatistics";

    public static AbstractXMLObjectParser PARSER = new AbstractXMLObjectParser() {
        @Override
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            TreeDataLikelihood dataLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);
            RepeatedMeasuresTraitDataModel traitModel =
                    (RepeatedMeasuresTraitDataModel) xo.getChild(RepeatedMeasuresTraitDataModel.class);
            WishartStatisticsWrapper treeWishartStatistics =
                    (WishartStatisticsWrapper) xo.getChild(WishartStatisticsWrapper.class);


            return new RepeatedMeasuresWishartStatistics(traitModel, dataLikelihood, treeWishartStatistics);
        }

        @Override
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new ElementRule(RepeatedMeasuresTraitDataModel.class),
                new ElementRule(TreeDataLikelihood.class),
                new ElementRule(WishartStatisticsWrapper.class)
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



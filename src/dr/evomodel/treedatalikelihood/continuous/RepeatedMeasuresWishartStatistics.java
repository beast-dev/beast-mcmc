package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.TreeTrait;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.model.MatrixParameterInterface;
import dr.math.distributions.WishartSufficientStatistics;
import dr.math.interfaces.ConjugateWishartStatisticsProvider;
import dr.math.matrixAlgebra.WrappedVector;

import static dr.evomodel.treedatalikelihood.preorder.AbstractRealizedContinuousTraitDelegate.REALIZED_TIP_TRAIT;

/**
 * @author Gabriel Hassler
 * @author Marc A. Suchard
 */
public class RepeatedMeasuresWishartStatistics implements ConjugateWishartStatisticsProvider {


    // TODO Make this class implement ConjugateWishartStatisticsProvider
    // TODO then can use PrecisionMatrixGibbsOperator and avoid code duplication

    private final RepeatedMeasuresTraitDataModel traitModel;
    private final WishartStatisticsWrapper treeWishartStatistics;
    private final TreeDataLikelihood treeLikelihood;
    private final TreeTrait tipTrait;
    private final String traitName;

    public RepeatedMeasuresWishartStatistics(RepeatedMeasuresTraitDataModel traitModel,
                                             TreeDataLikelihood treeLikelihood,
                                             WishartStatisticsWrapper treeWishartStatistics) {
        this.traitModel = traitModel;
        this.treeWishartStatistics = treeWishartStatistics;
        this.treeLikelihood = treeLikelihood;
        this.traitName = traitModel.getTraitName();

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
        // TODO -- data augmentation (?), then compute outer products

        double[] tipValues = (double[]) tipTrait.getTrait(treeLikelihood.getTree(), null);
        if (DEBUG) {
            System.err.println("tipValues: " + new WrappedVector.Raw(tipValues));
        } // TODO This "should" work, but probably doesn't...


        return null;
    }

    private static final boolean DEBUG = true;
}

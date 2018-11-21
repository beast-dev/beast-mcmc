package dr.evomodel.treedatalikelihood.continuous;

import dr.inference.model.MatrixParameterInterface;
import dr.math.distributions.WishartSufficientStatistics;
import dr.math.interfaces.ConjugateWishartStatisticsProvider;

/**
 * @author Gabriel Hassler
 * @author Marc A. Suchard
 */
public class RepeatedMeasuresWishartStatistics implements ConjugateWishartStatisticsProvider {


    // TODO Make this class implement ConjugateWishartStatisticsProvider
    // TODO then can use PrecisionMatrixGibbsOperator and avoid code duplication

    private final RepeatedMeasuresTraitDataModel traitModel;
    private final WishartStatisticsWrapper treeWishartStatistics;

    public RepeatedMeasuresWishartStatistics(RepeatedMeasuresTraitDataModel traitModel,
                                             WishartStatisticsWrapper treeWishartStatistics) {
        this.traitModel = traitModel;
        this.treeWishartStatistics = treeWishartStatistics;
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

    WishartSufficientStatistics getRepeatedMeasuresStatistics() {
        // TODO -- probably does not need to involve data augmentation, as I see in RepeatedMeasuresTraitSimulator
        
        return null;
    }
}

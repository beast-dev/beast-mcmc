package dr.evomodel.continuous.hmc;

import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.inference.model.Parameter;

/**
 * @author Marc A. Suchard
 */
public class CubicOrderTreePrecisionTraitProductProvider extends TreePrecisionTraitProductProvider {

    public CubicOrderTreePrecisionTraitProductProvider(TreeDataLikelihood treeDataLikelihood,
                                                       ContinuousDataLikelihoodDelegate likelihoodDelegate) {
        super(treeDataLikelihood, likelihoodDelegate);
    }

    @Override
    public double[] getProduct(Parameter vector) {
        return expensiveProduct(vector);
    }

    @Override
    public double[] getMassVector() {

        double[][] treeTraitVariance = likelihoodDelegate.getTreeTraitVariance();
        final int dim = treeTraitVariance.length;

        double[] mass = new double[dim];
        for (int i = 0; i < dim; ++i) {
            mass[i] = treeTraitVariance[i][i];
        }

        return mass;

    }

    @Override
    public double getTimeScale() {


        double[][] treeTraitVariance = likelihoodDelegate.getTreeTraitVariance();
        final int dim = treeTraitVariance.length;

        double max = Double.MIN_VALUE;
        for (int i = 0; i < dim; ++i) {
            max = Math.max(max, treeTraitVariance[i][i]);
        }

        return Math.sqrt(max);
    }
}

package dr.evomodel.continuous.hmc;

import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.inference.model.Model;
import dr.inference.model.Parameter;

/**
 * @author Marc A. Suchard
 */
public class CubicOrderTreePrecisionTraitProductProvider extends TreePrecisionTraitProductProvider {

    public CubicOrderTreePrecisionTraitProductProvider(TreeDataLikelihood treeDataLikelihood,
                                                       ContinuousDataLikelihoodDelegate likelihoodDelegate) {
        super(treeDataLikelihood, likelihoodDelegate);
        this.isPrecisionKnown = false;
    }

    @Override
    public double[] getProduct(Parameter vector) {
        return expensiveProduct(vector, getTreeTraitPrecision());
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

    private double[][] getTreeTraitPrecision() {

        if (!DO_CACHE) {
            treeTraitPrecision = likelihoodDelegate.getTreeTraitPrecision();
        } else if (!isPrecisionKnown) {
            treeTraitPrecision = likelihoodDelegate.getTreeTraitPrecision();
            isPrecisionKnown = true;
        }

        return treeTraitPrecision;
    }

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == tree) {
            isPrecisionKnown = false;
        }
    }

    @Override
    protected void storeState() {
        if (DO_CACHE) {
            final int dim = treeTraitPrecision.length;
            if (savedTreeTraitPrecision == null) {
                savedTreeTraitPrecision = new double[dim][dim];
            }

            for (int i = 0; i < dim; ++i) {
                System.arraycopy(treeTraitPrecision[i], 0, savedTreeTraitPrecision[i], 0, dim);
            }

            savedIsPrecisionKnown = isPrecisionKnown;
        }
    }

    @Override
    protected void restoreState() {
        if (DO_CACHE) {
            double[][] tmp = treeTraitPrecision;
            treeTraitPrecision = savedTreeTraitPrecision;
            savedTreeTraitPrecision = tmp;

            isPrecisionKnown = savedIsPrecisionKnown;
        }
    }

    private static final boolean DO_CACHE = true;

    private boolean isPrecisionKnown;
    private boolean savedIsPrecisionKnown;
    private double[][] treeTraitPrecision;
    private double[][] savedTreeTraitPrecision;
}

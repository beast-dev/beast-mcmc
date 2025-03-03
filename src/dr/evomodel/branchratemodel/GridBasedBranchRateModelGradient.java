package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evomodel.treedatalikelihood.discrete.DiscreteTraitBranchRateGradient;
import dr.inference.model.Parameter;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;

import java.util.Arrays;

public class GridBasedBranchRateModelGradient extends DiscreteTraitBranchRateGradient {

    private GridBasedBranchRateModel branchRateModel;

    public GridBasedBranchRateModelGradient(String traitName,
                                            TreeDataLikelihood treeDataLikelihood,
                                            BeagleDataLikelihoodDelegate likelihoodDelegate,
                                            Parameter rateParameter,
                                            GridBasedBranchRateModel gridBasedBranchRateModel,
                                            boolean useHessian) {
        super(traitName, treeDataLikelihood, likelihoodDelegate, rateParameter, useHessian);
        this.branchRateModel = gridBasedBranchRateModel;
    }

    @Override
    public double[] getGradientLogDensity() {

        long startTime;
        if (COUNT_TOTAL_OPERATIONS) {
            startTime = System.nanoTime();
        }

        double[] result = new double[branchRateModel.getGridPoints().getDimension() + 1];
        Arrays.fill(result, 0.0);

        double[] gradient = super.getGradientLogDensity();
//        int v = 0;
        for (int i = 0; i < tree.getNodeCount(); ++i) {
            final NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node)) {
                final int destinationIndex = getParameterIndexFromNode(node); //TODO  check this
                final double branchlengthDerivative = gradient[destinationIndex];
                for (int gridIndex = 0; gridIndex < branchRateModel.getGridPoints().getDimension() + 1; ++gridIndex) {
                    double suffStat = branchRateModel.getSufficientStatistics(gridIndex + (branchRateModel.getGridPoints().getDimension() + 1) * node.getNumber());
                    if (suffStat != 0) {
                        result[gridIndex] += branchlengthDerivative * suffStat;
                    }
                }
//                v++;
            }
        }
        if (COUNT_TOTAL_OPERATIONS) {
            ++getGradientLogDensityCount;
            long endTime = System.nanoTime();
            totalGradientTime += (endTime - startTime) / 1000000;
        }

        return result;
    }

    protected static final boolean COUNT_TOTAL_OPERATIONS = true;
    private long getGradientLogDensityCount = 0;
    private long totalGradientTime = 0;
}
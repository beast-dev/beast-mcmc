package dr.evomodel.branchratemodel;

import dr.evolution.tree.NodeRef;
import dr.evomodel.treedatalikelihood.discrete.DiscreteTraitBranchRateGradient;
import dr.inference.model.Parameter;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;

import java.util.Arrays;

public class GridBasedBranchRateModelGradient extends DiscreteTraitBranchRateGradient {

    private final GridBasedBranchRateModel branchRateModel;
    private final int nRates;

    public GridBasedBranchRateModelGradient(String traitName,
                                            TreeDataLikelihood treeDataLikelihood,
                                            BeagleDataLikelihoodDelegate likelihoodDelegate,
                                            Parameter rateParameter,
                                            GridBasedBranchRateModel gridBasedBranchRateModel,
                                            boolean useHessian) {
        super(traitName, treeDataLikelihood, likelihoodDelegate, rateParameter, useHessian);
        this.branchRateModel = gridBasedBranchRateModel;
        this.nRates = branchRateModel.getGridPoints().getDimension() + 1;
    }

    @Override
    public double[] getGradientLogDensity() {

        long startTime;
        if (COUNT_TOTAL_OPERATIONS) {
            startTime = System.nanoTime();
        }

        double[] result = new double[nRates];
        Arrays.fill(result, 0.0);

        double[] gradient = super.getGradientLogDensity();
//        int v = 0;
        int minGridIndex = 0;
        for (int i = 0; i < tree.getNodeCount(); ++i) {
            int nodeIndex = branchRateModel.getOrderedNodesIndexes(i);
            final NodeRef node = tree.getNode(i);
            if (!tree.isRoot(node)) {
                final int destinationIndex = getParameterIndexFromNode(node); // TODO  check this
                final double branchlengthDerivative = gradient[destinationIndex];
                double tChild = tree.getNodeHeight(tree.getNode(nodeIndex));
                double tParent = tree.getNodeHeight(tree.getParent(tree.getNode(nodeIndex)));

                while (minGridIndex < nRates - 1 && branchRateModel.getGridPoint(minGridIndex) < tChild) {
                    minGridIndex++;
                }
                int gridIndex = minGridIndex;

                if (gridIndex < nRates - 1 && branchRateModel.getGridPoint(gridIndex) < tParent) {
                    while (gridIndex < nRates - 1 && branchRateModel.getGridPoint(gridIndex) < tParent) {
                        double suffStat = branchRateModel.getSufficientStatistic(gridIndex + (nRates) * node.getNumber());
                        result[gridIndex] += branchlengthDerivative * suffStat;
                        gridIndex++;
                    }
                }
                double suffStat = branchRateModel.getSufficientStatistic(gridIndex + (nRates) * node.getNumber());
                result[gridIndex] += branchlengthDerivative * suffStat;
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
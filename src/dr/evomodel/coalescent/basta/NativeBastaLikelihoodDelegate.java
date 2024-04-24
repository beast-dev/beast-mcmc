package dr.evomodel.coalescent.basta;

import dr.evolution.alignment.PatternList;
import dr.evolution.tree.Tree;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.inference.model.Parameter;

import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class NativeBastaLikelihoodDelegate extends BastaLikelihoodDelegate.AbstractBastaLikelihoodDelegate {

    private final NativeBastaJniWrapper jni;

    public NativeBastaLikelihoodDelegate(String name,
                                         Tree tree,
                                         int stateCount,
                                         boolean transpose) {
        super(name, tree, stateCount, transpose);
        jni = NativeBastaJniWrapper.getBastaJniWrapper();
    }

    @Override
    String getStamp() { return "native"; }


    @Override
    protected void computeBranchIntervalOperations(List<Integer> intervalStarts,
                                                   List<BranchIntervalOperation> branchIntervalOperations) {
        if (branchIntervalOperations != null) {
            for (BranchIntervalOperation operation : branchIntervalOperations) {
                System.err.println(operation.toString());
            }
        }
    }

    @Override
    protected void computeTransitionProbabilityOperations(List<TransitionMatrixOperation> matrixOperations) {
        if (matrixOperations != null) {
            for (TransitionMatrixOperation operation : matrixOperations) {
                System.err.println(operation.toString());
            }
        }
    }

    @Override
    protected double computeCoalescentIntervalReduction(List<Integer> intervalStarts,
                                                        List<BranchIntervalOperation> branchIntervalOperations) {
        if (intervalStarts != null) {
            for (int start : intervalStarts) {
                System.err.println(start);
            }
        }

        return 0.0;
    }

    @Override
    protected void computeBranchIntervalOperationsGrad(List<Integer> intervalStarts, List<TransitionMatrixOperation> matrixOperations, List<BranchIntervalOperation> branchIntervalOperations) {

    }

    @Override
    protected void computeTransitionProbabilityOperationsGrad(List<TransitionMatrixOperation> matrixOperations) {

    }

    @Override
    protected double[][] computeCoalescentIntervalReductionGrad(List<Integer> intervalStarts, List<BranchIntervalOperation> branchIntervalOperations) {
        return new double[0][];
    }

    @Override
    protected double[] computeCoalescentIntervalReductionGradPopSize(List<Integer> intervalStarts, List<BranchIntervalOperation> branchIntervalOperations) {
        return new double[0];
    }
}

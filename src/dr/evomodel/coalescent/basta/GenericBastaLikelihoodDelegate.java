package dr.evomodel.coalescent.basta;

import dr.evolution.tree.Tree;
import dr.evomodel.substmodel.EigenDecomposition;
import dr.math.matrixAlgebra.WrappedVector;

import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class GenericBastaLikelihoodDelegate extends BastaLikelihoodDelegate.AbstractBastaLikelihoodDelegate {

    private static final boolean PRINT_COMMANDS = true;

    private final double[] partials;
    private final double[] matrices;
    private final double[] sizes;
    private final double[] coalescent;

    private final double[] temp;

    private final EigenDecomposition[] decompositions; // TODO flatten?

    public GenericBastaLikelihoodDelegate(String name,
                                          Tree tree,
                                          int stateCount) {
        super(name, tree, stateCount);

        this.partials = new double[maxNumCoalescentIntervals * tree.getNodeCount() * stateCount]; // TODO much too large
        this.matrices = new double[maxNumCoalescentIntervals * stateCount * stateCount]; // TODO much too small (except for strict-clock
        this.coalescent = new double[maxNumCoalescentIntervals];
        this.sizes = new double[2 * stateCount];
        this.decompositions = new EigenDecomposition[1];

        this.temp = new double[stateCount * stateCount];
    }

    @Override
    protected double computeBranchIntervalOperations(List<BranchIntervalOperation> branchIntervalOperations) {

        for (BranchIntervalOperation operation : branchIntervalOperations) { // TODO execute parallel by subIntervalNumber or executionOrder
            peelPartials(
                    partials, operation.outputBuffer,
                    operation.inputBuffer1, operation.inputBuffer2,
                    matrices,
                    operation.inputMatrix1, operation.inputMatrix2,
                    coalescent, operation.subIntervalNumber,
                    stateCount);

            if (PRINT_COMMANDS) {
                System.err.println(operation + " " +
                        new WrappedVector.Raw(partials, operation.outputBuffer * stateCount, stateCount));
            }
        }

        return 0.0;
    }

    @Override
    protected double computeTransitionProbabilityOperations(List<TransitionMatrixOperation> matrixOperations) {

        for (TransitionMatrixOperation operation : matrixOperations) { // TODO execute in parallel
            computeTransitionProbabilities(
                    operation.time,
                    matrices, operation.outputBuffer * stateCount * stateCount,
                    decompositions[operation.decompositionBuffer],
                    stateCount, temp);

            if (PRINT_COMMANDS) {
                System.err.println(operation + " " +
                        new WrappedVector.Raw(matrices,
                                operation.outputBuffer * stateCount * stateCount, stateCount * stateCount));
            }
        }

        return 0.0;
    }

    @Override
    protected double computeCoalescentIntervalReduction(List<Integer> intervalStarts,
                                                        List<BranchIntervalOperation> branchIntervalOperations) {
        for (int start : intervalStarts) { // TODO execute in parallel
            System.err.println(start);
        }

        return 0.0;
    }

    @Override
    public void setPartials(int index, double[] partials) {
        assert partials.length == stateCount;

        System.arraycopy(partials, 0, this.partials, index * stateCount, stateCount);
    }

    @Override
    public void setEigenDecomposition(int index, EigenDecomposition decomposition) {
        decompositions[index] = decomposition;
    }

    @Override
    public void setPopulationSizes(int index, double[] sizes) {
        assert sizes.length == stateCount;
        assert index == 0; // TODO generalize

        System.arraycopy(sizes, 0, this.sizes, index * stateCount, stateCount);
    }

    private static void peelPartials(double[] partials,
                                     int resultOffset,
                                     int leftPartialOffset, int rightPartialOffset,
                                     double[] matrices,
                                     int leftMatrixOffset, int rightMatrixOffset,
                                     double[] probability,
                                     int probabilityOffset,
                                     int stateCount) {

        resultOffset *= stateCount;

        // Handle left
        leftPartialOffset *= stateCount;
        leftMatrixOffset *= stateCount * stateCount;

        for (int i = 0; i < stateCount; ++i) {
            double sum = 0.0;
            for (int j = 0; j < stateCount; ++j) {
                sum += matrices[leftMatrixOffset + i * stateCount + j] *
                        partials[leftPartialOffset + j];
            }
            partials[resultOffset + i] = sum;
        }
        if (rightPartialOffset >= 0) {
            // Handle right
            rightPartialOffset *= stateCount;
            rightMatrixOffset *= stateCount * stateCount;

            for (int i = 0; i < stateCount; ++i) {
                double sum = 0.0;
                for (int j = 0; j < stateCount; ++j) {
                    sum += matrices[rightMatrixOffset + i * stateCount + j] *
                            partials[rightPartialOffset + j];
                }
                partials[resultOffset + i] *= sum;
            }

            probability[probabilityOffset] = 1.0; // TODO
        }
    }

    private static void computeTransitionProbabilities(double distance,
                                                       double[] matrix,
                                                       int matrixOffset,
                                                       EigenDecomposition eigen,
                                                       int stateCount,
                                                       double[] iexp) {

        assert matrix.length >= matrixOffset + stateCount * stateCount;
        assert iexp.length == stateCount * stateCount;

        boolean real = eigen.getEigenValues().length == stateCount;

        double[] Evec = eigen.getEigenVectors();
        double[] Eval = eigen.getEigenValues();
        double[] Ievc = eigen.getInverseEigenVectors();

        for (int i = 0; i < stateCount; i++) {

            if (real || Eval[stateCount + i] == 0) {
                // 1x1 block
                double temp = Math.exp(distance * Eval[i]);
                for (int j = 0; j < stateCount; j++) {
                    iexp[i * stateCount + j] = Ievc[i * stateCount + j] * temp;
                }
            } else {
                // 2x2 conjugate block
                // If A is 2x2 with complex conjugate pair eigenvalues a +/- bi, then
                // exp(At) = exp(at)*( cos(bt)I + \frac{sin(bt)}{b}(A - aI)).
                int i2 = i + 1;
                double b = Eval[stateCount + i];
                double expat = Math.exp(distance * Eval[i]);
                double expatcosbt = expat * Math.cos(distance * b);
                double expatsinbt = expat * Math.sin(distance * b);

                for (int j = 0; j < stateCount; j++) {
                    iexp[i * stateCount + j] = expatcosbt * Ievc[i * stateCount + j] +
                            expatsinbt * Ievc[i2 * stateCount + j];
                    iexp[i2 * stateCount + j] = expatcosbt * Ievc[i2 * stateCount + j] -
                            expatsinbt * Ievc[i * stateCount + j];
                }
                i++; // processed two conjugate rows
            }
        }

        int u = matrixOffset;
        for (int i = 0; i < stateCount; i++) {
            for (int j = 0; j < stateCount; j++) {
                double temp = 0.0;
                for (int k = 0; k < stateCount; k++) {
                    temp += Evec[i * stateCount + k] * iexp[k * stateCount + j];
                }
                matrix[u] = Math.abs(temp);
                u++;
            }
        }
    }
}

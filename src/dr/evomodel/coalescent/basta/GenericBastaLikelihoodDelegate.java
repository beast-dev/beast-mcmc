package dr.evomodel.coalescent.basta;

import dr.evolution.tree.Tree;
import dr.evomodel.substmodel.EigenDecomposition;
import dr.math.matrixAlgebra.WrappedVector;

import java.util.Arrays;
import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class GenericBastaLikelihoodDelegate extends BastaLikelihoodDelegate.AbstractBastaLikelihoodDelegate {

    private final double[] partials;
    private final double[] matrices;
    private final double[] sizes;
    private final double[] coalescent;

    private final double[] e;
    private final double[] f;
    private final double[] g;
    private final double[] h;

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

        this.e = new double[maxNumCoalescentIntervals * stateCount];
        this.f = new double[maxNumCoalescentIntervals * stateCount];
        this.g = new double[maxNumCoalescentIntervals * stateCount];
        this.h = new double[maxNumCoalescentIntervals * stateCount];

        this.temp = new double[stateCount * stateCount];
    }

    @Override
    protected void computeBranchIntervalOperations(List<Integer> intervalStarts,
                                                   List<BranchIntervalOperation> branchIntervalOperations) {

        Arrays.fill(coalescent, 0.0);

        for (int interval = 0; interval < intervalStarts.size() - 1; ++interval) { // execute in series by intervalNumber
            // TODO try grouping by executionOrder (unclear if more efficient, same total #)
            int start = intervalStarts.get(interval);
            int end = intervalStarts.get(interval + 1);

            computeInnerBranchIntervalOperations(branchIntervalOperations, start, end);
        }
    }

    protected void computeInnerBranchIntervalOperations(List<BranchIntervalOperation> branchIntervalOperations,
                                                      int start, int end) {

        for (int i = start; i < end; ++i) {
            BranchIntervalOperation operation = branchIntervalOperations.get(i);

            peelPartials(
                    partials, operation.outputBuffer,
                    operation.inputBuffer1, operation.inputBuffer2,
                    matrices,
                    operation.inputMatrix1, operation.inputMatrix2,
                    operation.accBuffer1, operation.accBuffer2,
                    coalescent, operation.intervalNumber,
                    sizes, 0,
                    stateCount);

            if (PRINT_COMMANDS) {
                System.err.println(operation + " " +
                        new WrappedVector.Raw(partials, operation.outputBuffer * stateCount, stateCount));
            }
        }
    }

    @Override
    protected void computeTransitionProbabilityOperations(List<TransitionMatrixOperation> matrixOperations) {
        computeInnerTransitionProbabilityOperations(matrixOperations, 0, matrixOperations.size(), temp);
    }

    protected void computeInnerTransitionProbabilityOperations(List<TransitionMatrixOperation> matrixOperations,
                                                               int start, int end, double[] temp) {
        for (int i = start; i < end; ++i) {
            TransitionMatrixOperation operation = matrixOperations.get(i);
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
    }

    @Override
    protected double computeCoalescentIntervalReduction(List<Integer> intervalStarts,
                                                        List<BranchIntervalOperation> branchIntervalOperations) {

        Arrays.fill(e, 0.0);
        Arrays.fill(f, 0.0);
        Arrays.fill(g, 0.0);
        Arrays.fill(h, 0.0);

        for (int interval = 0; interval < intervalStarts.size() - 1; ++interval) { // TODO execute in parallel (no race conditions)
            int start = intervalStarts.get(interval);
            int end = intervalStarts.get(interval + 1);

            for (int i = start; i < end; ++i) { // TODO execute in parallel (has race conditions)
                BranchIntervalOperation operation = branchIntervalOperations.get(i);
                reduceWithinInterval(e, f, g, h, partials,
                        operation.inputBuffer1, operation.inputBuffer2,
                        operation.accBuffer1, operation.accBuffer2,
                        operation.intervalNumber,
                        stateCount);

            }
        }

        double logL = 0.0;
        for (int i = 0; i < intervalStarts.size() - 1; ++i) { // TODO execute in parallel
            BranchIntervalOperation operation = branchIntervalOperations.get(intervalStarts.get(i));

            logL += reduceAcrossIntervals(e, f, g, h,
                    operation.intervalNumber, operation.intervalLength,
                    sizes, coalescent, stateCount);
        }

        if (PRINT_COMMANDS) {
            System.err.println("logL = " + logL + "\n");
        }

        return logL;
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

        System.arraycopy(sizes, 0, this.sizes, index * stateCount, stateCount);
    }

    private static void peelPartials(double[] partials,
                                     int resultOffset,
                                     int leftPartialOffset, int rightPartialOffset,
                                     double[] matrices,
                                     int leftMatrixOffset, int rightMatrixOffset,
                                     int leftAccOffset, int rightAccOffset,
                                     double[] probability, int probabilityOffset,
                                     double[] sizes, int sizesOffset,
                                     int stateCount) {

        resultOffset *= stateCount;

        // Handle left
        leftPartialOffset *= stateCount;
        leftMatrixOffset *= stateCount * stateCount;

        for (int i = 0; i < stateCount; ++i) {
            double sum = 0.0;
            for (int j = 0; j < stateCount; ++j) {
                sum += matrices[leftMatrixOffset + i * stateCount + j] * partials[leftPartialOffset + j];
            }
            partials[resultOffset + i] = sum;
        }

        if (rightPartialOffset >= 0) {
            // Handle right
            rightPartialOffset *= stateCount;
            rightMatrixOffset *= stateCount * stateCount;

            leftAccOffset *= stateCount;
            rightAccOffset *= stateCount;

            sizesOffset *= sizesOffset * stateCount;

            double prob = 0.0;
            for (int i = 0; i < stateCount; ++i) {
                double right = 0.0;
                for (int j = 0; j < stateCount; ++j) {
                    right += matrices[rightMatrixOffset + i * stateCount + j] * partials[rightPartialOffset + j];
                }
                // entry = left * right * size
                double left = partials[resultOffset + i];
                double entry = left * right / sizes[sizesOffset + i];

                partials[resultOffset + i] = entry;
                partials[leftAccOffset + i] = left;
                partials[rightAccOffset + i] = right;

                prob += entry;
            }

            for (int i = 0; i < stateCount; ++i) {
                partials[resultOffset + i] /= prob;
            }

            probability[probabilityOffset] = prob;
        }

        // TODO rescale?
    }

    private static void computeTransitionProbabilities(double distance,
                                                       double[] matrix,
                                                       int matrixOffset,
                                                       EigenDecomposition eigen,
                                                       int stateCount,
                                                       double[] iexp) {

        assert matrix.length >= matrixOffset * stateCount * stateCount;
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

    private static double reduceAcrossIntervals(double[] e, double[] f, double[] g, double[] h,
                                                int interval, double length,
                                                double[] sizes, double[] coalescent,
                                                int stateCount) {

        int offset = interval * stateCount;

        double sum = 0.0;
        for (int k = 0; k < stateCount; ++k) {
            sum += (e[offset + k] * e[offset + k] - f[offset + k] +
                    g[offset + k] * g[offset + k] - h[offset + k]) / sizes[k];
        }

        double logL = -length * sum / 4;

        double prob = coalescent[interval];
        if (prob != 0.0) {
            logL += Math.log(prob);
        }

        return logL;
    }

    private static void reduceWithinInterval(double[] e, double[] f, double[] g, double[] h,
                                             double[] partials,
                                             int startBuffer1, int startBuffer2,
                                             int endBuffer1, int endBuffer2,
                                             int interval, int stateCount) {

        interval *= stateCount;

        startBuffer1 *= stateCount;
        endBuffer1 *= stateCount;

        for (int i = 0; i < stateCount; ++i) {
            double startP = partials[startBuffer1 + i];
            e[interval + i] += startP;
            f[interval + i] += startP * startP;

            double endP = partials[endBuffer1 + i];
            g[interval + i] += endP;
            h[interval + i] += endP * endP;
        }

        if (startBuffer2 >= 0) {
            startBuffer2 *= stateCount;
            endBuffer2 *= stateCount;

            for (int i = 0; i < stateCount; ++i) {
                double startP = partials[startBuffer2 + i];
                e[interval + i] += startP;
                f[interval + i] += startP * startP;

                double endP = partials[endBuffer2 + i];
                g[interval + i] += endP;
                h[interval + i] += endP * endP;
            }
        }
    }
}

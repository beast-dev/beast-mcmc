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

    final static class InternalStorage {

        final double[] partials;
        final double[] matrices;
        final double[] sizes;
        final double[] coalescent;

        final double[] e;
        final double[] f;
        final double[] g;
        final double[] h;

        final EigenDecomposition[] decompositions; // TODO flatten?

        public InternalStorage(int maxNumCoalescentIntervals, int treeNodeCount, int stateCount) {

            this.partials = new double[maxNumCoalescentIntervals * (treeNodeCount + 1) * stateCount]; // TODO much too large
            this.matrices = new double[maxNumCoalescentIntervals * stateCount * stateCount]; // TODO much too small (except for strict-clock)
            this.coalescent = new double[maxNumCoalescentIntervals];
            this.sizes = new double[2 * stateCount];

            this.e = new double[maxNumCoalescentIntervals * stateCount];
            this.f = new double[maxNumCoalescentIntervals * stateCount];
            this.g = new double[maxNumCoalescentIntervals * stateCount];
            this.h = new double[maxNumCoalescentIntervals * stateCount];

            this.decompositions = new EigenDecomposition[1];
        }
    }

    static class GradientInternalStorage {

        // TODO greatly simplify
        private final double[][][] partials;
        private final double[][][] matrices;
        private final double[][][] coalescent;
        private final double[][][] e;
        private final double[][][] f;
        private final double[][][] g;
        private final double[][][] h;
        private final double[][] partialsGradPopSize;
        private final double[][] coalescentGradPopSize;
        private final double[][] eGradPopSize;
        private final double[][] fGradPopSize;
        private final double[][] gGradPopSize;
        private final double[][] hGradPopSize;

        GradientInternalStorage(int maxNumCoalescentIntervals, int treeNodeCount, int stateCount, boolean transpose) {

            // TODO note this should be split into separate substitutionModel / popSizes contains

            this.partials = new double[stateCount][stateCount][maxNumCoalescentIntervals * (treeNodeCount + 1) * stateCount];
            if (!transpose) {
                this.matrices = new double[stateCount][stateCount][maxNumCoalescentIntervals * stateCount * stateCount];
            } else {
                this.matrices = null;
            }
            this.coalescent = new double[stateCount][stateCount][maxNumCoalescentIntervals];

            this.e = new double[stateCount][stateCount][maxNumCoalescentIntervals * stateCount];
            this.f = new double[stateCount][stateCount][maxNumCoalescentIntervals * stateCount];
            this.g = new double[stateCount][stateCount][maxNumCoalescentIntervals * stateCount];
            this.h = new double[stateCount][stateCount][maxNumCoalescentIntervals * stateCount];

            this.partialsGradPopSize = new double[stateCount][maxNumCoalescentIntervals * (treeNodeCount + 1) * stateCount];
            this.coalescentGradPopSize = new double[stateCount][maxNumCoalescentIntervals];

            this.eGradPopSize = new double[stateCount][maxNumCoalescentIntervals * stateCount];
            this.fGradPopSize = new double[stateCount][maxNumCoalescentIntervals * stateCount];
            this.gGradPopSize = new double[stateCount][maxNumCoalescentIntervals * stateCount];
            this.hGradPopSize = new double[stateCount][maxNumCoalescentIntervals * stateCount];
        }
    }

    private final InternalStorage storage;
    private GradientInternalStorage gradientStorage;
    private final double[] temp;

    public GenericBastaLikelihoodDelegate(String name,
                                          Tree tree,
                                          int stateCount,
                                          boolean transpose) {
        super(name, tree, stateCount, transpose);

        this.storage = new InternalStorage(maxNumCoalescentIntervals, tree.getNodeCount(), stateCount);
        this.gradientStorage = null; //new GradientInternalStorage(maxNumCoalescentIntervals, tree.getNodeCount(), stateCount);

        this.temp = new double[stateCount * stateCount];
    }

    public InternalStorage getInternalStorage() {
        return storage;
    }

    @Override
    protected void allocateGradientMemory() {
        if (gradientStorage == null) {
            gradientStorage = new GradientInternalStorage(maxNumCoalescentIntervals, tree.getNodeCount(),
                    stateCount, transpose);
        }
    }

    @Override
    protected void computeBranchIntervalOperations(List<Integer> intervalStarts,
                                                   List<BranchIntervalOperation> branchIntervalOperations,
                                                   List<TransitionMatrixOperation> matrixOperations,
                                                   Mode mode) {

        Arrays.fill(storage.coalescent, 0.0);

        for (int interval = 0; interval < intervalStarts.size() - 1; ++interval) { // execute in series by intervalNumber
            // TODO try grouping by executionOrder (unclear if more efficient, same total #)
            int start = intervalStarts.get(interval);
            int end = intervalStarts.get(interval + 1);

            computeInnerBranchIntervalOperations(branchIntervalOperations, matrixOperations, start, end, mode);
        }

        if (PRINT_COMMANDS) {
            System.err.println("coalescent: " + new WrappedVector.Raw(storage.coalescent));
        }
    }

    protected void computeInnerBranchIntervalOperations(List<BranchIntervalOperation> branchIntervalOperations,
                                                        List<TransitionMatrixOperation> matrixOperations,
                                                        int start, int end,
                                                        Mode mode) {

        for (int i = start; i < end; ++i) {
            BranchIntervalOperation operation = branchIntervalOperations.get(i);

            if (mode == Mode.LIKELIHOOD) {
                peelPartials(
                        storage.partials, operation.outputBuffer,
                        operation.inputBuffer1, operation.inputBuffer2,
                        storage.matrices,
                        operation.inputMatrix1, operation.inputMatrix2,
                        operation.accBuffer1, operation.accBuffer2,
                        storage.coalescent, operation.intervalNumber,
                        storage.sizes, 0,
                        stateCount);

            } else if (mode == Mode.GRADIENT) {

                TransitionMatrixOperation matrixOperation = matrixOperations.get(operation.intervalNumber);

                peelPartials(
                        storage.partials, operation.outputBuffer,
                        operation.inputBuffer1, operation.inputBuffer2,
                        storage.matrices,
                        operation.inputMatrix1, operation.inputMatrix2,
                        operation.accBuffer1, operation.accBuffer2,
                        storage.coalescent, operation.intervalNumber,
                        storage.sizes, 0,
                        stateCount);

                peelPartialsGrad(
                        storage.partials, matrixOperation.time, operation.outputBuffer,
                        operation.inputBuffer1, operation.inputBuffer2,
                        storage.matrices,
                        operation.inputMatrix1, operation.inputMatrix2,
                        operation.accBuffer1, operation.accBuffer2,
                        storage.coalescent, operation.intervalNumber,
                        storage.sizes, 0,
                        stateCount);
            }

            if (PRINT_COMMANDS) {
                System.err.println(operation + " " +
                        new WrappedVector.Raw(storage.partials, operation.outputBuffer * stateCount, stateCount));
            }
        }
    }

    @Override
    protected void computeTransitionProbabilityOperations(List<TransitionMatrixOperation> matrixOperations,
                                                          Mode mode) {
        computeInnerTransitionProbabilityOperations(matrixOperations, 0, matrixOperations.size(), mode, temp);
    }

    protected void computeInnerTransitionProbabilityOperations(List<TransitionMatrixOperation> matrixOperations,
                                                               int start, int end, Mode mode, double[] temp) {
        for (int i = start; i < end; ++i) {
            TransitionMatrixOperation operation = matrixOperations.get(i);

            if (mode == Mode.LIKELIHOOD) {
                computeTransitionProbabilities(
                        operation.time,
                        storage.matrices, operation.outputBuffer * stateCount * stateCount,
                        storage.decompositions[operation.decompositionBuffer],
                        stateCount, temp);
            } else if (mode == Mode.GRADIENT) {
                computeTransitionProbabilitiesGrad(operation.time,
                        operation.outputBuffer * stateCount * stateCount);
            }

            if (PRINT_COMMANDS) {
                System.err.println(operation + " " +
                        new WrappedVector.Raw(storage.matrices,
                                operation.outputBuffer * stateCount * stateCount, stateCount * stateCount));
            }
        }
    }

    @Override
    String getStamp() { return "generic"; }

    interface Dispatch {

        void clearStorage();

        void reduceWithinInterval(BranchIntervalOperation operation);

        void reduceAcrossIntervals(BranchIntervalOperation operation, double[] out);

    }

    class LikelihoodDispatch implements Dispatch {

        @Override
        public void clearStorage() {
            Arrays.fill(storage.e, 0.0);
            Arrays.fill(storage.f, 0.0);
            Arrays.fill(storage.g, 0.0);
            Arrays.fill(storage.h, 0.0);
        }

        @Override
        public void reduceWithinInterval(BranchIntervalOperation operation) {
            GenericBastaLikelihoodDelegate.reduceWithinInterval(storage.e, storage.f, storage.g, storage.h,
                    storage.partials,
                    operation.inputBuffer1, operation.inputBuffer2,
                    operation.accBuffer1, operation.accBuffer2,
                    operation.intervalNumber,
                    stateCount);
        }

        @Override
        public void reduceAcrossIntervals(BranchIntervalOperation operation, double[] out) {
            GenericBastaLikelihoodDelegate.reduceAcrossIntervals(storage.e, storage.f, storage.g, storage.h,
                    operation.intervalNumber, operation.intervalLength,
                    storage.sizes, storage.coalescent, out, stateCount);
        }
    }

    class MigrationRateGradientDispatch implements Dispatch {

        @Override
        public void clearStorage() { }

        @Override
        public void reduceWithinInterval(BranchIntervalOperation operation) {
            reduceWithinIntervalGrad(storage.partials, gradientStorage.partials,
                    operation.inputBuffer1, operation.inputBuffer2,
                    operation.accBuffer1, operation.accBuffer2,
                    operation.intervalNumber,
                    stateCount);
        }

        @Override
        public void reduceAcrossIntervals(BranchIntervalOperation operation, double[] out) {
            reduceAcrossIntervalsGrad(storage.e, storage.f, storage.g, storage.h,
                    operation.intervalNumber, operation.intervalLength,
                    storage.sizes, storage.coalescent, stateCount, out);
        }
    }

    class PopulationSizesGradientDispatch implements Dispatch {

        @Override
        public void clearStorage() {
            Arrays.stream(gradientStorage.eGradPopSize).forEach(a -> Arrays.fill(a, 0));
            Arrays.stream(gradientStorage.fGradPopSize).forEach(a -> Arrays.fill(a, 0));
            Arrays.stream(gradientStorage.gGradPopSize).forEach(a -> Arrays.fill(a, 0));
            Arrays.stream(gradientStorage.hGradPopSize).forEach(a -> Arrays.fill(a, 0));
        }

        @Override
        public void reduceWithinInterval(BranchIntervalOperation operation) {
            reduceWithinIntervalGrad(storage.partials, gradientStorage.partials,
                    operation.inputBuffer1, operation.inputBuffer2,
                    operation.accBuffer1, operation.accBuffer2,
                    operation.intervalNumber,
                    stateCount);
        }

        @Override
        public void reduceAcrossIntervals(BranchIntervalOperation operation, double[] out) {
            reduceAcrossIntervalsGradPopSize(storage.e, storage.f, storage.g, storage.h,
                    operation.intervalNumber, operation.intervalLength,
                    storage.sizes, storage.coalescent, out, stateCount);
        }
    }


    @Override
    protected void computeCoalescentIntervalReduction(List<Integer> intervalStarts,
                                                        List<BranchIntervalOperation> branchIntervalOperations,
                                                        double[] out,
                                                        Mode mode,
                                                        StructuredCoalescentLikelihoodGradient.WrtParameter wrt) {

        final Dispatch dispatch;

        if (mode == Mode.LIKELIHOOD) {
            dispatch = new LikelihoodDispatch();
        } else if (mode == Mode.GRADIENT) {
            if (wrt == StructuredCoalescentLikelihoodGradient.WrtParameter.MIGRATION_RATE) {
                dispatch = new MigrationRateGradientDispatch();
            } else if (wrt == StructuredCoalescentLikelihoodGradient.WrtParameter.POPULATION_SIZE) {
                dispatch = new PopulationSizesGradientDispatch();
            } else {
                dispatch = null;
            }
        } else {
            throw new RuntimeException("Not yet implemented");
        }

        dispatchComputeCoalescentIntervalReduction(intervalStarts, branchIntervalOperations, out, dispatch);
    }

    protected void dispatchComputeCoalescentIntervalReduction(List<Integer> intervalStarts,
                                                        List<BranchIntervalOperation> branchIntervalOperations,
                                                        double[] out,
                                                        Dispatch dispatch) {

        dispatch.clearStorage();

        for (int interval = 0; interval < intervalStarts.size() - 1; ++interval) { // TODO execute in parallel (no race conditions)
            int start = intervalStarts.get(interval);
            int end = intervalStarts.get(interval + 1);

            for (int i = start; i < end; ++i) { // TODO execute in parallel (has race conditions)
                BranchIntervalOperation operation = branchIntervalOperations.get(i);
                dispatch.reduceWithinInterval(operation);
            }
        }

        for (int i = 0; i < intervalStarts.size() - 1; ++i) { // TODO execute in parallel
            BranchIntervalOperation operation = branchIntervalOperations.get(intervalStarts.get(i));
            dispatch.reduceAcrossIntervals(operation, out);
        }
    }

    private void peelPartialsGrad(double[] partials,
                                  double distance, int resultOffset,
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
        leftAccOffset *= stateCount;

        for (int a = 0; a < stateCount; ++a) {
            for (int b = 0; b < stateCount; ++b) {
                for (int i = 0; i < stateCount; ++i) {
                    double sum = 0.0;
                    if (transpose && i == b) {
                        sum += partials[leftAccOffset + a] * distance;
                    }
                    if (!transpose) {
                        for (int j = 0; j < stateCount; ++j) {
                            sum += gradientStorage.matrices[a][b][leftMatrixOffset + i * stateCount + j] * partials[leftPartialOffset + j];
                        }
                    }
                    for (int j = 0; j < stateCount; ++j) {
                        sum += matrices[leftMatrixOffset + i * stateCount + j] * gradientStorage.partials[a][b][leftPartialOffset + j];
                    }
                    gradientStorage.partials[a][b][resultOffset + i] = sum;

                    //throw new RuntimeException("Function should not depend on `transpose`");
                }
            }
        }

        for (int a = 0; a < stateCount; ++a) {
            for (int i = 0; i < stateCount; ++i) {
                double sum = 0.0;
                for (int j = 0; j < stateCount; ++j) {
                    sum += matrices[leftMatrixOffset + i * stateCount + j] * gradientStorage.partialsGradPopSize[a][leftPartialOffset + j];
                }
                gradientStorage.partialsGradPopSize[a][resultOffset + i] = sum;
            }
        }

        if (rightPartialOffset >= 0) {
            // Handle right
            rightPartialOffset *= stateCount;
            rightMatrixOffset *= stateCount * stateCount;

            rightAccOffset *= stateCount;
            // TODO: check bug?
            sizesOffset *= sizesOffset * stateCount;

            for (int a = 0; a < stateCount; ++a) {
                for (int b = 0; b < stateCount; ++b) {
                    double J = probability[probabilityOffset];

                    // first half
                    double partial_J_ab = 0.0;
                    for (int i = 0; i < stateCount; ++i) {
                        double rightGrad = 0.0;
                        if (transpose && i == b) {
                            rightGrad += partials[rightAccOffset + a] * distance;
                        }

                        if (!transpose) {
                            for (int j = 0; j < stateCount; ++j) {
                                rightGrad += gradientStorage.matrices[a][b][rightMatrixOffset + i * stateCount + j] * partials[rightPartialOffset + j];
                            }
                        }
                        for (int j = 0; j < stateCount; ++j) {
                            rightGrad += matrices[rightMatrixOffset + i * stateCount + j] * gradientStorage.partials[a][b][rightPartialOffset + j];
                        }
                        double leftGrad = gradientStorage.partials[a][b][resultOffset + i];
                        double left = partials[leftAccOffset + i];
                        double right = partials[rightAccOffset + i];

                        double entry = (leftGrad * right + rightGrad * left) / sizes[sizesOffset + i];
                        partial_J_ab += entry;

                        gradientStorage.partials[a][b][resultOffset + i] = entry / J;
                        gradientStorage.partials[a][b][leftAccOffset + i] = leftGrad;
                        gradientStorage.partials[a][b][rightAccOffset + i] = rightGrad;

                        //throw new RuntimeException("Function should not depend on `transpose`");
                    }
                    // second half
                    for (int i = 0; i < stateCount; ++i) {
                        double entry = partials[resultOffset + i];
                        gradientStorage.partials[a][b][resultOffset + i] -= partial_J_ab * entry / J;
                    }
                    gradientStorage.coalescent[a][b][probabilityOffset] = partial_J_ab;
                }
            }

            for (int a = 0; a < stateCount; ++a) {
                double J = probability[probabilityOffset];
                // first half
                double partial_J_ab_PopSize = 0.0;
                for (int i = 0; i < stateCount; ++i) {
                    double rightGradPopSize = 0.0;
                    for (int j = 0; j < stateCount; ++j) {
                        rightGradPopSize += matrices[rightMatrixOffset + i * stateCount + j] *  gradientStorage.partialsGradPopSize[a][rightPartialOffset + j];
                    }
                    double leftGradPopSize = gradientStorage.partialsGradPopSize[a][resultOffset + i];
                    double left = partials[leftAccOffset + i];
                    double right = partials[rightAccOffset + i];

                    double entry = (leftGradPopSize * right + rightGradPopSize * left) / sizes[sizesOffset + i];
                    if (i == a){
                        entry += left * right;
                    }
                    partial_J_ab_PopSize += entry;

                    gradientStorage.partialsGradPopSize[a][resultOffset + i] = entry / J;
                    gradientStorage.partialsGradPopSize[a][leftAccOffset + i] = leftGradPopSize;
                    gradientStorage.partialsGradPopSize[a][rightAccOffset + i] = rightGradPopSize;
                }
                // second half
                for (int i = 0; i < stateCount; ++i) {
                    double entry = partials[resultOffset + i];
                    gradientStorage.partialsGradPopSize[a][resultOffset + i] -= partial_J_ab_PopSize * entry / J;
                }
                gradientStorage.coalescentGradPopSize[a][probabilityOffset] = partial_J_ab_PopSize;
            }
        }
    }

    private void computeTransitionProbabilitiesGrad(double distance, int matrixOffset) {
        for (int a = 0; a < stateCount; a++) {
            for (int b = 0; b < stateCount; b++) {
                for (int c = 0; c < stateCount; c++) {
                    for (int d = 0; d < stateCount; d++) { // TODO MAS: last loop unnecessary (also S^4 storage is unnecessary)
                        if (d == b) {
                            // TODO MAS: should these be cached at all? why not generate on the fly (t * matrices[])
                            gradientStorage.matrices[a][b][matrixOffset + c*stateCount + b] =  distance * storage.matrices[matrixOffset + c*stateCount + a];
                        } else {
                            gradientStorage.matrices[a][b][matrixOffset + c*stateCount + d] = 0; // TODO MAS: avoid caching (many) zeros
                        }
                    }
                }
            }
        }
    }

    private void reduceAcrossIntervalsGrad(double[] e, double[] f, double[] g, double[] h,
                                           int interval, double length,
                                           double[] sizes, double[] coalescent,
                                           int stateCount, double[] out) {

        int offset = interval * stateCount;

        for (int a = 0; a < stateCount; a++) {
            for (int b = 0; b < stateCount; b++) {
                double sum = 0.0;
                for (int k = 0; k < stateCount; ++k) {
                    sum += (2 * e[offset + k] * gradientStorage.e[a][b][offset + k] - gradientStorage.f[a][b][offset + k] +
                            2 * g[offset + k] * gradientStorage.g[a][b][offset + k] - gradientStorage.h[a][b][offset + k]) / sizes[k];
                }

                double element = -length * sum / 4;

                double J = coalescent[interval];
                if (J != 0.0) {
                    element += gradientStorage.coalescent[a][b][interval] / J;
                }

                out[a * stateCount + b] += element;
            }
        }
    }

    private void reduceAcrossIntervalsGradPopSize(double[] e, double[] f, double[] g, double[] h,
                                                 int interval, double length,
                                                 double[] sizes, double[] coalescent,
                                                 double[] out,
                                                 int stateCount) {

        int offset = interval * stateCount;
//        double[] grad = new double[stateCount];

        for (int a = 0; a < stateCount; a++) {
                double sum = 0.0;
                for (int k = 0; k < stateCount; ++k) {
                    sum += (2 * e[offset + k] * gradientStorage.eGradPopSize[a][offset + k] - gradientStorage.fGradPopSize[a][offset + k] +
                            2 * g[offset + k] * gradientStorage.gGradPopSize[a][offset + k] - gradientStorage.hGradPopSize[a][offset + k]) / sizes[k];
                    if (k == a) {
                        sum += (e[offset + k] * e[offset + k]) - f[offset + k] +
                                (g[offset + k] * g[offset + k]) - h[offset + k];
                    }
                }

                double element = -length * sum / 4;

                double J = coalescent[interval];
                if (J != 0.0) {
                    element += gradientStorage.coalescentGradPopSize[a][interval] / J;
                }

                out[a] += element;
            }
    }

    private void reduceWithinIntervalGrad(double[] partials, double[][][] partialsGrad,
                                          int startBuffer1, int startBuffer2,
                                          int endBuffer1, int endBuffer2,
                                          int interval, int stateCount)  {
        interval *= stateCount;

        startBuffer1 *= stateCount;
        endBuffer1 *= stateCount;

        for (int a = 0; a < stateCount; ++a) {
            for (int b = 0; b < stateCount; ++b) {
                for (int i = 0; i < stateCount; ++i) {
                    double startPGrad = partialsGrad[a][b][startBuffer1 + i];
                    double startP = partials[startBuffer1 + i];
                    gradientStorage.e[a][b][interval + i] += startPGrad;
                    gradientStorage.f[a][b][interval + i] += 2 * startP * startPGrad;

                    double endPGrad = partialsGrad[a][b][endBuffer1 + i];
                    double endP = partials[endBuffer1 + i];
                    gradientStorage.g[a][b][interval + i] += endPGrad;
                    gradientStorage.h[a][b][interval + i] += 2 * endP * endPGrad;
                }
            }
        }

        for (int a = 0; a < stateCount; ++a) {
                for (int i = 0; i < stateCount; ++i) {
                    double startPGradPopSize = gradientStorage.partialsGradPopSize[a][startBuffer1 + i];
                    double startP = partials[startBuffer1 + i];
                    gradientStorage.eGradPopSize[a][interval + i] += startPGradPopSize;
                    gradientStorage.fGradPopSize[a][interval + i] += 2 * startP * startPGradPopSize;

                    double endPGradPopSize = gradientStorage.partialsGradPopSize[a][endBuffer1 + i];
                    double endP = partials[endBuffer1 + i];
                    gradientStorage.gGradPopSize[a][interval + i] += endPGradPopSize;
                    gradientStorage.hGradPopSize[a][interval + i] += 2 * endP * endPGradPopSize;
                }
            }

        if (startBuffer2 >= 0) {
            startBuffer2 *= stateCount;
            endBuffer2 *= stateCount;
            for (int a = 0; a < stateCount; ++a) {
                for (int b = 0; b < stateCount; ++b) {
                    for (int i = 0; i < stateCount; ++i) {
                        double startPGrad = partialsGrad[a][b][startBuffer2 + i];
                        double startP = partials[startBuffer2 + i];
                        gradientStorage.e[a][b][interval + i] += startPGrad;
                        gradientStorage.f[a][b][interval + i] += 2 * startP * startPGrad;

                        double endPGrad = partialsGrad[a][b][endBuffer2 + i];
                        double endP = partials[endBuffer2 + i];
                        gradientStorage.g[a][b][interval + i] += endPGrad;
                        gradientStorage.h[a][b][interval + i] += 2 * endP * endPGrad;
                    }
                }
            }

            for (int a = 0; a < stateCount; ++a) {
                    for (int i = 0; i < stateCount; ++i) {
                        double startPGradPopSize = gradientStorage.partialsGradPopSize[a][startBuffer2 + i];
                        double startP = partials[startBuffer2 + i];
                        gradientStorage.eGradPopSize[a][interval + i] += startPGradPopSize;
                        gradientStorage.fGradPopSize[a][interval + i] += 2 * startP * startPGradPopSize;

                        double endPGradPopSize = gradientStorage.partialsGradPopSize[a][endBuffer2 + i];
                        double endP = partials[endBuffer2 + i];
                        gradientStorage.gGradPopSize[a][interval + i] += endPGradPopSize;
                        gradientStorage.hGradPopSize[a][interval + i] += 2 * endP * endPGradPopSize;
                    }

            }
        }
    }

    @Override
    public void setPartials(int index, double[] partials) {
        assert partials.length == stateCount;

        System.arraycopy(partials, 0, storage.partials, index * stateCount, stateCount);
    }

    @Override
    public void updateEigenDecomposition(int index, EigenDecomposition decomposition, boolean flip) {

        if (transpose) {
            decomposition = decomposition.transpose();
        }

        storage.decompositions[index] = decomposition;
    }

    @Override
    public void updatePopulationSizes(int index, double[] sizes, boolean flip) {
        assert sizes.length == stateCount;

        System.arraycopy(sizes, 0, storage.sizes, index * stateCount, stateCount);
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

    private static void reduceAcrossIntervals(double[] e, double[] f, double[] g, double[] h,
                                                int interval, double length,
                                                double[] sizes, double[] coalescent,
                                                double[] out,
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

        out[0] += logL;
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


package dr.evomodel.coalescent.basta;

import beagle.*;
import beagle.basta.BeagleBasta;
import beagle.basta.BastaFactory;
import dr.evolution.tree.Tree;
import dr.evomodel.substmodel.EigenDecomposition;
import dr.evomodel.treedatalikelihood.BufferIndexHelper;
import dr.inference.model.CompoundLikelihood;
import dr.inference.model.Likelihood;
import dr.math.matrixAlgebra.WrappedVector;

import java.util.Arrays;
import java.util.List;
import java.util.logging.Logger;

import static beagle.basta.BeagleBasta.BASTA_OPERATION_SIZE;
import static dr.evomodel.treedatalikelihood.BeagleFunctionality.parseSystemPropertyIntegerArray;

/**
 * @author Marc A. Suchard
 */
public class BeagleBastaLikelihoodDelegate extends BastaLikelihoodDelegate.AbstractBastaLikelihoodDelegate {

    private static final int COALESCENT_PROBABILITY_INDEX = 0;

    private final BeagleBasta beagle;

    private final BufferIndexHelper eigenBufferHelper;
    private final OffsetBufferIndexHelper populationSizesBufferHelper;
    private static final String RESOURCE_ORDER_PROPERTY = "beagle.resource.order";
    private static final String PREFERRED_FLAGS_PROPERTY = "beagle.preferred.flags";
    private static final String REQUIRED_FLAGS_PROPERTY = "beagle.required.flags";
    int currentPartialsCount;
    int currentIntervalsCount;
    private static int instanceCount = 0;
    private static List<Integer> resourceOrder = null;
    private static List<Integer> preferredOrder = null;
    private static List<Integer> requiredOrder = null;
    private int currentOutputBuffer;
    private int maxOutputBuffer;
    private boolean updateStorage;
    private int currentTransitionMatrixBuffer = 0;
    private int storedTransitionMatrixBuffer = 0;

    public BeagleBastaLikelihoodDelegate(String name,
                                         Tree tree,
                                         int stateCount,
                                         boolean transpose) {
        super(name, tree, stateCount, transpose);

        this.currentPartialsCount = 3 * tree.getNodeCount();
        this.currentIntervalsCount = tree.getNodeCount();

        int coalescentBufferCount = 5; // E, F, G, H, probabilities
        if (resourceOrder == null) {
            resourceOrder = parseSystemPropertyIntegerArray(RESOURCE_ORDER_PROPERTY);
        }
        if (preferredOrder == null) {
            preferredOrder = parseSystemPropertyIntegerArray(PREFERRED_FLAGS_PROPERTY);
        }
        if (requiredOrder == null) {
            requiredOrder = parseSystemPropertyIntegerArray(REQUIRED_FLAGS_PROPERTY);
        }

        long requirementFlags = 0L;
        requirementFlags |= BeagleFlag.EIGEN_COMPLEX.getMask();

        int[] resourceList = null;
        long preferenceFlags = 0;

        if (resourceOrder.size() > 0) {
            // added the zero on the end so that a CPU is selected if requested resource fails
            resourceList = new int[]{resourceOrder.get(instanceCount % resourceOrder.size()), 0};
            if (resourceList[0] > 0) {
                preferenceFlags |= BeagleFlag.PROCESSOR_GPU.getMask(); // Add preference weight against CPU
            }
        }

        if (preferredOrder.size() > 0) {
            preferenceFlags = preferredOrder.get(instanceCount % preferredOrder.size());
        }

        if (requiredOrder.size() > 0) {
            requirementFlags = requiredOrder.get(instanceCount % requiredOrder.size());
        }

        if (!BeagleFlag.PRECISION_SINGLE.isSet(preferenceFlags)) {
            // if single precision not explicitly set then prefer double
            preferenceFlags |= BeagleFlag.PRECISION_DOUBLE.getMask();
        }

        if ((resourceList == null &&
                (BeagleFlag.PROCESSOR_GPU.isSet(preferenceFlags) ||
                        BeagleFlag.FRAMEWORK_CUDA.isSet(preferenceFlags) ||
                        BeagleFlag.FRAMEWORK_OPENCL.isSet(preferenceFlags)))
                ||
                (resourceList != null && resourceList[0] > 0)) {
            // non-CPU implementations don't have SSE so remove default preference for SSE
            // when using non-CPU preferences or prioritising non-CPU resource
            preferenceFlags &= ~BeagleFlag.VECTOR_SSE.getMask();
            preferenceFlags &= ~BeagleFlag.THREADING_CPP.getMask();
        }

        beagle = BastaFactory.loadBastaInstance(0, coalescentBufferCount, maxNumCoalescentIntervals,
                currentPartialsCount, 0, stateCount,
                1, 2, 2 * currentIntervalsCount, 1,
                1, resourceList, preferenceFlags, requirementFlags);

        eigenBufferHelper = new BufferIndexHelper(1, 0);
        populationSizesBufferHelper = new OffsetBufferIndexHelper(1, 0, 0);

        beagle.setCategoryRates(new double[] { 1.0 });

        final Logger logger = Logger.getLogger("dr.evomodel");

        logger.info("\nCreating BeagleBastaLikelihoodDelegate");

        InstanceDetails instanceDetails = beagle.getDetails();
        ResourceDetails resourceDetails = null;

        if (instanceDetails != null) {
            resourceDetails = BeagleFactory.getResourceDetails(instanceDetails.getResourceNumber());
            if (resourceDetails != null) {
                StringBuilder sb = new StringBuilder("  Using BEAGLE BASTA resource ");
                sb.append(resourceDetails.getNumber()).append(": ");
                sb.append(resourceDetails.getName()).append("\n");
                if (resourceDetails.getDescription() != null) {
                    String[] description = resourceDetails.getDescription().split("\\|");
                    for (String desc : description) {
                        if (desc.trim().length() > 0) {
                            sb.append("    ").append(desc.trim()).append("\n");
                        }
                    }
                }
                sb.append("    with instance flags: ").append(instanceDetails.toString());
                logger.info(sb.toString());
            } else {
                logger.info("  Error retrieving BEAGLE resource for instance: " + instanceDetails.toString());
            }
        } else {
            logger.info("  No external BEAGLE resources available, or resource list/requirements not met, using Java implementation");
        }
    }


    public void resize(int newNumPartials, int newNumCoalescentIntervals) {
        updateStorage = false;
        if (newNumPartials > currentPartialsCount) {
            this.currentPartialsCount = newNumPartials + 1;
            updateStorage = true;
        }

        if (newNumCoalescentIntervals > currentIntervalsCount) {
            this.currentIntervalsCount = newNumCoalescentIntervals;
            updateStorage = true;
        }

        if (updateStorage) {
            beagle.allocateCoalescentBuffers(5, currentIntervalsCount, currentPartialsCount, 0);
        }
    }

    @Override
    protected void allocateGradientMemory() {
        // Do nothing
    }

    @Override
    protected void computeBranchIntervalOperations(List<Integer> intervalStarts,
                                                   List<BranchIntervalOperation> branchIntervalOperations,
                                                   List<TransitionMatrixOperation> matrixOperations,
                                                   Mode mode, BastaLikelihood likelihood) {

        int[] operations = new int[branchIntervalOperations.size() * BASTA_OPERATION_SIZE]; // TODO instantiate once
        int[] intervals = new int[intervalStarts.size()]; // TODO instantiate once
        double[] lengths = new double[intervalStarts.size() - 1]; // TODO instantiate once

        vectorizeBranchIntervalOperations(intervalStarts, branchIntervalOperations, operations, intervals, lengths);
        updateStorage(maxOutputBuffer, maxNumCoalescentIntervals, likelihood);
        int populationSizeIndex = populationSizesBufferHelper.getOffsetIndex(0);

        if (mode == Mode.LIKELIHOOD) {
            beagle.updateBastaPartials(operations, branchIntervalOperations.size(),
                    intervals, intervalStarts.size(),
                    populationSizeIndex, COALESCENT_PROBABILITY_INDEX);
            // TODO do dipatch inside BEAGLE by passing mode.getModeAsInt()
        } else if (mode == Mode.GRADIENT) {
            beagle.updateBastaPartialsGrad(operations, branchIntervalOperations.size(),
                    intervals, intervalStarts.size(),
                    populationSizeIndex, COALESCENT_PROBABILITY_INDEX);
        } else {
            throw new RuntimeException("not yet implemented");
        }

        if (PRINT_COMMANDS) {
            double[] partials = new double[stateCount];
            for (BranchIntervalOperation operation : branchIntervalOperations) {
                getPartials(operation.outputBuffer, partials);
                System.err.println(operation + " " + new WrappedVector.Raw(partials));
            }

            double[] coalescent = new double[maxNumCoalescentIntervals];
            beagle.getBastaBuffer(COALESCENT_PROBABILITY_INDEX, coalescent);
            System.err.println("coalescent: " + new WrappedVector.Raw(coalescent));
        }
    }

    @Override
    String getStamp() { return "beagle"; }

    public Beagle getBeagleInstance() { return beagle; }

    @Override
    protected void computeTransitionProbabilityOperations(List<TransitionMatrixOperation> matrixOperations,
                                                          Mode mode) {

        int[] transitionMatrixIndices = new int[matrixOperations.size()]; // TODO instantiate once
        double[] branchLengths = new double[matrixOperations.size()]; // TODO instantiate once

        vectorizeTransitionMatrixOperations(matrixOperations, transitionMatrixIndices, branchLengths);

        int eigenIndex = eigenBufferHelper.getOffsetIndex(0);

        if (mode == Mode.LIKELIHOOD) {
            beagle.updateTransitionMatrices(eigenIndex, transitionMatrixIndices, null, null,
                    branchLengths, matrixOperations.size());
        } else if (mode == Mode.GRADIENT) {
            beagle.updateTransitionMatricesGrad(transitionMatrixIndices,
                    branchLengths, matrixOperations.size());
        }

        if (PRINT_COMMANDS) {
            double[] matrix = new double[stateCount * stateCount];
            for (TransitionMatrixOperation operation : matrixOperations) {
                getTransitionMatrices(operation.outputBuffer, matrix);
                System.err.println(operation + " " + new WrappedVector.Raw(matrix));
            }
        }
    }

    @Override
    protected void computeCoalescentIntervalReduction(List<Integer> intervalStarts,
                                                        List<BranchIntervalOperation> branchIntervalOperations,
                                                        double[] out,
                                                        Mode mode,
                                                        StructuredCoalescentLikelihoodGradient.WrtParameter wrt) {

        // TODO vectorize once per likelihood-calculation
        int[] operations = new int[branchIntervalOperations.size() * BASTA_OPERATION_SIZE]; // TODO instantiate once
        int[] intervals = new int[intervalStarts.size()]; // TODO instantiate once
        double[] lengths = new double[intervalStarts.size() - 1]; // TODO instantiate once

        vectorizeBranchIntervalOperations(intervalStarts, branchIntervalOperations, operations, intervals, lengths);

        int populationSizeIndex = populationSizesBufferHelper.getOffsetIndex(0);

        beagle.accumulateBastaPartials(operations, branchIntervalOperations.size(),
                intervals, intervalStarts.size(), lengths,
                populationSizeIndex, COALESCENT_PROBABILITY_INDEX,
                out);
    }

//    @Override
    protected double[] computeCoalescentIntervalReductionGrad(List<Integer> intervalStarts,
                                                              List<BranchIntervalOperation> branchIntervalOperations,
                                                              double[] gradient) {
        int[] operations = new int[branchIntervalOperations.size() * BASTA_OPERATION_SIZE]; // TODO instantiate once
        int[] intervals = new int[intervalStarts.size()]; // TODO instantiate once
        double[] lengths = new double[intervalStarts.size() - 1]; // TODO instantiate once

        vectorizeBranchIntervalOperations(intervalStarts, branchIntervalOperations, operations, intervals, lengths);

        int populationSizeIndex = populationSizesBufferHelper.getOffsetIndex(0);

        double[] out = new double[stateCount*stateCount];
        beagle.accumulateBastaPartialsGrad(operations, branchIntervalOperations.size(),
                intervals, intervalStarts.size(), lengths,
                populationSizeIndex, COALESCENT_PROBABILITY_INDEX,
                out);

        return out;
    }

//    @Override
    protected double[] computeCoalescentIntervalReductionGradPopSize(List<Integer> intervalStarts,
                                                                     List<BranchIntervalOperation> branchIntervalOperations,
                                                                     double[] out) {
        return new double[0];
    }

    @Override
    public void setPartials(int index, double[] partials) {
        beagle.setPartials(index, partials);
    }

    @Override
    public void getPartials(int index, double[] partials) {
        assert index >= 0;
        assert partials != null;
        assert partials.length >= stateCount;

        beagle.getPartials(index, Beagle.NONE, partials);
    }

    @Override
    public void getTransitionMatrices(int index, double[] matrix) {
        assert index >= 0;
        assert matrix != null;
        assert matrix.length >= stateCount * stateCount;

        beagle.getTransitionMatrix(index, matrix);
    }

    @Override
    public void updateEigenDecomposition(int index, EigenDecomposition decomposition, boolean flip) {
        if (flip) {
            eigenBufferHelper.flipOffset(0);
        }

        if (transpose) {
            decomposition = decomposition.transpose();
        }

        beagle.setEigenDecomposition(
                eigenBufferHelper.getOffsetIndex(0),
                decomposition.getEigenVectors(),
                decomposition.getInverseEigenVectors(),
                decomposition.getEigenValues());
    }

    @Override
    public void updatePopulationSizes(int index, double[] sizes, boolean flip) {
        if (flip) {
            populationSizesBufferHelper.flipOffset(0);
        }

        beagle.setStateFrequencies(populationSizesBufferHelper.getOffsetIndex(0), sizes);
    }

    @Override
    public void updateStorage(int maxBufferCount, int treeNodeCount, BastaLikelihood likelihood) {
        int newNumPartials = maxBufferCount + 1;
        resize(newNumPartials, maxNumCoalescentIntervals);
        if (likelihood != null && updateStorage) {
            likelihood.setTipData();
        }
    }

    @Override
    public void storeState() {
        //populationSizesBufferHelper.storeState();
        //eigenBufferHelper.storeState();
        storedTransitionMatrixBuffer = currentTransitionMatrixBuffer;
    }

    @Override
    public void restoreState() {
        //populationSizesBufferHelper.restoreState();
       // eigenBufferHelper.restoreState();
        currentTransitionMatrixBuffer = storedTransitionMatrixBuffer;
    }


    public void flipTransitionMatrixBuffer(List<TransitionMatrixOperation> matrixOperations) {
        currentTransitionMatrixBuffer = (currentTransitionMatrixBuffer == 0 ? matrixOperations.size() : 0);
    }

    private void vectorizeTransitionMatrixOperations(List<TransitionMatrixOperation> matrixOperations,
                                                     int[] transitionMatrixIndices,
                                                     double[] branchLengths) {
        int k = 0;
        int currentTMBuffer = currentTransitionMatrixBuffer;
        for (TransitionMatrixOperation op : matrixOperations) {
            transitionMatrixIndices[k] = op.outputBuffer + currentTMBuffer; // TODO double-buffer
            branchLengths[k] = op.time;
            ++k;
        }
    }

    int tipCount = -1;
    int[] map;
    int used;

    private static final boolean CACHE_FRIENDLY = true;

    int map(int buffer) {
        if (buffer < tipCount) {
            return buffer;
        } else {
            if (map[buffer] == -1) {
                map[buffer] = used;
                ++used;
            }
            return map[buffer];
        }
    }


   public void vectorizeBranchIntervalOperations(List<Integer> intervalStarts,
                                                   List<BranchIntervalOperation> branchIntervalOperations,
                                                   int[] operations,
                                                   int[] intervals,
                                                   double[] lengths) {

        if (CACHE_FRIENDLY) {
            this.tipCount = tree.getExternalNodeCount();
            if (this.map == null) {
                this.map = new int[maxNumCoalescentIntervals * (tree.getNodeCount() + 1)];
            }
            Arrays.fill(this.map, -1);
            this.used = this.tipCount;
        }

        // TODO double-buffer
        int k = 0;
        int currentTMBuffer = currentTransitionMatrixBuffer;
        for (BranchIntervalOperation op : branchIntervalOperations) {

            if (CACHE_FRIENDLY) {
                operations[k] = map(op.outputBuffer);
                operations[k + 1] = map(op.inputBuffer1);
                operations[k + 2] = op.inputMatrix1 + currentTMBuffer;
                operations[k + 3] = map(op.inputBuffer2);
                operations[k + 4] = op.inputMatrix2 + currentTMBuffer;
                operations[k + 5] = map(op.accBuffer1);
                operations[k + 6] = map(op.accBuffer2);
                operations[k + 7] = op.intervalNumber;
            } else {
                operations[k] = op.outputBuffer;
                operations[k + 1] = op.inputBuffer1;
                operations[k + 2] = op.inputMatrix1 + currentTMBuffer;
                operations[k + 3] = op.inputBuffer2;
                operations[k + 4] = op.inputMatrix2 + currentTMBuffer;
                operations[k + 5] = op.accBuffer1;
                operations[k + 6] = op.accBuffer2;
                operations[k + 7] = op.intervalNumber;
            }
            currentOutputBuffer = operations[k + 6];

            if (currentOutputBuffer > maxOutputBuffer) {
                maxOutputBuffer = currentOutputBuffer;
            }

            k += BASTA_OPERATION_SIZE;
        }
        
        int i = 0;
        for (int end = intervalStarts.size() - 1; i < end; ++i) {
            int start = intervalStarts.get(i);
            intervals[i] = start;
            lengths[i] = branchIntervalOperations.get(start).intervalLength;
        }
        intervals[i] = intervalStarts.get(i);
    }

    private boolean releaseSingleton = true;

    private void releaseBeagle() throws Throwable {
        if (beagle != null && releaseSingleton) {
            beagle.finalize();
            releaseSingleton = false;
        }
    }

    public static void releaseBeagleBastaLikelihoodDelegate(BastaLikelihood treeDataLikelihood) throws Throwable {
        BastaLikelihoodDelegate likelihoodDelegate = treeDataLikelihood.getLikelihoodDelegate();
        if (likelihoodDelegate instanceof BeagleBastaLikelihoodDelegate) {
            BeagleBastaLikelihoodDelegate delegate = (BeagleBastaLikelihoodDelegate) likelihoodDelegate;
            delegate.releaseBeagle();
        }
    }

    public static void releaseAllBeagleBastaInstances() throws Throwable {
        for (Likelihood likelihood : dr.inference.model.Likelihood.FULL_LIKELIHOOD_SET) {
            if (likelihood instanceof BastaLikelihood) {
                releaseBeagleBastaLikelihoodDelegate((BastaLikelihood) likelihood);
            } else if (likelihood instanceof CompoundLikelihood) {
                for (Likelihood likelihood2: ((CompoundLikelihood) likelihood).getLikelihoods()) {
                    if (likelihood2 instanceof BastaLikelihood) {
                        releaseBeagleBastaLikelihoodDelegate((BastaLikelihood) likelihood2);
                    }
                }
            }
        }
    }

    static class OffsetBufferIndexHelper extends BufferIndexHelper {

        public OffsetBufferIndexHelper(int maxIndexValue, int minIndexValue, int bufferSetNumber) {
            super(maxIndexValue, minIndexValue, bufferSetNumber);
        }

        @Override
        protected int computeOffset(int offset) { return offset; }
    }
}

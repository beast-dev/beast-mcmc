package dr.evomodel.coalescent.basta;

import beagle.Beagle;
import beagle.basta.BeagleBasta;
import beagle.basta.BastaFactory;
import dr.evolution.tree.Tree;
import dr.evomodel.substmodel.EigenDecomposition;
import dr.evomodel.treedatalikelihood.BufferIndexHelper;

import java.util.List;

/**
 * @author Marc A. Suchard
 */
public class BeagleBastaLikelihoodDelegate extends BastaLikelihoodDelegate.AbstractBastaLikelihoodDelegate {

    private final BeagleBasta beagle;

    private final BufferIndexHelper eigenBufferHelper;
    private final OffsetBufferIndexHelper populationSizesBufferHelper;

    public BeagleBastaLikelihoodDelegate(String name,
                                         Tree tree,
                                         int stateCount) {
        super(name, tree, stateCount);

        beagle = BastaFactory.loadBastaInstance(1, 1, 1, 16,
                1, 1, 1, 1,
                1, null, 0L, 0L);

        eigenBufferHelper = new BufferIndexHelper(1, 0);
        populationSizesBufferHelper = new OffsetBufferIndexHelper(1, 0, 0);

        double[] tmp = new double[16];
        beagle.setPartials(0, tmp);
    }

    @Override
    protected void computeBranchIntervalOperations(List<Integer> intervalStarts,
                                                   List<BranchIntervalOperation> branchIntervalOperations) {

    }

    @Override
    protected void computeTransitionProbabilityOperations(List<TransitionMatrixOperation> matrixOperations) {

    }

    @Override
    protected double computeCoalescentIntervalReduction(List<Integer> intervalStarts,
                                                        List<BranchIntervalOperation> branchIntervalOperations) {
        return 0;
    }

    @Override
    public void setPartials(int index, double[] partials) {
        beagle.setPartials(index, partials);
    }

    @Override
    public void getPartials(int index, double[] partials) {
        assert index >= 0;
        assert partials != null;

        beagle.getPartials(index, Beagle.NONE, partials);
    }

    @Override
    public void updateEigenDecomposition(int index, EigenDecomposition decomposition, boolean flip) {
        if (flip) {
            eigenBufferHelper.flipOffset(0);
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

        beagle.setPartials(populationSizesBufferHelper.getOffsetIndex(0),
                sizes);
    }

    class OffsetBufferIndexHelper extends BufferIndexHelper {

        public OffsetBufferIndexHelper(int maxIndexValue, int minIndexValue, int bufferSetNumber) {
            super(maxIndexValue, minIndexValue, bufferSetNumber);
        }

        @Override
        protected int computeOffset(int offset) { return offset; }
    }
}

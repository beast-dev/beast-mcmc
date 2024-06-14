package dr.evomodel.treedatalikelihood.action;

import beagle.Beagle;
import dr.evolution.tree.Tree;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.treedatalikelihood.EvolutionaryProcessDelegate;
import dr.evomodel.treedatalikelihood.PreOrderSettings;

import java.util.ArrayList;
import java.util.HashMap;

public class ActionSubstitutionModelDelegate implements EvolutionaryProcessDelegate {

    private final Tree tree;
    private final BranchModel branchModel;
    private final int nodeCount;

    private final int stateCount;

    private final HashMap<SubstitutionModel, Integer> eigenIndexMap;

    public ActionSubstitutionModelDelegate (Tree tree,
                                            BranchModel branchModel,
                                            int nodeCount) {
        this.tree = tree;
        this.branchModel = branchModel;
        this.nodeCount = nodeCount;
        this.stateCount = branchModel.getRootFrequencyModel().getFrequencyCount();
        this.eigenIndexMap = new HashMap<>();
        for (int i = 0; i < getSubstitutionModelCount(); i++) {
            eigenIndexMap.put(branchModel.getSubstitutionModels().get(i), i);
        }
    }

    @Override
    public boolean canReturnComplexDiagonalization() {
        return false;
    }

    @Override
    public int getEigenBufferCount() {
        return branchModel.getSubstitutionModels().size();
    }

    @Override
    public int getMatrixBufferCount() {
        return nodeCount - 1;
    }

    @Override
    public int getInfinitesimalMatrixBufferIndex(int branchIndex) {
        return getMatrixBufferCount() + eigenIndexMap.get(getSubstitutionModelForBranch(branchIndex));
    }

    @Override
    public int getInfinitesimalSquaredMatrixBufferIndex(int branchIndex) {
        return getMatrixBufferCount() + getEigenBufferCount() + getInfinitesimalMatrixBufferIndex(branchIndex);
    }

    @Override
    public int getFirstOrderDifferentialMatrixBufferIndex(int branchIndex) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public int getSecondOrderDifferentialMatrixBufferIndex(int branchIndex) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public void cacheInfinitesimalMatrix(Beagle beagle, int bufferIndex, double[] differentialMatrix) {
        beagle.setDifferentialMatrix(getInfinitesimalMatrixBufferIndex(bufferIndex), differentialMatrix);
    }

    @Override
    public void cacheInfinitesimalSquaredMatrix(Beagle beagle, int bufferIndex, double[] differentialMatrix) {
        beagle.setDifferentialMatrix(getInfinitesimalSquaredMatrixBufferIndex(bufferIndex), differentialMatrix);
    }

    @Override
    public void cacheFirstOrderDifferentialMatrix(Beagle beagle, int branchIndex, double[] differentialMassMatrix) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public int getCachedMatrixBufferCount(PreOrderSettings settings) {
        int matrixBufferCount = getInfinitesimalMatrixBufferCount(settings) + getDifferentialMassMatrixBufferCount(settings);
        return matrixBufferCount;
    }

    private int getInfinitesimalMatrixBufferCount(PreOrderSettings settings) {
        if (settings.branchRateDerivative) {
            return 2 * getEigenBufferCount();
        } else {
            return 0;
        }
    }

    private int getDifferentialMassMatrixBufferCount(PreOrderSettings settings) {
        if (settings.branchInfinitesimalDerivative) {
            return 2 * (nodeCount - 1);
        } else {
            return 0;
        }
    }
    @Override
    public int getSubstitutionModelCount() {
        return branchModel.getSubstitutionModels().size();
    }

    @Override
    public SubstitutionModel getSubstitutionModel(int index) {
        return branchModel.getSubstitutionModels().get(index);
    }

    @Override
    public SubstitutionModel getSubstitutionModelForBranch(int branchIndex) {
        BranchModel.Mapping mapping = branchModel.getBranchModelMapping(tree.getNode(branchIndex));
        int[] order = mapping.getOrder();

        if (order.length > 1) {
            throw new RuntimeException("Not yet implemented");
        }

        return getSubstitutionModel(order[0]);
    }

    @Override
    public int getEigenIndex(int bufferIndex) {
        throw new RuntimeException("Eigen index should not be called for actions.");
    }

    @Override
    public int getMatrixIndex(int branchIndex) {
        return branchIndex;
    }

    @Override
    public double[] getRootStateFrequencies() {
        return branchModel.getRootFrequencyModel().getFrequencies();
    }

    @Override
    public void updateSubstitutionModels(Beagle beagle, boolean flipBuffers) {

        double[] tmpQ = new double[stateCount * stateCount];
        for (int i = 0; i < getSubstitutionModelCount(); i++) {
            SubstitutionModel substitutionModel = getSubstitutionModel(i);
            substitutionModel.getInfinitesimalMatrix(tmpQ);

            ArrayList rowIndices = new ArrayList<Integer>();
            ArrayList colIndices = new ArrayList<Integer>();
            ArrayList values = new ArrayList<Double>();
            int nonZeros = 0;
            for (int row = 0; row < stateCount; row++) {
                for (int col = 0; col < stateCount; col++) {
                    if (tmpQ[row * stateCount + col] != 0) {
                        rowIndices.add(row);
                        colIndices.add(col);
                        values.add(tmpQ[row * stateCount + col]);
                        nonZeros++;
                    }
                }
            }
            int[] inRowIndices = new int[rowIndices.size()];
            int[] inColIndices = new int[colIndices.size()];
            double[] inValues = new double[values.size()];
            for (int j = 0; j < nonZeros; j++) {
                inRowIndices[j] = (int) rowIndices.get(j);
                inColIndices[j] = (int) colIndices.get(j);
                inValues[j] = (double) values.get(j);
            }

            beagle.setSparseMatrix(i, inRowIndices, inColIndices, inValues, nonZeros);
        }

    }

    private int getEigenIndexForBranch(int branchIndex) {
        SubstitutionModel substitutionModel = getSubstitutionModelForBranch(branchIndex);
        return eigenIndexMap.get(substitutionModel);
    }

    @Override
    public void updateTransitionMatrices(Beagle beagle, int[] branchIndices, double[] edgeLengths, int updateCount, boolean flipBuffers) {
        for (int i = 0; i < updateCount; i++) {
            beagle.updateTransitionMatrices(getEigenIndexForBranch(branchIndices[i]), new int[]{branchIndices[i]}, null, null, new double[]{edgeLengths[i]}, 1);
        }
    }

    @Override
    public void flipTransitionMatrices(int[] branchIndices, int updateCount) {
        throw new RuntimeException("Transition matrix is not used for actions.");
    }

    @Override
    public void storeState() {

    }

    @Override
    public void restoreState() {

    }
}

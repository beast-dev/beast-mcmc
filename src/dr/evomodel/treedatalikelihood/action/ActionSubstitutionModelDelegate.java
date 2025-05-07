package dr.evomodel.treedatalikelihood.action;

import beagle.Beagle;
import dr.evolution.tree.Tree;
import dr.evomodel.branchmodel.BranchModel;
import dr.evomodel.substmodel.ActionEnabledSubstitution;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.treedatalikelihood.EvolutionaryProcessDelegate;
import dr.evomodel.treedatalikelihood.PreOrderSettings;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class ActionSubstitutionModelDelegate implements EvolutionaryProcessDelegate {

    private final Tree tree;
    private final BranchModel branchModel;
    private final int nodeCount;

    private final List<ActionEnabledSubstitution> substitutionModels;

    private final int stateCount;

    private final HashMap<SubstitutionModel, Integer> eigenIndexMap;

    public ActionSubstitutionModelDelegate (Tree tree,
                                            BranchModel branchModel,
                                            int nodeCount) {
        this.tree = tree;
        this.branchModel = branchModel;
        this.substitutionModels = getSubstitutionModels(branchModel);
        this.nodeCount = nodeCount;
        this.stateCount = branchModel.getRootFrequencyModel().getFrequencyCount();
        this.eigenIndexMap = new HashMap<>();
        for (int i = 0; i < getSubstitutionModelCount(); i++) {
            eigenIndexMap.put(substitutionModels.get(i), i);
        }
    }

    private List<ActionEnabledSubstitution> getSubstitutionModels(BranchModel branchModel) {
        List<ActionEnabledSubstitution> substitutionModels = new ArrayList<>();
        for (SubstitutionModel substitutionModel : branchModel.getSubstitutionModels()) {
            if (substitutionModel instanceof ActionEnabledSubstitution) {
                substitutionModels.add((ActionEnabledSubstitution) substitutionModel);
            } else {
                substitutionModels.add(new ActionEnabledSubstitution.ActionEnabledSubstitutionWrap("original.substitution.model", substitutionModel));
            }
        }
        return substitutionModels;
    }

    @Override
    public boolean canReturnComplexDiagonalization() {
        return false;
    }

    @Override
    public int getEigenBufferCount() {
        return substitutionModels.size();
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
        final int offset = getMatrixBufferCount();
        int bufferIndex = offset + branchIndex;
        return bufferIndex;
    }

    @Override
    public int getSecondOrderDifferentialMatrixBufferIndex(int branchIndex) {
        return getFirstOrderDifferentialMatrixBufferIndex(branchIndex) + nodeCount - 1;
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
        beagle.setDifferentialMatrix(getFirstOrderDifferentialMatrixBufferIndex(branchIndex), differentialMassMatrix);
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
        return substitutionModels.size();
    }

    @Override
    public SubstitutionModel getSubstitutionModel(int index) {
        return substitutionModels.get(index);
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

        for (int i = 0; i < getSubstitutionModelCount(); i++) {
            ActionEnabledSubstitution substitutionModel = (ActionEnabledSubstitution) getSubstitutionModel(i);

            int[] rowIndices = new int[substitutionModel.getNonZeroEntryCount()];
            int[] colIndices = new int[substitutionModel.getNonZeroEntryCount()];
            double[] values = new double[substitutionModel.getNonZeroEntryCount()];

            substitutionModel.getNonZeroEntries(rowIndices, colIndices, values);
            beagle.setSparseMatrix(i, rowIndices, colIndices, values, substitutionModel.getNonZeroEntryCount());
        }

    }

    public void cacheSparseDifferentialMatrix(Beagle beagle, int substitutionIndex) {
        ActionEnabledSubstitution substitutionModel = (ActionEnabledSubstitution) getSubstitutionModel(substitutionIndex);

        int[] rowIndices = new int[substitutionModel.getNonZeroEntryCount()];
        int[] colIndices = new int[substitutionModel.getNonZeroEntryCount()];
        double[] values = new double[substitutionModel.getNonZeroEntryCount()];

        substitutionModel.getNonZeroEntries(rowIndices, colIndices, values);
        beagle.setSparseDifferentialMatrix(getInfinitesimalMatrixBufferIndex(substitutionIndex), rowIndices, colIndices, values, substitutionModel.getNonZeroEntryCount());
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

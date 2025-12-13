package dr.evomodel.continuous;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeAttributeProvider;
import dr.inference.model.*;
import org.ejml.data.DenseMatrix64F;

import static dr.evomodel.treedatalikelihood.hmc.AbstractPrecisionGradient.flatten;

/**
 * Multivariate OU "elastic" model:
 *     dX_t = -A (X_t - μ) dt + Σ dW_t
 *
 * Encodes A (the strength-of-selection matrix) via a MatrixParameterInterface.
 * Supports multiple parametrizations through a strategy pattern:
 *   - Diagonal: A = diag(λ_i)
 *   - Decomposed: A = R Λ R^{-1} from CompoundEigenMatrix
 *   - General: A = R Λ R^{-1} via EJML eigen-decomposition
 *   - Block Diagonal: A = R D R^{-1} with block structure
 *
 * @author Marc A. Suchard
 * @author Paul Bastide
 */
public class MultivariateElasticModel extends AbstractModel implements TreeAttributeProvider {

    // ============================================================
    // CONSTANTS
    // ============================================================

    private static final String ELASTIC_PROCESS = "multivariateElasticModel";
    private static final String ELASTIC_TREE_ATTRIBUTE = "strengthOfSelection";
    public static final double LOG2PI = Math.log(2 * Math.PI);

    // ============================================================
    // CORE STATE
    // ============================================================

    private final int dim;
    private final MatrixParameterInterface strengthOfSelectionParam;
    private final SelectionMatrixStrategy strategy;

    // Configuration flag
    private boolean useApproximateExponentialGradient = false;

    private BasisRepresentation currentBasis;
    private BasisRepresentation savedBasis;
    private boolean basisKnown = false;

    // ============================================================
    // CONSTRUCTION
    // ============================================================

    /**
     * Main constructor.
     */
    public MultivariateElasticModel(MatrixParameterInterface strengthOfSelectionMatrixParameter) {
        super(ELASTIC_PROCESS);

        this.strengthOfSelectionParam = strengthOfSelectionMatrixParameter;
        this.dim = strengthOfSelectionMatrixParameter.getRowDimension();

        validateSquareMatrix(strengthOfSelectionMatrixParameter, dim);

        this.strategy = SelectionMatrixStrategyFactory.create(strengthOfSelectionMatrixParameter);

        addVariable(strengthOfSelectionMatrixParameter);
    }

    /**
     * Empty constructor for XML / BEAST-style instantiation.
     */
    public MultivariateElasticModel() {
        super(ELASTIC_PROCESS);
        this.dim = 0;
        this.strengthOfSelectionParam = null;
        this.strategy = null;
    }

    private void validateSquareMatrix(MatrixParameterInterface param, int dim) {
        if (dim != param.getColumnDimension()) {
            throw new IllegalArgumentException(
                    "Strength of Selection matrix must be square. Got " +
                            dim + "x" + param.getColumnDimension());
        }
    }

    // ============================================================
    // PUBLIC API - Matrix Access
    // ============================================================

    public MatrixParameterInterface getMatrixParameterInterface() {
        return strengthOfSelectionParam;
    }

    public MatrixParameterInterface getStrengthOfSelectionMatrixParameter() {
        return strengthOfSelectionParam;
    }

    public double[][] getStrengthOfSelectionMatrix() {
        if (strengthOfSelectionParam == null) {
            return null;
        }
        return strengthOfSelectionParam.getParameterAsMatrix();
    }

    public double[] getStrengthOfSelectionMatrixAsVector() {
        return flatten(getStrengthOfSelectionMatrix());
    }

    public int getDimension() {
        return dim;
    }

    public int getNumberOfParameters() {
        return strategy.getNumberOfParameters();
    }

    // ============================================================
    // PUBLIC API - Basis Representation
    // ============================================================

    /**
     * Eigenvalues in the working basis.
     * For block-diagonal, returns the diagonal elements of D.
     */
    public double[] getBasisEigenValues() {
        return getBasis().getEigenvalues();
    }

    /**
     * Flattened basis matrix R (row-major).
     */
    public double[] getBasisEigenVectors() {
        return getBasis().getEigenvectors();
    }

    public double[] getBasisRotations() {
        return getBasis().getRotations();
    }

    public double[] getBasisD() {
        return getBasis().getValuesD();
    }

    /**
     * Basis matrix R used to represent A in (quasi-)diagonal form.
     */
    public DenseMatrix64F getBasisR() {
        return getBasis().getR();
    }

    /**
     * Inverse of basis matrix R.
     */
    public DenseMatrix64F getBasisRinv() {
        return getBasis().getRinv();
    }

    // ============================================================
    // PUBLIC API - Block-Diagonal Specific
    // ============================================================

    /**
     * Starting indices of each block in the block-diagonal structure.
     * Only available for block-diagonal parametrization.
     *
     * @throws IllegalStateException if not block-diagonal
     */
    public int[] getBlockStarts() {
        return getBlockStructure().getBlockStarts();
    }

    /**
     * Sizes (1 or 2) of each block in the block-diagonal structure.
     * Only available for block-diagonal parametrization.
     *
     * @throws IllegalStateException if not block-diagonal
     */
    public int[] getBlockSizes() {
        return getBlockStructure().getBlockSizes();
    }

    private BlockStructure getBlockStructure() {
        BlockStructure blockStructure = getBasis().getBlockStructure();
        if (blockStructure == null) {
            return null;
        }
        return blockStructure;
    }
    public boolean hasBlockStructure() {
        return getBlockStructure() != null;
    }

    // ============================================================
    // PUBLIC API - Type Queries
    // ============================================================

    public boolean isDiagonal() {
        return strategy != null && strategy.isDiagonal();
    }

    public boolean isSymmetric() {
        return strategy != null && strategy.isSymmetric();
    }

    public boolean isBlockDiag() {
        return strategy != null && strategy.isBlockDiagonal();
    }

    // ============================================================
    // PUBLIC API - Configuration
    // ============================================================

    /**
     * If true, delegates may use an approximate gradient of the matrix exponential
     * (e.g. d exp ≈ -t exp * J) instead of exact differentiation.
     */
    public void setApproximateExponentialGradient(boolean approximate) {
        this.useApproximateExponentialGradient = approximate;
    }

    public boolean isApproximateExponentialGradient() {
        return useApproximateExponentialGradient;
    }

    // ============================================================
    // INTERNAL - Basis Computation
    // ============================================================

    private BasisRepresentation getBasis() {
        if (!basisKnown) {
            currentBasis = strategy.computeBasis(strengthOfSelectionParam, dim);
            basisKnown = true;
        }
        return currentBasis;
    }

    // ============================================================
    // MODEL LIFECYCLE
    // ============================================================

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        //nothing
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        basisKnown = false;
    }

    @Override
    protected void storeState() {
        savedBasis = currentBasis;
    }

    @Override
    protected void restoreState() {
        currentBasis = savedBasis;
        basisKnown = (currentBasis != null);  // or just = true if we guarantee non-null
    }

    @Override
    protected void acceptState() {
        // Immutable basis representations don't need acceptance
    }

    // ============================================================
    // TREE ATTRIBUTES
    // ============================================================

    @Override
    public String[] getTreeAttributeLabel() {
        return new String[]{ELASTIC_TREE_ATTRIBUTE};
    }

    @Override
    public String[] getAttributeForTree(Tree tree) {
        if (strengthOfSelectionParam != null) {
            return new String[]{strengthOfSelectionParam.toString()};
        }
        return new String[]{"null"};
    }

    // ============================================================
    // INNER CLASSES - Strategy Pattern
    // ============================================================

    /**
     * Strategy interface for computing basis representations from different
     * matrix parametrizations.
     */
    interface SelectionMatrixStrategy {
        BasisRepresentation computeBasis(MatrixParameterInterface param, int dim);
        boolean isDiagonal();
        boolean isSymmetric();
        boolean isBlockDiagonal();
        default int getNumberOfParameters() {
            return -1; // Unknown
        }
    }

    /**
     * Factory for creating appropriate strategy based on parameter type.
     */
    static class SelectionMatrixStrategyFactory {
        static SelectionMatrixStrategy create(MatrixParameterInterface param) {
            if (param instanceof BlockDiagonalCosSinMatrixParameter) {
                return new BlockDiagonalStrategy();
            } else if (param instanceof DiagonalMatrix) {
                return new DiagonalStrategy();
            } else if (param instanceof CompoundEigenMatrix) {
                return new DecomposedStrategy(param.isConstrainedSymmetric());
            } else {
                return new GeneralEigenStrategy(param.isConstrainedSymmetric());
            }
        }
    }
}
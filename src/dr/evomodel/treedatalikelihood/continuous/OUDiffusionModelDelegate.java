package dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.continuous.MultivariateElasticModel;
import dr.evomodel.treedatalikelihood.continuous.cdi.ContinuousDiffusionIntegrator;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import dr.inference.model.Model;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;
import dr.evomodel.treedatalikelihood.continuous.cdi.*;

import java.util.List;

/**
 * OU diffusion model delegate with branch-specific drift and constant diffusion.
 *
 * This version offloads all actualization- and Fréchet-related logic to
 * pluggable OUActualizationStrategy implementations:
 *   - DiagonalOUActualizationStrategy
 *   - SymmetricOUActualizationStrategy
 *   - BlockSchurOUActualizationStrategy
 *   - GeneralOUActualizationStrategy
 *
 * Responsibilities of this class:
 *  - glue between tree / CDI and OU strategies,
 *  - model-change plumbing,
 *  - root / fixed-root semantics for displacement wrt root,
 *  - CDI wiring: stationary variance & OU transition matrices.
 */
/**
 * A simple OU diffusion model delegate with branch-specific drift and constant diffusion
 *
 * @author Marc A. Suchard
 * @author Paul Bastide
 */
public class OUDiffusionModelDelegate
        extends AbstractDriftDiffusionModelDelegate
        implements OUDelegateContext {

    // ----------------------------------------------------------------------
    // Core models
    // ----------------------------------------------------------------------

    private final MultivariateElasticModel   elasticModel;
    private final MultivariateDiffusionModel diffusionModel;

    // Strategy handling all OU actualization & adjoint math
    private final OUActualizationStrategy actualizationStrategy;

    // ----------------------------------------------------------------------
    // Constructors
    // ----------------------------------------------------------------------

    public OUDiffusionModelDelegate(Tree tree,
                                         MultivariateDiffusionModel diffusionModel,
                                         List<BranchRateModel> branchRateModels,
                                         MultivariateElasticModel elasticModel) {
        this(tree, diffusionModel, branchRateModels, elasticModel, 0);
    }

    private OUDiffusionModelDelegate(Tree tree,
                                          MultivariateDiffusionModel diffusionModel,
                                          List<BranchRateModel> branchRateModels,
                                          MultivariateElasticModel elasticModel,
                                          int partitionNumber) {
        super(tree, diffusionModel, branchRateModels, partitionNumber);
        this.elasticModel   = elasticModel;
        this.diffusionModel = diffusionModel;
        addModel(elasticModel);

        // ------------------------------------------------------------------
        // [Legacy] Strategy selection (all OU logic lives in these classes)
        // ------------------------------------------------------------------
        if (elasticModel.isDiagonal()) {
            this.actualizationStrategy = new OUActualizationStrategyDiagonal(this);
        } else if (elasticModel.isBlockDiag()) {
            this.actualizationStrategy = null; // Use backpropagation for the block-diagonal case
//                    new OUActualizationStrategyBlockDecomposition(this);
//            throw new RuntimeException("Use backpropagation for the block-diagonal case.");
        } else {
            this.actualizationStrategy = new OUHurwitzActualizationStrategy(this);
        }
    }

    // ----------------------------------------------------------------------
    // OUDelegateContext implementation
    // ----------------------------------------------------------------------

    @Override
    public Tree getTree() {
        return tree;
    }

    @Override
    public MultivariateElasticModel getElasticModel() {
        return elasticModel;
    }

    @Override
    public int getDim() {
        return dim;
    }

    public double[] getOptimalValues(NodeRef node) {
        return getDriftRate(node);
    }
    // In OUDiffusionModelDelegate

    public DenseMatrix64F getSelectionMatrixAsEJML() {
        int dim = getDim();
        double[][] Sraw = elasticModel.getStrengthOfSelectionMatrix(); // already exists
        DenseMatrix64F S = new DenseMatrix64F(dim, dim);
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                S.set(i, j, Sraw[i][j]);
            }
        }
        return S;
    }

    public DenseMatrix64F getDiffusionCovarianceAsEJML() {
        int dim = getDim();
        double[] temp = diffusionModel.getPrecisionMatrixAsVector();

        DenseMatrix64F precision = new DenseMatrix64F(dim, dim);
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                precision.set(i, j, temp[i * dim + j]);
            }
        }

        DenseMatrix64F Sigma = new DenseMatrix64F(dim, dim);
        CommonOps.invert(precision, Sigma);
        return Sigma;
    }


// Inside OUDiffusionModelDelegate

    /**
     * Stationary covariance Σ_stat (dim × dim) from CDI, in the ORIGINAL basis.
     */
    public DenseMatrix64F getStationaryCovarianceFromCDI(ContinuousDiffusionIntegrator cdi) {
        int dim = getDim();
        int precisionIndex = getEigenBufferOffsetIndex(0);

        double[] flat = new double[dim * dim];
        ((SafeMultivariateActualizedWithDriftIntegrator) cdi).getDiffusionStationaryVariance(precisionIndex, flat);

        DenseMatrix64F SigmaStat = new DenseMatrix64F(dim, dim);
        int k = 0;
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                SigmaStat.set(i, j, flat[k++]);
            }
        }
        return SigmaStat;
    }

    /**
     * Branch-specific actualization A(t) = Q(t) for the branch leading to 'node'.
     */
    public DenseMatrix64F getBranchTransitionMatrix(NodeRef node,
                                                    ContinuousDiffusionIntegrator cdi) {
        int dim = getDim();
        int bufferIndex = getMatrixBufferOffsetIndex(node.getNumber());

        double[] flat = new double[dim * dim];
        cdi.getBranchActualization(bufferIndex, flat); // already exists in CDI

        DenseMatrix64F A = new DenseMatrix64F(dim, dim);
        int k = 0;
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                A.set(i, j, flat[k++]);
            }
        }
        return A;
    }

    /**
     * Branch-specific variance V(t) for the branch leading to 'node'.
     * For non-integrated OU, dimTrait == dimProcess.
     */
    public DenseMatrix64F getBranchVarianceMatrix(NodeRef node,
                                                  ContinuousDiffusionIntegrator cdi) {
        int dim = getDim();
        int bufferIndex = getMatrixBufferOffsetIndex(node.getNumber());

        double[] flat = new double[dim * dim];
        ((SafeMultivariateActualizedWithDriftIntegrator) cdi).getBranchVariance(bufferIndex, flat);

        DenseMatrix64F V = new DenseMatrix64F(dim, dim);
        int k = 0;
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                V.set(i, j, flat[k++]);
            }
        }
        return V;
    }

    /** Global S in EJML form (you already basically had this in your cache). */
    public DenseMatrix64F getSelectionMatrixAsDense() {
        double[] flat = elasticModel.getStrengthOfSelectionMatrixAsVector();
        int dim = getDim();
        DenseMatrix64F S = new DenseMatrix64F(dim, dim);
        int k = 0;
        for (int i = 0; i < dim; ++i) {
            for (int j = 0; j < dim; ++j) {
                S.set(i, j, flat[k++]);
            }
        }
        return S;
    }

    /** Global μ(node) as a column vector. */
    public DenseMatrix64F getOptimalValueAsDense(NodeRef node) {
        double[] temp = getOptimalValues(node);
        int dim = getDim();
        DenseMatrix64F mu = new DenseMatrix64F(dim, 1);
        for (int i = 0; i < dim; ++i) {
            mu.set(i, 0, temp[i]);
        }
        return mu;
    }

    @Override
    public int getMatrixBufferOffsetIndex(int nodeIndex) {
        return super.getMatrixBufferOffsetIndex(nodeIndex);
    }

    @Override
    public int getEigenBufferOffsetIndex(int eigenIndex) {
        return super.getEigenBufferOffsetIndex(eigenIndex);
    }

    @Override
    public double[] getDriftRate(NodeRef node) {
        return super.getDriftRate(node);
    }

    // ----------------------------------------------------------------------
    // Model change plumbing
    // ----------------------------------------------------------------------

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        if (model == elasticModel) {
            // Eigen system / basis may have changed
            fireModelChanged(); //TODO add model inside here

        } else {
            super.handleModelChangedEvent(model, object, index);
        }

    }

    // ----------------------------------------------------------------------
    // Feature flags / accessors
    // ----------------------------------------------------------------------

    @Override
    public boolean hasDrift() {
        return true;
    }

    @Override
    public boolean hasActualization() {
        return true;
    }

    @Override
    public boolean hasDiagonalActualization() {
        // Still useful for BEAST internals, but no branching logic here.
        return elasticModel.isDiagonal();
    }

    public boolean hasBlockDiagActualization() {
        int[] blockStarts = elasticModel.getBlockStarts();
        int[] blockSizes  = elasticModel.getBlockSizes();
        return blockStarts != null && blockSizes != null && blockStarts.length > 0;
    }

    public MultivariateDiffusionModel getDiffusionModel() {
        return diffusionModel;
    }

    public boolean isSymmetric() {
        return elasticModel.isSymmetric();
    }

    public double[][] getStrengthOfSelection() {
        return elasticModel.getStrengthOfSelectionMatrix();
    }

    /**
     * In non-BLOCK_SCHUR, these are actual eigenvalues; in BLOCK_SCHUR, this
     * returns the diagonal entries of the real Schur form (1x1/2x2 blocks).
     */
    public double[] getEigenValuesStrengthOfSelection() {
        return elasticModel.getBasisEigenValues();
    }

    /**
     * Flattened basis matrix R (EJML row-major layout).
     */
    public double[] getEigenVectorsStrengthOfSelection() {
        return elasticModel.getBasisEigenVectors();
    }

    // ----------------------------------------------------------------------
    // CDI wiring
    // ----------------------------------------------------------------------

    @Override
    public void setDiffusionModels(ContinuousDiffusionIntegrator cdi, boolean flip) {

        super.setDiffusionModels(cdi, flip);

        final int precisionIndex = getEigenBufferOffsetIndex(0);
            cdi.setDiffusionStationaryVariance(
                    precisionIndex,
                    elasticModel.getBasisD(),
                    elasticModel.getBasisRotations()
            );
    }

    @Override
    public void updateDiffusionMatrices(ContinuousDiffusionIntegrator cdi,
                                        int[] branchIndices,
                                        double[] edgeLengths,
                                        int updateCount,
                                        boolean flip) {

        int[] probabilityIndices = new int[updateCount];

        for (int i = 0; i < updateCount; i++) {
            if (flip) {
                flipMatrixBufferOffset(branchIndices[i]);
            }
            probabilityIndices[i] = getMatrixBufferOffsetIndex(branchIndices[i]);
        }

        cdi.updateOrnsteinUhlenbeckDiffusionMatrices(
                getEigenBufferOffsetIndex(0),
                probabilityIndices,
                edgeLengths,
                getDriftRates(branchIndices, updateCount),
                elasticModel.getBasisD(),
                elasticModel.getBasisRotations(),
                updateCount);
    }

    // ----------------------------------------------------------------------
    // Gradients: variance wrt variance parameters
    // ----------------------------------------------------------------------

    @Override
    public DenseMatrix64F getGradientVarianceWrtVariance(NodeRef node,
                                                         ContinuousDiffusionIntegrator cdi,
                                                         ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                         DenseMatrix64F gradient) {
        // Root variance handled by base class; only non-root goes through strategies.
        if (tree.isRoot(node)) {
            return super.getGradientVarianceWrtVariance(node, cdi, likelihoodDelegate, gradient);
        }
        return actualizationStrategy.gradientVarianceWrtVariance(node, cdi, likelihoodDelegate, gradient);
    }

    // ----------------------------------------------------------------------
    // Gradients: variance wrt attenuation / selection
    // ----------------------------------------------------------------------

    DenseMatrix64F getGradientVarianceWrtAttenuation(NodeRef node,
                                                     ContinuousDiffusionIntegrator cdi,
                                                     BranchSufficientStatistics statistics,
                                                     DenseMatrix64F dSigma) {
        final int nodeIndex = node.getNumber();
        if (tree.isRoot(node)) {
            throw new IllegalArgumentException(
                    "Gradient wrt attenuation is not available for the root (node " + nodeIndex + ").");
        }
        return actualizationStrategy.gradientVarianceWrtAttenuation(node, cdi, statistics, dSigma);
    }

    // ----------------------------------------------------------------------
    // Gradients: displacement wrt drift / attenuation
    // ----------------------------------------------------------------------

    @Override
    DenseMatrix64F getGradientDisplacementWrtDrift(NodeRef node,
                                                   ContinuousDiffusionIntegrator cdi,
                                                   ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                   DenseMatrix64F gradient) {
        // Root vs non-root semantics are handled in getGradientDisplacementWrtRoot;
        // here we just implement the branch-local rule via strategy.
        return actualizationStrategy.gradientDisplacementWrtDrift(node, cdi, gradient);
    }

    DenseMatrix64F getGradientDisplacementWrtAttenuation(NodeRef node,
                                                         ContinuousDiffusionIntegrator cdi,
                                                         BranchSufficientStatistics statistics,
                                                         DenseMatrix64F gradient) {
        if (tree.isRoot(node)) {
            throw new IllegalArgumentException(
                    "Gradient wrt attenuation is not available for the root displacement.");
        }
        return actualizationStrategy.gradientDisplacementWrtAttenuation(node, cdi, statistics, gradient);
    }

    // ----------------------------------------------------------------------
    // Gradients: displacement wrt root
    //   (this is about root / fixed-root semantics only; actualization itself
    //    is handled by the strategy via rootGradient)
    // ----------------------------------------------------------------------

    @Override
    public double[] getGradientDisplacementWrtRoot(NodeRef node,
                                                   ContinuousDiffusionIntegrator cdi,
                                                   ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                   DenseMatrix64F gradient) {
        boolean fixedRoot =
                likelihoodDelegate.getRootProcessDelegate().getPseudoObservations()
                        == Double.POSITIVE_INFINITY;

        // Case 1: fixed root, child of root → push gradient down using branch actualization.
        if (fixedRoot && tree.isRoot(tree.getParent(node))) {
            return actualizationStrategy.rootGradient(getMatrixBufferOffsetIndex(node.getNumber()), cdi, gradient);
        }

        // Case 2: free root: gradient at the root itself is direct.
        if (!fixedRoot && tree.isRoot(node)) {
            return gradient.getData();
        }

        // Case 3: no contribution to root from other nodes.
        return new double[gradient.getNumRows()];
    }

    // ----------------------------------------------------------------------
    // Heritability-style summaries
    // ----------------------------------------------------------------------

    @Override
    public void getMeanTipVariances(final double priorSampleSize,
                                    final double[] treeLengths,
                                    final DenseMatrix64F traitVariance,
                                    final DenseMatrix64F varSum) {
        // Delegate does not care about diagonal vs general; strategy chooses basis.
        actualizationStrategy.meanTipVariances(priorSampleSize, treeLengths, traitVariance, varSum);
    }


    /**
     * NEW METHODS FOR OUDiffusionModelDelegate
     *
     * These implement direct backpropagation, matching the working Python code.
     * Add these to OUDiffusionModelDelegate class.
     */

    /**
     * Compute variance gradient using DIRECT BACKPROP.
     *
     * @param node Current node
     * @param cdi Continuous diffusion integrator
     * @param statistics Branch sufficient statistics
     * @param dL_dJ Gradient w.r.t. precision (from root)
     * @param dL_deta Gradient w.r.t. displacement (from root)
     * @param dL_dc Gradient w.r.t. constant term
     * @return Packed gradient w.r.t. S (variance contribution only)
     */
    public DenseMatrix64F getGradientVarianceWrtAttenuationDirectBackprop(
            NodeRef node,
            ContinuousDiffusionIntegrator cdi,
            BranchSufficientStatistics statistics,
            DenseMatrix64F dL_dJ,
            DenseMatrix64F dL_deta,
            double dL_dc) {

        // Delegate to the appropriate strategy based on whether S is diagonal
        return actualizationStrategy.gradientVarianceDirectBackprop(
                node, cdi, statistics, dL_dJ, dL_deta, dL_dc);
    }

    /**
     * Compute displacement gradient using DIRECT BACKPROP.
     *
     * @param node Current node
     * @param cdi Continuous diffusion integrator
     * @param statistics Branch sufficient statistics
     * @param dL_dJ Gradient w.r.t. precision (from root)
     * @param dL_deta Gradient w.r.t. displacement (from root)
     * @param dL_dc Gradient w.r.t. constant term
     * @return Packed gradient w.r.t. S (displacement contribution only)
     */
    public DenseMatrix64F getGradientDisplacementWrtAttenuationDirectBackprop(
            NodeRef node,
            ContinuousDiffusionIntegrator cdi,
            BranchSufficientStatistics statistics,
            DenseMatrix64F dL_dJ,
            DenseMatrix64F dL_deta,
            double dL_dc) {

        // Delegate to the appropriate strategy
        return actualizationStrategy.gradientDisplacementDirectBackprop(
                node, cdi, statistics, dL_dJ, dL_deta, dL_dc);
    }


    @Override
    public ContinuousDiffusionIntegrator createIntegrator(
            PrecisionType precisionType,
            int numTraits,
            int dimTrait,
            int dimProcess,
            int partialBufferCount,
            int matrixBufferCount,
            boolean allowSingular) {

        // SCALAR case – usually not what OU uses, but keep something safe.
        if (precisionType == PrecisionType.SCALAR) {
            return new ContinuousDiffusionIntegrator.Basic(
                    precisionType,
                    numTraits,
                    dimTrait,
                    dimProcess,
                    partialBufferCount,
                    matrixBufferCount
            );
        }

        // FULL precision, OU-specific branches (this is essentially the
        // old logic that lived in ContinuousDataLikelihoodDelegate, but
        // now localized here).

        if (this instanceof IntegratedOUDiffusionModelDelegate) {
            // Integrated OU: actualized + drift, half-dim trick
            // (this is exactly what you had before)
            return new SafeMultivariateActualizedWithDriftIntegrator(
                    precisionType,
                    numTraits,
                    dimTrait,
                    dimProcess,
                    partialBufferCount,
                    matrixBufferCount,
                    isSymmetric()
            );
        }

        if (hasDiagonalActualization()) {
            return new SafeMultivariateDiagonalActualizedWithDriftIntegrator(
                    precisionType,
                    numTraits,
                    dimTrait,
                    dimProcess,
                    partialBufferCount,
                    matrixBufferCount
            );
        }

        if (hasBlockDiagActualization()) {
            return new SafeMultivariateBlockDiagonalActualizedWithDriftIntegrator(
                    precisionType,
                    numTraits,
                    dimTrait,
                    dimProcess,
                    partialBufferCount,
                    matrixBufferCount
            );
        }

        // Generic OU case: SafeMultivariateActualizedWithDriftIntegrator
        return new SafeMultivariateActualizedWithDriftIntegrator(
                precisionType,
                numTraits,
                dimTrait,
                dimProcess,
                partialBufferCount,
                matrixBufferCount,
                isSymmetric()
        );
    }
}

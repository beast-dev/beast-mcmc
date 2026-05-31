package test.dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.continuous.ou.DenseSelectionMatrixParameterization;
import dr.evomodel.continuous.ou.MatrixExponentialUtils;
import dr.evomodel.continuous.ou.OUProcessModel;
import dr.evomodel.continuous.ou.SelectionMatrixParameterization;
import dr.evomodel.continuous.ou.SelectionMatrixParameterizationFactory;
import dr.evomodel.continuous.ou.canonical.CanonicalBranchWorkspace;
import dr.evomodel.continuous.ou.canonical.CanonicalPreparedBranchHandle;
import dr.evomodel.continuous.ou.blockdiagonal.BlockDiagonalBranchGradientWorkspace;
import dr.evomodel.continuous.ou.blockdiagonal.BlockDiagonalPreparedBranchBasis;
import dr.evomodel.continuous.ou.blockdiagonal.BlockDiagonalPreparedBranchHandle;
import dr.evomodel.continuous.ou.blockdiagonal.BlockDiagonalSelectionMatrixParameterization;
import dr.evomodel.continuous.ou.blockdiagonal.OrthogonalBlockDiagonalSelectionMatrixParameterization;
import dr.evomodel.treedatalikelihood.continuous.backprop.BlockDiagonalFrechetHelper;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalLocalTransitionAdjoints;
import dr.inference.model.BlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.MatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;

import java.util.Arrays;

public class BlockDiagonalSelectionMatrixParameterizationTest extends ContinuousTraitTest {

    private static final double TOL = 1.0e-9;

    public BlockDiagonalSelectionMatrixParameterizationTest(final String name) {
        super(name);
    }

    public void testGeneralBlockForwardMomentsMatchDenseVanLoan() {
        final BlockDiagonalPolarStableMatrixParameter blockParameter = buildBlockParameter();
        final SelectionMatrixParameterization selection =
                SelectionMatrixParameterizationFactory.create(blockParameter);

        assertTrue(selection instanceof BlockDiagonalSelectionMatrixParameterization);
        assertFalse(selection instanceof OrthogonalBlockDiagonalSelectionMatrixParameterization);

        final BlockDiagonalSelectionMatrixParameterization blockSelection =
                (BlockDiagonalSelectionMatrixParameterization) selection;
        final DenseSelectionMatrixParameterization denseSelection =
                new DenseSelectionMatrixParameterization(blockParameter);
        final MatrixParameter diffusionMatrix = makeMatrixParameter("Q.general.block", new double[][]{
                {1.20, 0.10, 0.05, 0.02},
                {0.10, 0.95, 0.04, 0.03},
                {0.05, 0.04, 1.15, 0.07},
                {0.02, 0.03, 0.07, 0.90}
        });
        final double dt = 0.37;
        final int dimension = 4;

        final double[] blockF = new double[dimension * dimension];
        final double[] denseF = new double[dimension * dimension];
        blockSelection.fillTransitionMatrixFlat(dt, blockF);
        denseSelection.fillTransitionMatrixFlat(dt, denseF);
        assertVectorEquals("transition matrix", denseF, blockF, TOL);

        final double[] blockV = new double[dimension * dimension];
        final double[] denseV = new double[dimension * dimension];
        blockSelection.fillTransitionCovarianceFlat(diffusionMatrix, dt, blockV);
        fillDenseTransitionCovarianceFlat(denseSelection, diffusionMatrix, dt, denseV);
        assertVectorEquals("transition covariance", denseV, blockV, TOL);

        final double[] stationaryMean = new double[]{0.3, -0.4, 0.2, 0.6};
        final CanonicalPreparedBranchHandle prepared = blockSelection.createPreparedBranchHandle();
        final CanonicalBranchWorkspace workspace = blockSelection.createBranchWorkspace();
        blockSelection.prepareBranch(dt, stationaryMean, prepared);

        final double[] preparedF = new double[dimension * dimension];
        final double[] preparedOffset = new double[dimension];
        final double[] preparedV = new double[dimension * dimension];
        assertTrue(blockSelection.fillTransitionMomentsPreparedFlat(
                prepared, diffusionMatrix, workspace, preparedF, preparedOffset, preparedV));
        assertVectorEquals("prepared transition matrix", blockF, preparedF, TOL);
        assertVectorEquals("prepared transition covariance", blockV, preparedV, TOL);

        final double[] expectedOffset = new double[dimension];
        fillTransitionOffset(blockF, stationaryMean, expectedOffset, dimension);
        assertVectorEquals("prepared transition offset", expectedOffset, preparedOffset, TOL);
    }

    public void testGeneralBlockNativeTransitionPullbackMatchesFiniteDifference() {
        final BlockDiagonalPolarStableMatrixParameter blockParameter = buildBlockParameter();
        final BlockDiagonalSelectionMatrixParameterization blockSelection =
                (BlockDiagonalSelectionMatrixParameterization)
                        SelectionMatrixParameterizationFactory.create(blockParameter);
        final int dimension = 4;
        final double dt = 0.29;
        final double[] stationaryMean = new double[]{0.3, -0.4, 0.2, 0.6};
        final double[] dLogL_dF = new double[]{
                0.20, -0.10, 0.05, 0.30,
                -0.25, 0.40, 0.12, -0.08,
                0.07, -0.18, 0.33, 0.09,
                -0.11, 0.24, -0.15, 0.28
        };
        final double[] dLogL_df = new double[]{0.11, -0.21, 0.17, 0.06};

        final NativeBlockGradient nativeGradient = accumulateNativeTransition(
                blockParameter, blockSelection, dt, stationaryMean, dLogL_dF, dLogL_df);
        final NativeBlockGradient finiteDifference = finiteDifferenceTransitionGradient(
                blockParameter, blockSelection, dt, stationaryMean, dLogL_dF, dLogL_df, dimension);
        assertNativeBlockGradientMatchesExpected(
                "transition", finiteDifference, nativeGradient, 2.0e-5);
    }

    public void testGeneralBlockNativeTransitionUsesSpecializedFrechetPath() {
        final BlockDiagonalPolarStableMatrixParameter blockParameter = buildBlockParameter();
        final BlockDiagonalSelectionMatrixParameterization blockSelection =
                (BlockDiagonalSelectionMatrixParameterization)
                        SelectionMatrixParameterizationFactory.create(blockParameter);
        final double dt = 0.29;
        final double[] stationaryMean = new double[]{0.3, -0.4, 0.2, 0.6};
        final double[] dLogL_dF = new double[]{
                0.20, -0.10, 0.05, 0.30,
                -0.25, 0.40, 0.12, -0.08,
                0.07, -0.18, 0.33, 0.09,
                -0.11, 0.24, -0.15, 0.28
        };
        final double[] dLogL_df = new double[]{0.11, -0.21, 0.17, 0.06};

        BlockDiagonalFrechetHelper.resetExactPlanInstrumentation();
        accumulateNativeTransition(blockParameter, blockSelection, dt, stationaryMean, dLogL_dF, dLogL_df);

        assertTrue("general block transition should use block Frechet adjoint",
                blockFrechetApplicationCount() > 0L);
    }

    public void testGeneralBlockPublicTransitionGradientUsesNativeParameterLayout() {
        final BlockDiagonalPolarStableMatrixParameter blockParameter = buildBlockParameter();
        final BlockDiagonalSelectionMatrixParameterization blockSelection =
                (BlockDiagonalSelectionMatrixParameterization)
                        SelectionMatrixParameterizationFactory.create(blockParameter);
        final int dimension = 4;
        final double dt = 0.29;
        final double[] stationaryMean = new double[]{0.3, -0.4, 0.2, 0.6};
        final double[] dLogL_dF = new double[]{
                0.20, -0.10, 0.05, 0.30,
                -0.25, 0.40, 0.12, -0.08,
                0.07, -0.18, 0.33, 0.09,
                -0.11, 0.24, -0.15, 0.28
        };
        final double[] dLogL_df = new double[]{0.11, -0.21, 0.17, 0.06};

        final NativeBlockGradient finiteDifference = finiteDifferenceTransitionGradient(
                blockParameter, blockSelection, dt, stationaryMean, dLogL_dF, dLogL_df, dimension);
        final double[] expected = assembleFullBlockGradient(blockParameter, finiteDifference, dimension);

        final double[] actual = new double[blockSelection.getSelectionGradientDimension()];
        BlockDiagonalFrechetHelper.resetExactPlanInstrumentation();
        blockSelection.accumulateGradientFromTransitionFlat(
                dt, stationaryMean, dLogL_dF, dLogL_df, actual);

        assertVectorEquals("public transition native parameter gradient", expected, actual, 2.0e-5);
        assertTrue("public transition should use block Frechet adjoint",
                blockFrechetApplicationCount() > 0L);
    }

    public void testGeneralBlockNativeCovariancePullbackMatchesFiniteDifference() {
        final BlockDiagonalPolarStableMatrixParameter blockParameter = buildBlockParameter();
        final BlockDiagonalSelectionMatrixParameterization blockSelection =
                (BlockDiagonalSelectionMatrixParameterization)
                        SelectionMatrixParameterizationFactory.create(blockParameter);
        final MatrixParameter diffusionMatrix = makeDiffusionMatrix();
        final int dimension = 4;
        final double dt = 0.37;
        final double[] dLogL_dV = new double[]{
                0.30, -0.12, 0.07, 0.04,
                -0.05, 0.25, -0.18, 0.10,
                0.02, -0.16, 0.40, -0.09,
                0.06, 0.13, -0.03, 0.22
        };

        final NativeBlockGradient nativeGradient = accumulateNativeCovariance(
                blockParameter, blockSelection, diffusionMatrix, dt, dLogL_dV);
        final NativeBlockGradient finiteDifference = finiteDifferenceCovarianceGradient(
                blockParameter, blockSelection, diffusionMatrix, dt, dLogL_dV, dimension);
        assertNativeBlockGradientMatchesExpected(
                "covariance", finiteDifference, nativeGradient, 2.0e-5);
    }

    public void testGeneralBlockNativeCovarianceUsesSpecializedFrechetPath() {
        final BlockDiagonalPolarStableMatrixParameter blockParameter = buildBlockParameter();
        final BlockDiagonalSelectionMatrixParameterization blockSelection =
                (BlockDiagonalSelectionMatrixParameterization)
                        SelectionMatrixParameterizationFactory.create(blockParameter);
        final MatrixParameter diffusionMatrix = makeDiffusionMatrix();
        final double dt = 0.37;
        final double[] dLogL_dV = new double[]{
                0.30, -0.12, 0.07, 0.04,
                -0.05, 0.25, -0.18, 0.10,
                0.02, -0.16, 0.40, -0.09,
                0.06, 0.13, -0.03, 0.22
        };

        BlockDiagonalFrechetHelper.resetExactPlanInstrumentation();
        accumulateNativeCovariance(blockParameter, blockSelection, diffusionMatrix, dt, dLogL_dV);

        assertTrue("general block covariance should use block Frechet adjoint",
                blockFrechetApplicationCount() > 0L);
    }

    public void testGeneralBlockNativeDiffusionPullbackMatchesDenseGradient() {
        final BlockDiagonalPolarStableMatrixParameter blockParameter = buildBlockParameter();
        final BlockDiagonalSelectionMatrixParameterization blockSelection =
                (BlockDiagonalSelectionMatrixParameterization)
                        SelectionMatrixParameterizationFactory.create(blockParameter);
        final MatrixParameter diffusionMatrix = makeDiffusionMatrix();
        final int dimension = 4;
        final double dt = 0.41;
        final double[] dLogL_dV = new double[]{
                0.21, 0.05, -0.07, 0.02,
                0.08, 0.31, 0.09, -0.11,
                -0.04, 0.12, 0.27, 0.14,
                0.03, -0.06, 0.10, 0.24
        };

        final OUProcessModel denseProcess = makeOuProcess(blockParameter, diffusionMatrix, dimension);
        final double[] expected = new double[dimension * dimension];
        denseProcess.accumulateDiffusionGradientFlat(dt, dLogL_dV, false, expected);

        final double[] actual = new double[dimension * dimension];
        blockSelection.accumulateDiffusionGradientFlat(
                diffusionMatrix, dt, dLogL_dV, false, actual);
        assertVectorEquals("diffusion gradient", expected, actual, 5.0e-8);
    }

    public void testGeneralBlockPreparedNativeAdjointsMatchFiniteDifference() {
        final BlockDiagonalPolarStableMatrixParameter blockParameter = buildBlockParameter();
        final BlockDiagonalSelectionMatrixParameterization blockSelection =
                (BlockDiagonalSelectionMatrixParameterization)
                        SelectionMatrixParameterizationFactory.create(blockParameter);
        final MatrixParameter diffusionMatrix = makeDiffusionMatrix();
        final int dimension = 4;
        final double dt = 0.31;
        final double[] stationaryMean = new double[]{0.3, -0.4, 0.2, 0.6};
        final double[] dLogL_dF = new double[]{
                -0.17, 0.06, 0.11, -0.03,
                0.14, -0.22, 0.05, 0.19,
                -0.08, 0.10, -0.13, 0.07,
                0.03, -0.16, 0.09, -0.24
        };
        final double[] dLogL_df = new double[]{-0.09, 0.15, 0.07, -0.12};
        final double[] dLogL_dV = new double[]{
                0.27, 0.04, -0.12, 0.08,
                -0.02, 0.36, 0.11, -0.06,
                -0.07, 0.13, 0.18, 0.05,
                0.03, -0.09, 0.02, 0.29
        };

        final NativeBlockGradient finiteDifference = finiteDifferenceTransitionGradient(
                blockParameter, blockSelection, dt, stationaryMean, dLogL_dF, dLogL_df, dimension);
        addNativeBlockGradient(
                finiteDifference,
                finiteDifferenceCovarianceGradient(
                        blockParameter, blockSelection, diffusionMatrix, dt, dLogL_dV, dimension));

        final CanonicalPreparedBranchHandle prepared = blockSelection.createPreparedBranchHandle();
        final CanonicalBranchWorkspace workspace = blockSelection.createBranchWorkspace();
        blockSelection.prepareBranch(dt, stationaryMean, prepared);
        blockSelection.prepareBranchCovariance(
                (BlockDiagonalPreparedBranchBasis) ((BlockDiagonalPreparedBranchHandle) prepared).getBasis(),
                diffusionMatrix,
                (BlockDiagonalBranchGradientWorkspace) workspace);

        final CanonicalLocalTransitionAdjoints adjoints = new CanonicalLocalTransitionAdjoints(dimension);
        System.arraycopy(dLogL_dF, 0, adjoints.dLogL_dF, 0, dLogL_dF.length);
        System.arraycopy(dLogL_df, 0, adjoints.dLogL_df, 0, dLogL_df.length);
        System.arraycopy(dLogL_dV, 0, adjoints.dLogL_dOmega, 0, dLogL_dV.length);

        final NativeBlockGradient nativeGradient = new NativeBlockGradient(
                blockParameter.getBlockDiagonalNParameters(),
                blockParameter.getCompressedDDimension(),
                dimension);
        blockSelection.accumulateNativeGradientFromAdjointsPreparedFlat(
                prepared,
                diffusionMatrix,
                adjoints,
                workspace,
                nativeGradient.compressed,
                nativeGradient.rotation);
        blockParameter.chainGradient(nativeGradient.compressed, nativeGradient.nativeBlock);

        assertNativeBlockGradientMatchesExpected(
                "prepared selection", finiteDifference, nativeGradient, 2.0e-5);
    }

    public void testGeneralBlockPreparedCovarianceReuseAcrossWorkspacesMatchesFiniteDifference() {
        final BlockDiagonalPolarStableMatrixParameter blockParameter = buildBlockParameter();
        final BlockDiagonalSelectionMatrixParameterization blockSelection =
                (BlockDiagonalSelectionMatrixParameterization)
                        SelectionMatrixParameterizationFactory.create(blockParameter);
        final MatrixParameter diffusionMatrix = makeDiffusionMatrix();
        final int dimension = 4;
        final double dt = 0.33;
        final double[] stationaryMean = new double[]{0.3, -0.4, 0.2, 0.6};
        final double[] dLogL_dV = new double[]{
                0.23, -0.10, 0.06, 0.11,
                -0.04, 0.32, -0.15, 0.07,
                0.05, -0.13, 0.37, -0.08,
                0.02, 0.09, -0.05, 0.26
        };

        final NativeBlockGradient finiteDifference = finiteDifferenceCovarianceGradient(
                blockParameter, blockSelection, diffusionMatrix, dt, dLogL_dV, dimension);

        final CanonicalPreparedBranchHandle prepared = blockSelection.createPreparedBranchHandle();
        final CanonicalBranchWorkspace prepareWorkspace = blockSelection.createBranchWorkspace();
        final CanonicalBranchWorkspace pullbackWorkspace = blockSelection.createBranchWorkspace();
        blockSelection.prepareBranch(dt, stationaryMean, prepared);
        blockSelection.prepareBranchCovariance(
                (BlockDiagonalPreparedBranchBasis) ((BlockDiagonalPreparedBranchHandle) prepared).getBasis(),
                diffusionMatrix,
                (BlockDiagonalBranchGradientWorkspace) prepareWorkspace);

        final CanonicalLocalTransitionAdjoints adjoints = new CanonicalLocalTransitionAdjoints(dimension);
        System.arraycopy(dLogL_dV, 0, adjoints.dLogL_dOmega, 0, dLogL_dV.length);

        final NativeBlockGradient nativeGradient = new NativeBlockGradient(
                blockParameter.getBlockDiagonalNParameters(),
                blockParameter.getCompressedDDimension(),
                dimension);
        blockSelection.accumulateNativeGradientFromAdjointsPreparedFlat(
                prepared,
                diffusionMatrix,
                adjoints,
                pullbackWorkspace,
                nativeGradient.compressed,
                nativeGradient.rotation);
        blockParameter.chainGradient(nativeGradient.compressed, nativeGradient.nativeBlock);

        assertNativeBlockGradientMatchesExpected(
                "prepared covariance reused in fresh workspace",
                finiteDifference,
                nativeGradient,
                2.0e-5);
    }

    public void testGeneralBlockPreparedDelayedDiffusionPullbackMatchesDenseGradient() {
        final BlockDiagonalPolarStableMatrixParameter blockParameter = buildBlockParameter();
        final BlockDiagonalSelectionMatrixParameterization blockSelection =
                (BlockDiagonalSelectionMatrixParameterization)
                        SelectionMatrixParameterizationFactory.create(blockParameter);
        final MatrixParameter diffusionMatrix = makeDiffusionMatrix();
        final int dimension = 4;
        final double dt = 0.27;
        final double[] stationaryMean = new double[]{0.3, -0.4, 0.2, 0.6};
        final double[] dLogL_dV = new double[]{
                0.19, -0.07, 0.04, 0.08,
                -0.03, 0.26, -0.11, 0.05,
                0.06, -0.14, 0.31, -0.02,
                0.09, 0.03, -0.08, 0.23
        };

        final OUProcessModel denseProcess = makeOuProcess(blockParameter, diffusionMatrix, dimension);
        final NativeBlockGradient finiteDifference = finiteDifferenceCovarianceGradient(
                blockParameter, blockSelection, diffusionMatrix, dt, dLogL_dV, dimension);
        final double[] expectedDiffusionGradient = new double[dimension * dimension];
        denseProcess.accumulateDiffusionGradientFlat(
                dt, dLogL_dV, false, expectedDiffusionGradient);

        final CanonicalPreparedBranchHandle prepared = blockSelection.createPreparedBranchHandle();
        final CanonicalBranchWorkspace workspace = blockSelection.createBranchWorkspace();
        blockSelection.prepareBranch(dt, stationaryMean, prepared);

        final CanonicalLocalTransitionAdjoints adjoints = new CanonicalLocalTransitionAdjoints(dimension);
        System.arraycopy(dLogL_dV, 0, adjoints.dLogL_dOmega, 0, dLogL_dV.length);

        final NativeBlockGradient nativeGradient = new NativeBlockGradient(
                blockParameter.getBlockDiagonalNParameters(),
                blockParameter.getCompressedDDimension(),
                dimension);
        final double[] dBasisDiffusionGradient = new double[dimension * dimension];
        assertTrue("general block path should support delayed diffusion rotation",
                blockSelection.supportsDelayedDiffusionGradientRotation());
        blockSelection.accumulateNativeSelectionAndDiffusionGradientFromAdjointsPreparedFlat(
                prepared,
                diffusionMatrix,
                adjoints,
                workspace,
                nativeGradient.compressed,
                nativeGradient.rotation,
                true,
                dBasisDiffusionGradient);
        blockParameter.chainGradient(nativeGradient.compressed, nativeGradient.nativeBlock);

        assertNativeBlockGradientMatchesExpected(
                "prepared delayed covariance selection",
                finiteDifference,
                nativeGradient,
                2.0e-5);

        final double[] actualDiffusionGradient = new double[dimension * dimension];
        blockSelection.finishDiffusionGradientFromDBasisFlat(
                dBasisDiffusionGradient,
                actualDiffusionGradient,
                workspace);
        assertVectorEquals(
                "prepared delayed diffusion gradient",
                expectedDiffusionGradient,
                actualDiffusionGradient,
                5.0e-8);
    }

    private static BlockDiagonalPolarStableMatrixParameter buildBlockParameter() {
        final MatrixParameter r = makeMatrixParameter("R.general.block", new double[][]{
                {1.00, 0.20, 0.00, 0.10},
                {0.10, 1.10, 0.15, 0.00},
                {0.00, 0.05, 0.90, 0.25},
                {0.05, 0.00, 0.10, 1.05}
        });
        final Parameter scalar = new Parameter.Default(0);
        final Parameter rho = new Parameter.Default("rho.general.block", new double[]{0.7, 0.9});
        final Parameter theta = new Parameter.Default("theta.general.block", new double[]{0.25, -0.35});
        final Parameter t = new Parameter.Default("t.general.block", new double[]{0.10, -0.08});
        return new BlockDiagonalPolarStableMatrixParameter("block.general", r, scalar, rho, theta, t);
    }

    private static MatrixParameter makeDiffusionMatrix() {
        return makeMatrixParameter("Q.general.block", new double[][]{
                {1.20, 0.10, 0.05, 0.02},
                {0.10, 0.95, 0.04, 0.03},
                {0.05, 0.04, 1.15, 0.07},
                {0.02, 0.03, 0.07, 0.90}
        });
    }

    private static OUProcessModel makeOuProcess(final MatrixParameterInterface driftMatrix,
                                                final MatrixParameter diffusionMatrix,
                                                final int dimension) {
        final Parameter mean = new Parameter.Default("mu.general.block", new double[]{0.3, -0.4, 0.2, 0.6});
        final MatrixParameter initialCovariance = makeMatrixParameter("P0.general.block", new double[][]{
                {1.1, 0.1, 0.0, 0.0},
                {0.1, 1.2, 0.1, 0.0},
                {0.0, 0.1, 1.3, 0.1},
                {0.0, 0.0, 0.1, 1.4}
        });
        return new OUProcessModel(
                "ou.general.block",
                dimension,
                driftMatrix,
                diffusionMatrix,
                mean,
                initialCovariance,
                OUProcessModel.CovarianceGradientMethod.STATIONARY_LYAPUNOV);
    }

    private static NativeBlockGradient accumulateNativeTransition(
            final BlockDiagonalPolarStableMatrixParameter blockParameter,
            final BlockDiagonalSelectionMatrixParameterization blockSelection,
            final double dt,
            final double[] stationaryMean,
            final double[] dLogL_dF,
            final double[] dLogL_df) {
        final NativeBlockGradient nativeGradient = new NativeBlockGradient(
                blockParameter.getBlockDiagonalNParameters(),
                blockParameter.getCompressedDDimension(),
                blockParameter.getRowDimension());
        blockSelection.accumulateNativeGradientFromTransitionFlat(
                dt,
                stationaryMean,
                dLogL_dF,
                dLogL_df,
                nativeGradient.compressed,
                nativeGradient.rotation);
        blockParameter.chainGradient(nativeGradient.compressed, nativeGradient.nativeBlock);
        return nativeGradient;
    }

    private static NativeBlockGradient accumulateNativeCovariance(
            final BlockDiagonalPolarStableMatrixParameter blockParameter,
            final BlockDiagonalSelectionMatrixParameterization blockSelection,
            final MatrixParameterInterface diffusionMatrix,
            final double dt,
            final double[] dLogL_dV) {
        final NativeBlockGradient nativeGradient = new NativeBlockGradient(
                blockParameter.getBlockDiagonalNParameters(),
                blockParameter.getCompressedDDimension(),
                blockParameter.getRowDimension());
        blockSelection.accumulateNativeGradientFromCovarianceStationaryFlat(
                diffusionMatrix,
                dt,
                dLogL_dV,
                nativeGradient.compressed,
                nativeGradient.rotation);
        blockParameter.chainGradient(nativeGradient.compressed, nativeGradient.nativeBlock);
        return nativeGradient;
    }

    private static long blockFrechetApplicationCount() {
        return BlockDiagonalFrechetHelper.getExactPlanApplicationCount()
                + BlockDiagonalFrechetHelper.getEqualDiagonalPlanApplicationCount();
    }

    private static NativeBlockGradient finiteDifferenceTransitionGradient(
            final BlockDiagonalPolarStableMatrixParameter blockParameter,
            final BlockDiagonalSelectionMatrixParameterization blockSelection,
            final double dt,
            final double[] stationaryMean,
            final double[] dLogL_dF,
            final double[] dLogL_df,
            final int dimension) {
        final NativeBlockGradient gradient = new NativeBlockGradient(
                blockParameter.getBlockDiagonalNParameters(),
                blockParameter.getCompressedDDimension(),
                dimension);
        final Objective objective = new Objective() {
            @Override
            public double value() {
                return transitionObjective(
                        blockSelection, dt, stationaryMean, dLogL_dF, dLogL_df, dimension);
            }
        };
        fillFiniteDifferenceNativeBlockGradient(blockParameter, objective, gradient.nativeBlock);
        fillFiniteDifferenceRotationGradient(blockParameter, objective, gradient.rotation, dimension);
        return gradient;
    }

    private static NativeBlockGradient finiteDifferenceCovarianceGradient(
            final BlockDiagonalPolarStableMatrixParameter blockParameter,
            final BlockDiagonalSelectionMatrixParameterization blockSelection,
            final MatrixParameterInterface diffusionMatrix,
            final double dt,
            final double[] dLogL_dV,
            final int dimension) {
        final NativeBlockGradient gradient = new NativeBlockGradient(
                blockParameter.getBlockDiagonalNParameters(),
                blockParameter.getCompressedDDimension(),
                dimension);
        final Objective objective = new Objective() {
            @Override
            public double value() {
                return covarianceObjective(blockSelection, diffusionMatrix, dt, dLogL_dV, dimension);
            }
        };
        fillFiniteDifferenceNativeBlockGradient(blockParameter, objective, gradient.nativeBlock);
        fillFiniteDifferenceRotationGradient(blockParameter, objective, gradient.rotation, dimension);
        return gradient;
    }

    private static void fillFiniteDifferenceNativeBlockGradient(
            final BlockDiagonalPolarStableMatrixParameter blockParameter,
            final Objective objective,
            final double[] out) {
        int index = 0;
        if (blockParameter.hasLeadingOneByOneBlock()) {
            out[index++] = finiteDifference(blockParameter.getScalarBlockParameter(), 0, objective);
        }
        for (int family = 0; family < blockParameter.getTwoByTwoParameterFamilyCount(); ++family) {
            final Parameter parameter = blockParameter.getTwoByTwoBlockParameter(family);
            for (int block = 0; block < blockParameter.getNum2x2Blocks(); ++block) {
                out[index++] = finiteDifference(parameter, block, objective);
            }
        }
    }

    private static void fillFiniteDifferenceRotationGradient(
            final BlockDiagonalPolarStableMatrixParameter blockParameter,
            final Objective objective,
            final double[] out,
            final int dimension) {
        final MatrixParameter rotation = blockParameter.getRotationMatrixParameter();
        for (int row = 0; row < dimension; ++row) {
            for (int col = 0; col < dimension; ++col) {
                out[row * dimension + col] = finiteDifference(rotation, row, col, objective);
            }
        }
    }

    private static double finiteDifference(final Parameter parameter,
                                           final int index,
                                           final Objective objective) {
        final double saved = parameter.getParameterValue(index);
        final double step = 1.0e-6;
        parameter.setParameterValue(index, saved + step);
        final double plus = objective.value();
        parameter.setParameterValue(index, saved - step);
        final double minus = objective.value();
        parameter.setParameterValue(index, saved);
        return (plus - minus) / (2.0 * step);
    }

    private static double finiteDifference(final MatrixParameter parameter,
                                           final int row,
                                           final int col,
                                           final Objective objective) {
        final double saved = parameter.getParameterValue(row, col);
        final double step = 1.0e-6;
        parameter.setParameterValue(row, col, saved + step);
        final double plus = objective.value();
        parameter.setParameterValue(row, col, saved - step);
        final double minus = objective.value();
        parameter.setParameterValue(row, col, saved);
        return (plus - minus) / (2.0 * step);
    }

    private static double transitionObjective(
            final BlockDiagonalSelectionMatrixParameterization blockSelection,
            final double dt,
            final double[] stationaryMean,
            final double[] dLogL_dF,
            final double[] dLogL_df,
            final int dimension) {
        final double[] transitionMatrix = new double[dimension * dimension];
        blockSelection.fillTransitionMatrixFlat(dt, transitionMatrix);
        double objective = 0.0;
        for (int row = 0; row < dimension; ++row) {
            final int rowOffset = row * dimension;
            double transformedMean = 0.0;
            for (int col = 0; col < dimension; ++col) {
                final double f = transitionMatrix[rowOffset + col];
                objective += dLogL_dF[rowOffset + col] * f;
                transformedMean += f * stationaryMean[col];
            }
            objective += dLogL_df[row] * (stationaryMean[row] - transformedMean);
        }
        return objective;
    }

    private static double covarianceObjective(
            final BlockDiagonalSelectionMatrixParameterization blockSelection,
            final MatrixParameterInterface diffusionMatrix,
            final double dt,
            final double[] dLogL_dV,
            final int dimension) {
        final double[] covariance = new double[dimension * dimension];
        blockSelection.fillTransitionCovarianceFlat(diffusionMatrix, dt, covariance);
        double objective = 0.0;
        for (int i = 0; i < covariance.length; ++i) {
            objective += dLogL_dV[i] * covariance[i];
        }
        return objective;
    }

    private static void addNativeBlockGradient(final NativeBlockGradient target,
                                               final NativeBlockGradient increment) {
        for (int i = 0; i < target.nativeBlock.length; ++i) {
            target.nativeBlock[i] += increment.nativeBlock[i];
        }
        for (int i = 0; i < target.rotation.length; ++i) {
            target.rotation[i] += increment.rotation[i];
        }
    }

    private static double[] assembleFullBlockGradient(
            final BlockDiagonalPolarStableMatrixParameter blockParameter,
            final NativeBlockGradient gradient,
            final int dimension) {
        final int nativeBlockDim = blockParameter.getBlockDiagonalNParameters();
        final double[] out = new double[nativeBlockDim + dimension * dimension];
        System.arraycopy(gradient.nativeBlock, 0, out, 0, nativeBlockDim);
        int index = nativeBlockDim;
        for (int col = 0; col < dimension; ++col) {
            for (int row = 0; row < dimension; ++row) {
                out[index++] = gradient.rotation[row * dimension + col];
            }
        }
        return out;
    }

    private static void assertNativeBlockGradientMatchesExpected(
            final String label,
            final NativeBlockGradient expected,
            final NativeBlockGradient nativeGradient,
            final double tolerance) {
        assertVectorEquals(label + " native block",
                expected.nativeBlock,
                nativeGradient.nativeBlock,
                tolerance);
        assertVectorEquals(label + " rotation",
                expected.rotation,
                nativeGradient.rotation,
                tolerance);
    }

    private interface Objective {
        double value();
    }

    private static final class NativeBlockGradient {
        final double[] nativeBlock;
        final double[] compressed;
        final double[] rotation;

        private NativeBlockGradient(final int nativeBlockDimension,
                                    final int compressedDimension,
                                    final int matrixDimension) {
            this.nativeBlock = new double[nativeBlockDimension];
            this.compressed = new double[compressedDimension];
            this.rotation = new double[matrixDimension * matrixDimension];
        }
    }

    private static void fillDenseTransitionCovarianceFlat(final SelectionMatrixParameterization selection,
                                                          final MatrixParameterInterface diffusionMatrix,
                                                          final double dt,
                                                          final double[] out) {
        final int d = selection.getDimension();
        final double[] a = new double[d * d];
        final double[] q = new double[d * d];
        final double[] vanLoan = new double[4 * d * d];
        final double[] vanLoanExp = new double[4 * d * d];

        selection.fillSelectionMatrixFlat(a);
        fillDenseMatrix(diffusionMatrix, q);

        for (int i = 0; i < d; ++i) {
            for (int j = 0; j < d; ++j) {
                vanLoan[i * 2 * d + j] = -dt * a[i * d + j];
                vanLoan[i * 2 * d + j + d] = dt * q[i * d + j];
                vanLoan[(i + d) * 2 * d + j] = 0.0;
                vanLoan[(i + d) * 2 * d + j + d] = dt * a[j * d + i];
            }
        }

        MatrixExponentialUtils.expmFlat(vanLoan, vanLoanExp, 2 * d);

        for (int i = 0; i < d; ++i) {
            final int iOff = i * d;
            for (int j = 0; j < d; ++j) {
                double sum = 0.0;
                for (int k = 0; k < d; ++k) {
                    sum += vanLoanExp[i * 2 * d + k + d] * vanLoanExp[j * 2 * d + k];
                }
                out[iOff + j] = sum;
            }
        }
        MatrixExponentialUtils.symmetrizeFlat(out, d);
    }

    private static void fillTransitionOffset(final double[] transitionMatrix,
                                             final double[] stationaryMean,
                                             final double[] out,
                                             final int dimension) {
        for (int i = 0; i < dimension; ++i) {
            double transformedMean = 0.0;
            final int rowOffset = i * dimension;
            for (int j = 0; j < dimension; ++j) {
                transformedMean += transitionMatrix[rowOffset + j] * stationaryMean[j];
            }
            out[i] = stationaryMean[i] - transformedMean;
        }
    }

    private static void fillDenseMatrix(final MatrixParameterInterface matrix, final double[] out) {
        final int dimension = matrix.getRowDimension();
        for (int i = 0; i < dimension; ++i) {
            for (int j = 0; j < dimension; ++j) {
                out[i * dimension + j] = matrix.getParameterValue(i, j);
            }
        }
    }

    private static MatrixParameter makeMatrixParameter(final String name,
                                                       final double[][] values) {
        final int rows = values.length;
        final int cols = values[0].length;
        final Parameter[] columns = new Parameter[cols];
        for (int j = 0; j < cols; ++j) {
            final double[] column = new double[rows];
            for (int i = 0; i < rows; ++i) {
                column[i] = values[i][j];
            }
            columns[j] = new Parameter.Default(name + ".col" + j, column);
        }
        return new MatrixParameter(name, columns);
    }

    private static void assertVectorEquals(final String label,
                                           final double[] expected,
                                           final double[] actual,
                                           final double tolerance) {
        assertEquals(label + " length", expected.length, actual.length);
        for (int i = 0; i < expected.length; ++i) {
            if (Math.abs(expected[i] - actual[i]) > tolerance) {
                fail(label + "[" + i + "] expected=" + expected[i]
                        + " actual=" + actual[i]
                        + " expectedVector=" + Arrays.toString(expected)
                        + " actualVector=" + Arrays.toString(actual));
            }
        }
    }
}

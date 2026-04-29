package test.dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.TreeTrait;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.branchratemodel.ArbitraryBranchRates;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.continuous.MultivariateElasticModel;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.BranchSpecificGradient;
import dr.evomodel.treedatalikelihood.continuous.BranchSpecificOptimaGradient;
import dr.evomodel.treedatalikelihood.continuous.ContinuousDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitDataModel;
import dr.evomodel.treedatalikelihood.continuous.ContinuousTraitGradientForBranch;
import dr.evomodel.treedatalikelihood.continuous.DiffusionProcessDelegate;
import dr.evomodel.treedatalikelihood.continuous.OUDiffusionModelDelegate;
import dr.evomodel.treedatalikelihood.continuous.TimeSeriesOUGaussianBranchTransitionProvider;
import dr.evomodel.treedatalikelihood.hmc.AbstractDiffusionGradient;
import dr.evomodel.treedatalikelihood.preorder.BranchConditionalDistributionDelegate;
import dr.evomodel.treedatalikelihood.preorder.BranchSufficientStatistics;
import dr.inference.model.GivensRotationMatrixParameter;
import dr.evomodel.treedatalikelihood.continuous.cdi.PrecisionType;
import dr.inference.hmc.JointGradient;
import dr.inference.model.CompoundParameter;
import dr.inference.model.DiagonalMatrix;
import dr.inference.model.MatrixParameter;
import dr.inference.model.OrthogonalBlockDiagonalPolarStableMatrixParameter;
import dr.inference.model.Parameter;
import dr.inference.timeseries.engine.gaussian.CanonicalLocalTransitionAdjoints;
import dr.inference.timeseries.gaussian.OrthogonalBlockDiagonalSelectionMatrixParameterization;
import dr.math.MultivariateFunction;
import dr.math.NumericalDerivative;
import junit.framework.Test;
import junit.framework.TestSuite;

import java.lang.reflect.Method;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import static dr.evomodel.treedatalikelihood.hmc.AbstractDiffusionGradient.ParameterDiffusionGradient.createDiagonalAttenuationGradient;

/**
 * Validates the opt-in forward bridge between {@link OUDiffusionModelDelegate}
 * and the shared time-series OU branch kernel on a small test tree.
 */
public class OUDiffusionKernelBridgeValidationTest extends ContinuousTraitTest {

    private static final String USE_BRIDGE_PROPERTY =
            "dr.evomodel.treedatalikelihood.continuous.useTimeSeriesOUBridge";
    private static final String VALIDATION_PROPERTY =
            "dr.evomodel.treedatalikelihood.continuous.validateTimeSeriesOUBridge";

    public OUDiffusionKernelBridgeValidationTest(final String name) {
        super(name);
    }

    public void testDelegateForwardPathMatchesKernelBridge_DiagonalOU() {
        final String previousValue = System.getProperty(VALIDATION_PROPERTY);
        System.setProperty(VALIDATION_PROPERTY, "true");
        try {
            final double logLikelihood = buildDiagonalOuTreeLikelihood().getLogLikelihood();
            assertTrue("Likelihood should be finite when kernel bridge validation is enabled",
                    Double.isFinite(logLikelihood));
        } finally {
            if (previousValue == null) {
                System.clearProperty(VALIDATION_PROPERTY);
            } else {
                System.setProperty(VALIDATION_PROPERTY, previousValue);
            }
        }
    }

    public void testDelegateForwardPathCanBeReplacedByKernelBridge_DiagonalOU() {
        final String previousUse = System.getProperty(USE_BRIDGE_PROPERTY);
        final String previousValidation = System.getProperty(VALIDATION_PROPERTY);
        System.setProperty(USE_BRIDGE_PROPERTY, "true");
        System.setProperty(VALIDATION_PROPERTY, "true");
        try {
            final double logLikelihood = buildDiagonalOuTreeLikelihood().getLogLikelihood();
            assertTrue("Likelihood should be finite when kernel bridge replacement is enabled",
                    Double.isFinite(logLikelihood));
        } finally {
            if (previousUse == null) {
                System.clearProperty(USE_BRIDGE_PROPERTY);
            } else {
                System.setProperty(USE_BRIDGE_PROPERTY, previousUse);
            }
            if (previousValidation == null) {
                System.clearProperty(VALIDATION_PROPERTY);
            } else {
                System.setProperty(VALIDATION_PROPERTY, previousValidation);
            }
        }
    }

    public void testCanonicalBranchGradientIsFinite_ForBranchSpecificOptima() {
        final BranchSpecificSetup setup = buildBranchSpecificOuTreeLikelihood();
        final ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient traitGradient =
                new ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient(
                        dimTrait,
                        treeModel,
                        setup.likelihoodDelegate,
                        Arrays.asList(ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter.WRT_BRANCH_SPECIFIC_DRIFT));
        final BranchSpecificOptimaGradient gradient =
                new BranchSpecificOptimaGradient("trait", setup.treeDataLikelihood, setup.likelihoodDelegate,
                        traitGradient, setup.branchOptima);

        final double[] analytic = gradient.getGradientLogDensity();
        final int numBranches = treeModel.getNodeCount() - 1;
        for (int branch = 0; branch < numBranches; ++branch) {
            if (treeModel.isExternal(treeModel.getNode(branch))) {
                continue;
            }
            for (int trait = 0; trait < dimTrait; ++trait) {
                final int index = branch + numBranches * trait;
                assertTrue("branch-specific optima gradient idx=" + index + " must be finite",
                        Double.isFinite(analytic[index]));
            }
        }
    }

    public void testCanonicalBranchGradientMatchesFD_ForDiagonalSelectionStrength() {
        final BranchSpecificSetup setup = buildBranchSpecificOuTreeLikelihood();
        final ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient traitGradient =
                new ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient(
                        dimTrait,
                        treeModel,
                        setup.likelihoodDelegate,
                        Arrays.asList(ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter.WRT_DIAGONAL_SELECTION_STRENGTH));
        final BranchSpecificGradient gradient =
                new BranchSpecificGradient("trait", setup.treeDataLikelihood, setup.likelihoodDelegate,
                        traitGradient, setup.attenuationParameter);

        final double[] analytic = gradient.getGradientLogDensity();
        final double[] numeric = numericalGradient(setup.treeDataLikelihood, setup.attenuationParameter, 1e-5);

        assertGradientClose("diagonal attenuation gradient", analytic, numeric, 1e-4, 2e-2);
    }

    public void testJointGradientDelegatesAttenuationGradient() {
        final BranchSpecificSetup setup = buildBranchSpecificOuTreeLikelihood();
        final ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient traitGradient =
                new ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient(
                        dimTrait,
                        treeModel,
                        setup.likelihoodDelegate,
                        Arrays.asList(ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter.WRT_DIAGONAL_SELECTION_STRENGTH));
        final BranchSpecificGradient branchSpecificGradient =
                new BranchSpecificGradient("trait", setup.treeDataLikelihood, setup.likelihoodDelegate,
                        traitGradient, setup.attenuationParameter);
        final AbstractDiffusionGradient.ParameterDiffusionGradient attenuationGradient =
                createDiagonalAttenuationGradient(branchSpecificGradient,
                        setup.treeDataLikelihood,
                        new DiagonalMatrix(setup.attenuationParameter));
        final JointGradient jointGradient =
                new JointGradient(Arrays.<dr.inference.hmc.GradientWrtParameterProvider>asList(attenuationGradient));

        final double[] joint = jointGradient.getGradientLogDensity();
        final double[] child = attenuationGradient.getGradientLogDensity();
        assertEquals("dimension", child.length, joint.length);
        for (int i = 0; i < child.length; ++i) {
            assertEquals("joint attenuation gradient idx=" + i, child[i], joint[i], 0.0);
        }
    }

    public void testTreeBridgePreservesOrthogonalBlockSelectionParameterization() {
        final int dimension = 4;
        final Parameter angles = new Parameter.Default(new double[]{0.2, -0.15, 0.1, 0.05, -0.08, 0.12});
        final GivensRotationMatrixParameter rotation =
                new GivensRotationMatrixParameter("R.block.tree", dimension, angles);
        final Parameter scalar = new Parameter.Default(0);
        final Parameter rho = new Parameter.Default(new double[]{0.8, 1.1});
        final Parameter theta = new Parameter.Default(new double[]{0.35, -0.2});
        final Parameter t = new Parameter.Default(new double[]{0.15, -0.1});
        final OrthogonalBlockDiagonalPolarStableMatrixParameter blockSelection =
                new OrthogonalBlockDiagonalPolarStableMatrixParameter(
                        "A.block.tree", rotation, scalar, rho, theta, t);

        final MultivariateElasticModel blockElastic = new MultivariateElasticModel(blockSelection);
        final MultivariateElasticModel denseElastic =
                new MultivariateElasticModel(makeMatrixParameter("A.dense.tree", blockSelection.getParameterAsMatrix()));
        final MultivariateDiffusionModel diffusion = new MultivariateDiffusionModel(
                makeMatrixParameter("Q.precision.tree", new double[][]{
                        {2.0, 0.1, 0.0, 0.0},
                        {0.1, 1.8, 0.05, 0.0},
                        {0.0, 0.05, 1.6, 0.08},
                        {0.0, 0.0, 0.08, 1.9}
                }));

        final TimeSeriesOUGaussianBranchTransitionProvider blockProvider =
                new TimeSeriesOUGaussianBranchTransitionProvider(blockElastic, diffusion);
        final TimeSeriesOUGaussianBranchTransitionProvider denseProvider =
                new TimeSeriesOUGaussianBranchTransitionProvider(denseElastic, diffusion);

        assertTrue("tree bridge should preserve orthogonal block parametrization",
                blockProvider.getProcessModel().getSelectionMatrixParameterization()
                        instanceof OrthogonalBlockDiagonalSelectionMatrixParameterization);

        final double dt = 0.17;
        final double[][] blockF = new double[dimension][dimension];
        final double[][] denseF = new double[dimension][dimension];
        final double[][] blockV = new double[dimension][dimension];
        final double[][] denseV = new double[dimension][dimension];

        blockProvider.fillBranchTransitionMatrix(dt, blockF);
        denseProvider.fillBranchTransitionMatrix(dt, denseF);
        blockProvider.fillBranchTransitionCovariance(dt, blockV);
        denseProvider.fillBranchTransitionCovariance(dt, denseV);

        assertMatricesClose(blockF, denseF, 1e-9);
        assertMatricesClose(blockV, denseV, 1e-9);
    }

    public void testCanonicalOrthogonalBlockLocalAdjointsMatchDenseTreeBridge() {
        final String previousUse = System.getProperty(USE_BRIDGE_PROPERTY);
        try {
            System.setProperty(USE_BRIDGE_PROPERTY, "true");
            final OrthogonalBlockTreeSetup setup = buildOrthogonalBlockTreeSetups();
            setup.block.treeDataLikelihood.getLogLikelihood();
            setup.dense.treeDataLikelihood.getLogLikelihood();

            final NodeRef node = getFirstInternalNonRootNode();
            final BranchSufficientStatistics blockStatistics = getBranchStatistics(setup.block, node);
            final BranchSufficientStatistics denseStatistics = getBranchStatistics(setup.dense, node);

            final CanonicalLocalTransitionAdjoints blockAdjoints = new CanonicalLocalTransitionAdjoints(dimTrait);
            final CanonicalLocalTransitionAdjoints denseAdjoints = new CanonicalLocalTransitionAdjoints(dimTrait);
            setup.block.delegate.fillCanonicalLocalAdjointsForBranch(blockStatistics, blockAdjoints);
            setup.dense.delegate.fillCanonicalLocalAdjointsForBranch(denseStatistics, denseAdjoints);

            assertVectorsClose(blockAdjoints.dLogL_dF, denseAdjoints.dLogL_dF, 1e-8);
            assertVectorsClose(blockAdjoints.dLogL_df, denseAdjoints.dLogL_df, 1e-8);
            assertVectorsClose(blockAdjoints.dLogL_dOmega, denseAdjoints.dLogL_dOmega, 1e-8);
        } finally {
            if (previousUse == null) {
                System.clearProperty(USE_BRIDGE_PROPERTY);
            } else {
                System.setProperty(USE_BRIDGE_PROPERTY, previousUse);
            }
        }
    }

    public void testCanonicalOrthogonalBlockSelectionGradientMatchesDenseTreeBridge() {
        final String previousUse = System.getProperty(USE_BRIDGE_PROPERTY);
        try {
            System.setProperty(USE_BRIDGE_PROPERTY, "true");
            final OrthogonalBlockTreeSetup setup = buildOrthogonalBlockTreeSetups();
            setup.block.treeDataLikelihood.getLogLikelihood();
            setup.dense.treeDataLikelihood.getLogLikelihood();

            final NodeRef node = getFirstInternalNonRootNode();
            final BranchSufficientStatistics blockStatistics = getBranchStatistics(setup.block, node);
            final BranchSufficientStatistics denseStatistics = getBranchStatistics(setup.dense, node);

            final double[] blockGradient = setup.block.delegate.getCanonicalGradientSelectionForBranch(
                    blockStatistics, node, setup.block.likelihoodDelegate.getIntegrator());
            final double[] denseGradient = setup.dense.delegate.getCanonicalGradientSelectionForBranch(
                    denseStatistics, node, setup.dense.likelihoodDelegate.getIntegrator());

            assertGradientClose("canonical block vs dense tree selection gradient",
                    blockGradient, denseGradient, 1e-8, 1e-8);
        } finally {
            if (previousUse == null) {
                System.clearProperty(USE_BRIDGE_PROPERTY);
            } else {
                System.setProperty(USE_BRIDGE_PROPERTY, previousUse);
            }
        }
    }

    public void testVarianceGradientPathIsSelectedFromSelectionParameterization() {
        final OrthogonalBlockTreeSetup setup = buildOrthogonalBlockTreeSetups();
        assertTrue("block OU should use the canonical variance-gradient bridge",
                setup.block.delegate.usesCanonicalVarianceGradient());
        assertFalse("dense OU should keep the legacy variance-gradient path",
                setup.dense.delegate.usesCanonicalVarianceGradient());
    }

    public void testCanonicalOrthogonalBlockGlobalNumericGradientIsFinite() {
        final String previousUse = System.getProperty(USE_BRIDGE_PROPERTY);
        try {
            System.setProperty(USE_BRIDGE_PROPERTY, "true");
            final OrthogonalBlockTreeSetup setup = buildOrthogonalBlockTreeSetups();
            final Parameter nativeParameter = setup.block.blockParameter.getParameter();

            final double[] numeric = numericalGradient(setup.block.treeDataLikelihood, nativeParameter, 1e-5);
            for (int i = 0; i < numeric.length; ++i) {
                assertTrue("numeric orthogonal block tree gradient idx=" + i + " must be finite",
                        Double.isFinite(numeric[i]));
            }
        } finally {
            if (previousUse == null) {
                System.clearProperty(USE_BRIDGE_PROPERTY);
            } else {
                System.setProperty(USE_BRIDGE_PROPERTY, previousUse);
            }
        }
    }

    public void testCanonicalOrthogonalBlockGlobalNumericGradientStressNearBoundaries() {
        final String previousUse = System.getProperty(USE_BRIDGE_PROPERTY);
        try {
            System.setProperty(USE_BRIDGE_PROPERTY, "true");
            final BoundaryCase[] cases = new BoundaryCase[]{
                    new BoundaryCase("small-rho-zero-theta", 0.2, 1.0e-3, 0.0, 1.0e-6, new double[]{0.0, 0.35, 1.56}),
                    new BoundaryCase("theta-near-edge", 0.4, 0.05, 1.45, 1.0e-3, new double[]{1.3, -1.2, 0.7}),
                    new BoundaryCase("negative-theta-edge", 1.1, 0.08, -1.45, -1.0e-3, new double[]{1.56, -0.9, 0.2})
            };

            for (final BoundaryCase boundaryCase : cases) {
                final OrthogonalBlockTreeSetup setup = buildOrthogonalBlockTreeSetups(
                        boundaryCase.scalar,
                        boundaryCase.rho,
                        boundaryCase.theta,
                        boundaryCase.t,
                        boundaryCase.angles);
                final Parameter nativeParameter = setup.block.blockParameter.getParameter();

                final double logLikelihood = setup.block.treeDataLikelihood.getLogLikelihood();
                assertTrue(boundaryCase.label + " log-likelihood must be finite", Double.isFinite(logLikelihood));

                final double[] numeric = numericalGradient(setup.block.treeDataLikelihood, nativeParameter, 5e-6);
                for (int i = 0; i < numeric.length; ++i) {
                    assertTrue(boundaryCase.label + " numeric gradient idx=" + i + " must be finite",
                            Double.isFinite(numeric[i]));
                }
            }
        } finally {
            if (previousUse == null) {
                System.clearProperty(USE_BRIDGE_PROPERTY);
            } else {
                System.setProperty(USE_BRIDGE_PROPERTY, previousUse);
            }
        }
    }

    private TreeDataLikelihood buildDiagonalOuTreeLikelihood() {
        return buildDiagonalOuTreeLikelihoodWithStrictOptima().treeDataLikelihood;
    }

    private BranchSpecificSetup buildDiagonalOuTreeLikelihoodWithStrictOptima() {
        final List<BranchRateModel> optimalTraitsModels = new ArrayList<BranchRateModel>();
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.1", new double[]{1.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.2", new double[]{2.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.3", new double[]{-2.0})));

        final Parameter attenuation = new Parameter.Default(new double[]{0.1, 100.0, 50.0});
        final DiagonalMatrix strengthOfSelectionMatrixParam =
                new DiagonalMatrix(attenuation);

        final DiffusionProcessDelegate diffusionProcessDelegate =
                new OUDiffusionModelDelegate(
                        treeModel,
                        diffusionModel,
                        optimalTraitsModels,
                        new MultivariateElasticModel(strengthOfSelectionMatrixParam));

        final ContinuousDataLikelihoodDelegate likelihoodDelegate = new ContinuousDataLikelihoodDelegate(
                treeModel,
                diffusionProcessDelegate,
                buildObservedDataModel(),
                rootPrior,
                rateTransformation,
                rateModel,
                false);

        return new BranchSpecificSetup(
                new TreeDataLikelihood(likelihoodDelegate, treeModel, rateModel),
                likelihoodDelegate,
                attenuation,
                null);
    }

    private BranchSpecificSetup buildBranchSpecificOuTreeLikelihood() {
        final int branchCount = treeModel.getNodeCount() - 1;
        final List<BranchRateModel> optimalTraitsModels = new ArrayList<BranchRateModel>();
        final CompoundParameter branchOptima = new CompoundParameter("branchOptima");
        final ArbitraryBranchRates.BranchRateTransform transform =
                ArbitraryBranchRates.make(false, true, false, null, null);

        for (int trait = 0; trait < dimTrait; ++trait) {
            final double[] values = new double[branchCount];
            for (int i = 0; i < branchCount; ++i) {
                values[i] = 0.4 * (trait + 1) - 0.1 * i;
            }
            final Parameter branchParameter = new Parameter.Default("optima." + trait, values);
            branchOptima.addParameter(branchParameter);
            optimalTraitsModels.add(new ArbitraryBranchRates(treeModel, branchParameter, transform, false));
        }

        final Parameter attenuation = new Parameter.Default(new double[]{0.25, 0.6, 0.9});
        final DiagonalMatrix strengthOfSelectionMatrixParam =
                new DiagonalMatrix(attenuation);

        final DiffusionProcessDelegate diffusionProcessDelegate =
                new OUDiffusionModelDelegate(
                        treeModel,
                        diffusionModel,
                        optimalTraitsModels,
                        new MultivariateElasticModel(strengthOfSelectionMatrixParam));

        final ContinuousDataLikelihoodDelegate likelihoodDelegate = new ContinuousDataLikelihoodDelegate(
                treeModel,
                diffusionProcessDelegate,
                buildObservedDataModel(),
                rootPrior,
                rateTransformation,
                rateModel,
                false);

        return new BranchSpecificSetup(
                new TreeDataLikelihood(likelihoodDelegate, treeModel, rateModel),
                likelihoodDelegate,
                attenuation,
                branchOptima);
    }

    private OrthogonalBlockTreeSetup buildOrthogonalBlockTreeSetups() {
        return buildOrthogonalBlockTreeSetups(
                1.3, 0.85, 0.25, -0.08, new double[]{0.2, -0.15, 0.1},
                buildObservedDataModel());
    }

    private OrthogonalBlockTreeSetup buildOrthogonalBlockTreeSetupsWithPartialTipMissing() {
        return buildOrthogonalBlockTreeSetups(
                1.3, 0.85, 0.25, -0.08, new double[]{0.2, -0.15, 0.1},
                buildPartiallyMissingDataModel());
    }

    private OrthogonalBlockTreeSetup buildOrthogonalBlockTreeSetups(final double scalarValue,
                                                                    final double rhoValue,
                                                                    final double thetaValue,
                                                                    final double tValue,
                                                                    final double[] angleValues) {
        return buildOrthogonalBlockTreeSetups(
                scalarValue, rhoValue, thetaValue, tValue, angleValues, buildObservedDataModel());
    }

    private OrthogonalBlockTreeSetup buildOrthogonalBlockTreeSetups(final double scalarValue,
                                                                    final double rhoValue,
                                                                    final double thetaValue,
                                                                    final double tValue,
                                                                    final double[] angleValues,
                                                                    final ContinuousTraitDataModel dataModel) {
        final int dimension = dimTrait;
        final Parameter angles = new Parameter.Default(angleValues);
        final GivensRotationMatrixParameter rotation =
                new GivensRotationMatrixParameter("R.block.tree", dimension, angles);
        final Parameter scalar = new Parameter.Default("scalar.block.tree", new double[]{scalarValue});
        final Parameter rho = new Parameter.Default("rho.block.tree", new double[]{rhoValue});
        final Parameter theta = new Parameter.Default("theta.block.tree", new double[]{thetaValue});
        final Parameter t = new Parameter.Default("t.block.tree", new double[]{tValue});
        final OrthogonalBlockDiagonalPolarStableMatrixParameter blockSelection =
                new OrthogonalBlockDiagonalPolarStableMatrixParameter(
                        "A.block.tree", rotation, scalar, rho, theta, t);

        final BranchSpecificSetup blockSetup = buildStrictOptimaOuTreeLikelihood(
                new MultivariateElasticModel(blockSelection), dataModel);
        final BranchSpecificSetup denseSetup = buildStrictOptimaOuTreeLikelihood(
                new MultivariateElasticModel(makeMatrixParameter("A.dense.block.tree",
                        blockSelection.getParameterAsMatrix())), dataModel);

        return new OrthogonalBlockTreeSetup(
                new TreeOuSetup(blockSetup.treeDataLikelihood,
                        blockSetup.likelihoodDelegate,
                        (OUDiffusionModelDelegate) blockSetup.likelihoodDelegate.getDiffusionProcessDelegate(),
                        blockSelection),
                new TreeOuSetup(denseSetup.treeDataLikelihood,
                        denseSetup.likelihoodDelegate,
                        (OUDiffusionModelDelegate) denseSetup.likelihoodDelegate.getDiffusionProcessDelegate(),
                        null));
    }

    private BranchSpecificSetup buildStrictOptimaOuTreeLikelihood(final MultivariateElasticModel elasticModel) {
        return buildStrictOptimaOuTreeLikelihood(elasticModel, buildObservedDataModel());
    }

    private BranchSpecificSetup buildStrictOptimaOuTreeLikelihood(final MultivariateElasticModel elasticModel,
                                                                  final ContinuousTraitDataModel dataModel) {
        final List<BranchRateModel> optimalTraitsModels = new ArrayList<BranchRateModel>();
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.1", new double[]{1.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.2", new double[]{2.0})));
        optimalTraitsModels.add(new StrictClockBranchRates(new Parameter.Default("rate.3", new double[]{-2.0})));

        final DiffusionProcessDelegate diffusionProcessDelegate =
                new OUDiffusionModelDelegate(
                        treeModel,
                        diffusionModel,
                        optimalTraitsModels,
                        elasticModel);

        final ContinuousDataLikelihoodDelegate likelihoodDelegate = new ContinuousDataLikelihoodDelegate(
                treeModel,
                diffusionProcessDelegate,
                dataModel,
                rootPrior,
                rateTransformation,
                rateModel,
                false);

        return new BranchSpecificSetup(
                new TreeDataLikelihood(likelihoodDelegate, treeModel, rateModel),
                likelihoodDelegate,
                null,
                null);
    }

    @SuppressWarnings("unchecked")
    private BranchSufficientStatistics getBranchStatistics(final TreeOuSetup setup,
                                                           final NodeRef node) {
        ensureBranchConditionalTrait(setup.likelihoodDelegate, "trait");
        final TreeTrait<List<BranchSufficientStatistics>> provider =
                setup.treeDataLikelihood.getTreeTrait(BranchConditionalDistributionDelegate.getName("trait"));
        assertNotNull("branch conditional trait provider", provider);
        final List<BranchSufficientStatistics> statistics = provider.getTrait(treeModel, node);
        assertNotNull("branch statistics", statistics);
        assertEquals("single trait expected", 1, statistics.size());
        return statistics.get(0);
    }

    private static void ensureBranchConditionalTrait(final ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                                     final String traitName) {
        try {
            final Method method = ContinuousDataLikelihoodDelegate.class
                    .getDeclaredMethod("addBranchConditionalDensityTrait", String.class);
            method.setAccessible(true);
            method.invoke(likelihoodDelegate, traitName);
        } catch (final Exception exception) {
            throw new RuntimeException("Failed to install branch conditional trait", exception);
        }
    }

    private NodeRef getFirstInternalNonRootNode() {
        for (int i = 0; i < treeModel.getNodeCount(); ++i) {
            final NodeRef node = treeModel.getNode(i);
            if (!treeModel.isRoot(node) && !treeModel.isExternal(node)) {
                return node;
            }
        }
        fail("expected an internal non-root node in the test tree");
        return null;
    }


    private static double[] numericalGradient(final TreeDataLikelihood likelihood,
                                              final Parameter parameter,
                                              final double step) {
        final double[] savedValues = parameter.getParameterValues();
        final MultivariateFunction function = new MultivariateFunction() {
            @Override
            public double evaluate(final double[] argument) {
                for (int i = 0; i < argument.length; ++i) {
                    parameter.setParameterValue(i, argument[i]);
                }
                likelihood.makeDirty();
                return likelihood.getLogLikelihood();
            }

            @Override
            public int getNumArguments() {
                return parameter.getDimension();
            }

            @Override
            public double getLowerBound(final int n) {
                return Double.NEGATIVE_INFINITY;
            }

            @Override
            public double getUpperBound(final int n) {
                return Double.POSITIVE_INFINITY;
            }
        };
        final double[] numeric = NumericalDerivative.gradient(function, savedValues.clone());
        for (int i = 0; i < savedValues.length; ++i) {
            parameter.setParameterValue(i, savedValues[i]);
        }
        likelihood.makeDirty();
        return numeric;
    }

    private static void assertGradientClose(final String label,
                                            final double[] analytic,
                                            final double[] numeric,
                                            final double absoluteTolerance,
                                            final double relativeTolerance) {
        assertEquals("dimension", analytic.length, numeric.length);
        for (int i = 0; i < analytic.length; ++i) {
            final double tolerance = Math.max(
                    absoluteTolerance,
                    relativeTolerance * Math.max(Math.abs(analytic[i]), Math.abs(numeric[i])));
            assertEquals(label + " idx=" + i, numeric[i], analytic[i], tolerance);
        }
    }

    private ContinuousTraitDataModel buildObservedDataModel() {
        return new ContinuousTraitDataModel(
                "observedDataModel",
                traitParameter,
                new boolean[traitParameter.getDimension()],
                true,
                dimTrait,
                PrecisionType.FULL);
    }

    private ContinuousTraitDataModel buildPartiallyMissingDataModel() {
        final boolean[] missingIndicators = new boolean[traitParameter.getDimension()];
        missingIndicators[3] = true;
        missingIndicators[4] = true;
        missingIndicators[5] = true;
        missingIndicators[7] = true;
        return new ContinuousTraitDataModel(
                "partiallyMissingDataModel",
                traitParameter,
                missingIndicators,
                true,
                dimTrait,
                PrecisionType.FULL);
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

    private static void assertMatricesClose(final double[][] actual,
                                            final double[][] expected,
                                            final double tolerance) {
        assertEquals("row dimension", expected.length, actual.length);
        for (int i = 0; i < expected.length; ++i) {
            assertEquals("column dimension", expected[i].length, actual[i].length);
            for (int j = 0; j < expected[i].length; ++j) {
                assertEquals("matrix[" + i + "][" + j + "]",
                        expected[i][j], actual[i][j], tolerance);
            }
        }
    }

    private static void assertVectorsClose(final double[] actual,
                                           final double[] expected,
                                           final double tolerance) {
        assertEquals("vector dimension", expected.length, actual.length);
        for (int i = 0; i < expected.length; ++i) {
            assertEquals("vector[" + i + "]", expected[i], actual[i], tolerance);
        }
    }

    private static final class BranchSpecificSetup {
        final TreeDataLikelihood treeDataLikelihood;
        final ContinuousDataLikelihoodDelegate likelihoodDelegate;
        final Parameter attenuationParameter;
        final CompoundParameter branchOptima;

        private BranchSpecificSetup(final TreeDataLikelihood treeDataLikelihood,
                                    final ContinuousDataLikelihoodDelegate likelihoodDelegate,
                                    final Parameter attenuationParameter,
                                    final CompoundParameter branchOptima) {
            this.treeDataLikelihood = treeDataLikelihood;
            this.likelihoodDelegate = likelihoodDelegate;
            this.attenuationParameter = attenuationParameter;
            this.branchOptima = branchOptima;
        }
    }

    private static final class TreeOuSetup {
        final TreeDataLikelihood treeDataLikelihood;
        final ContinuousDataLikelihoodDelegate likelihoodDelegate;
        final OUDiffusionModelDelegate delegate;
        final OrthogonalBlockDiagonalPolarStableMatrixParameter blockParameter;

        private TreeOuSetup(final TreeDataLikelihood treeDataLikelihood,
                            final ContinuousDataLikelihoodDelegate likelihoodDelegate,
                            final OUDiffusionModelDelegate delegate,
                            final OrthogonalBlockDiagonalPolarStableMatrixParameter blockParameter) {
            this.treeDataLikelihood = treeDataLikelihood;
            this.likelihoodDelegate = likelihoodDelegate;
            this.delegate = delegate;
            this.blockParameter = blockParameter;
        }
    }

    private static final class OrthogonalBlockTreeSetup {
        final TreeOuSetup block;
        final TreeOuSetup dense;

        private OrthogonalBlockTreeSetup(final TreeOuSetup block,
                                         final TreeOuSetup dense) {
            this.block = block;
            this.dense = dense;
        }
    }

    private static final class BoundaryCase {
        final String label;
        final double scalar;
        final double rho;
        final double theta;
        final double t;
        final double[] angles;

        private BoundaryCase(final String label,
                             final double scalar,
                             final double rho,
                             final double theta,
                             final double t,
                             final double[] angles) {
            this.label = label;
            this.scalar = scalar;
            this.rho = rho;
            this.theta = theta;
            this.t = t;
            this.angles = angles;
        }
    }

    public static Test suite() {
        return new TestSuite(OUDiffusionKernelBridgeValidationTest.class);
    }
}

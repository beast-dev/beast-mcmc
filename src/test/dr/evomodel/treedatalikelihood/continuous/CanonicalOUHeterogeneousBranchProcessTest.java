package test.dr.evomodel.treedatalikelihood.continuous;

import dr.evomodel.continuous.MultivariateElasticModel;
import dr.evomodel.treedatalikelihood.continuous.CanonicalOUBranchProcessSpecification;
import dr.evomodel.treedatalikelihood.continuous.canonical.adapter.CanonicalOUTreeLikelihoodIntegrator;
import dr.evomodel.treedatalikelihood.continuous.canonical.adapter.HomogeneousCanonicalOUBranchTransitionProvider;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.inference.model.DiagonalMatrix;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public final class CanonicalOUHeterogeneousBranchProcessTest
        extends CanonicalOUMessagePasserTestSupport {

    public CanonicalOUHeterogeneousBranchProcessTest(final String name) {
        super(name);
    }

    public void testProviderDispatchesBranchesToConfiguredOUProcess() {
        final MultivariateElasticModel slow = diagonalElasticModel("hetero.dispatch.slow", 0.7, 1.0, 1.4);
        final MultivariateElasticModel fast = diagonalElasticModel("hetero.dispatch.fast", 2.2, 2.5, 3.0);
        final int childUsingFastProcess = firstNonRootNode();
        final Parameter assignment =
                nodeAssignment("hetero.dispatch.assignment", childUsingFastProcess, 1);
        final CanonicalOUBranchProcessSpecification specification =
                specification(assignment, slow, fast);

        final Parameter stationaryMean =
                new Parameter.Default("hetero.dispatch.mu", MU_VALUES.clone());
        final HomogeneousCanonicalOUBranchTransitionProvider heterogeneous =
                new HomogeneousCanonicalOUBranchTransitionProvider(
                        treeModel,
                        specification,
                        diffusionModel,
                        stationaryMean,
                        rateModel,
                        rateTransformation);
        final HomogeneousCanonicalOUBranchTransitionProvider homogeneousFast =
                new HomogeneousCanonicalOUBranchTransitionProvider(
                        treeModel,
                        fast,
                        diffusionModel,
                        stationaryMean,
                        rateModel,
                        rateTransformation);

        assertEquals("configured branch process index",
                1, heterogeneous.getBranchProcessIndex(childUsingFastProcess));

        final CanonicalGaussianTransition actual = new CanonicalGaussianTransition(dimTrait);
        final CanonicalGaussianTransition expected = new CanonicalGaussianTransition(dimTrait);
        heterogeneous.fillCanonicalTransition(childUsingFastProcess, actual);
        homogeneousFast.fillCanonicalTransition(childUsingFastProcess, expected);

        assertTransitionEquals("heterogeneous branch transition", expected, actual);
    }

    public void testLikelihoodRespondsToBranchProcessAssignment() {
        final MultivariateElasticModel slow = diagonalElasticModel("hetero.like.slow", 0.7, 1.0, 1.4);
        final MultivariateElasticModel fast = diagonalElasticModel("hetero.like.fast", 2.2, 2.5, 3.0);
        final Parameter stationaryMean =
                new Parameter.Default("hetero.like.mu", MU_VALUES.clone());

        final CanonicalOUTreeLikelihoodIntegrator allSlow =
                new CanonicalOUTreeLikelihoodIntegrator(
                        treeModel,
                        specification(nodeAssignment("hetero.like.assignment.slow", -1, 0), slow, fast),
                        diffusionModel,
                        dataModel,
                        rootPrior,
                        stationaryMean,
                        rateModel,
                        rateTransformation);
        final CanonicalOUTreeLikelihoodIntegrator mixed =
                new CanonicalOUTreeLikelihoodIntegrator(
                        treeModel,
                        specification(nodeAssignment("hetero.like.assignment.mixed", firstNonRootNode(), 1), slow, fast),
                        diffusionModel,
                        dataModel,
                        rootPrior,
                        stationaryMean,
                        rateModel,
                        rateTransformation);

        final double allSlowLogLikelihood = allSlow.calculateLogLikelihood();
        final double mixedLogLikelihood = mixed.calculateLogLikelihood();

        assertTrue("branch-process assignment should change the likelihood",
                Math.abs(allSlowLogLikelihood - mixedLogLikelihood) > 1.0e-8);
    }

    public void testHeterogeneousProcessGradientsWorkForSharedSelectionParameter() {
        final MatrixParameter sharedSelection =
                new MatrixParameter("hetero.gradient.shared.A", dimTrait, dimTrait);
        fillMatrixParameter(sharedSelection, A_VALUES);
        final MultivariateElasticModel sharedElastic =
                new MultivariateElasticModel(sharedSelection);
        final CanonicalOUTreeLikelihoodIntegrator integrator =
                new CanonicalOUTreeLikelihoodIntegrator(
                        treeModel,
                        specification(
                                alternatingAssignment("hetero.gradient.shared.assignment"),
                                sharedElastic,
                                sharedElastic),
                        diffusionModel,
                        dataModel,
                        rootPrior,
                        new Parameter.Default("hetero.gradient.shared.mu", MU_VALUES.clone()),
                        rateModel,
                        rateTransformation);

        final double[] gradA = new double[dimTrait * dimTrait];
        final double[] gradQ = new double[dimTrait * dimTrait];
        final double[] gradMu = new double[dimTrait];
        integrator.computeSelectionGradient(gradA);
        integrator.computeDiffusionGradient(gradQ);
        integrator.computeStationaryMeanGradient(gradMu);

        assertFinite("selection gradient", gradA);
        assertFinite("diffusion gradient", gradQ);
        assertFinite("stationary mean gradient", gradMu);
    }

    public void testHeterogeneousIndependentSelectionGradientRejectsCleanly() {
        final MatrixParameter selection0 =
                new MatrixParameter("hetero.gradient.independent.A0", dimTrait, dimTrait);
        final MatrixParameter selection1 =
                new MatrixParameter("hetero.gradient.independent.A1", dimTrait, dimTrait);
        fillMatrixParameter(selection0, A_VALUES);
        fillMatrixParameter(selection1, new double[][]{
                {1.6, 0.2, 0.1},
                {0.2, 1.9, 0.0},
                {0.1, 0.0, 2.4}
        });

        final CanonicalOUTreeLikelihoodIntegrator integrator =
                new CanonicalOUTreeLikelihoodIntegrator(
                        treeModel,
                        specification(
                                alternatingAssignment("hetero.gradient.independent.assignment"),
                                new MultivariateElasticModel(selection0),
                                new MultivariateElasticModel(selection1)),
                        diffusionModel,
                        dataModel,
                        rootPrior,
                        new Parameter.Default("hetero.gradient.independent.mu", MU_VALUES.clone()),
                        rateModel,
                        rateTransformation);

        try {
            integrator.computeSelectionGradient(new double[dimTrait * dimTrait]);
            fail("Independent heterogeneous selection gradients should reject cleanly.");
        } catch (UnsupportedOperationException expected) {
            assertTrue(expected.getMessage().contains("heterogeneous processes require active branch processes"));
        }
    }

    private CanonicalOUBranchProcessSpecification specification(
            final Parameter assignment,
            final MultivariateElasticModel... elasticModels) {
        final List<MultivariateElasticModel> models =
                new ArrayList<MultivariateElasticModel>(Arrays.asList(elasticModels));
        return new CanonicalOUBranchProcessSpecification(models, assignment);
    }

    private Parameter nodeAssignment(final String id,
                                     final int selectedChildNode,
                                     final int selectedProcess) {
        final double[] values = new double[treeModel.getNodeCount()];
        for (int nodeIndex = 0; nodeIndex < values.length; ++nodeIndex) {
            values[nodeIndex] = 0.0;
        }
        if (selectedChildNode >= 0) {
            values[selectedChildNode] = selectedProcess;
        }
        return new Parameter.Default(id, values);
    }

    private Parameter alternatingAssignment(final String id) {
        final double[] values = new double[treeModel.getNodeCount()];
        final int rootIndex = treeModel.getRoot().getNumber();
        for (int nodeIndex = 0; nodeIndex < values.length; ++nodeIndex) {
            values[nodeIndex] = nodeIndex == rootIndex ? 0.0 : nodeIndex % 2;
        }
        return new Parameter.Default(id, values);
    }

    private int firstNonRootNode() {
        final int rootIndex = treeModel.getRoot().getNumber();
        for (int nodeIndex = 0; nodeIndex < treeModel.getNodeCount(); ++nodeIndex) {
            if (nodeIndex != rootIndex) {
                return nodeIndex;
            }
        }
        throw new IllegalStateException("Tree has no non-root branch.");
    }

    private MultivariateElasticModel diagonalElasticModel(final String id,
                                                          final double... values) {
        return new MultivariateElasticModel(
                new DiagonalMatrix(new Parameter.Default(id + ".diag", values)));
    }

    private void assertFinite(final String label, final double[] values) {
        for (int i = 0; i < values.length; ++i) {
            assertTrue(label + "[" + i + "] is finite", Double.isFinite(values[i]));
        }
    }

    private void assertTransitionEquals(final String label,
                                        final CanonicalGaussianTransition expected,
                                        final CanonicalGaussianTransition actual) {
        assertArrayEquals(label + " precisionXX", expected.precisionXX, actual.precisionXX);
        assertArrayEquals(label + " precisionXY", expected.precisionXY, actual.precisionXY);
        assertArrayEquals(label + " precisionYX", expected.precisionYX, actual.precisionYX);
        assertArrayEquals(label + " precisionYY", expected.precisionYY, actual.precisionYY);
        assertArrayEquals(label + " informationX", expected.informationX, actual.informationX);
        assertArrayEquals(label + " informationY", expected.informationY, actual.informationY);
        assertEquals(label + " logNormalizer", expected.logNormalizer, actual.logNormalizer, 1.0e-12);
    }

    private void assertArrayEquals(final String label,
                                   final double[] expected,
                                   final double[] actual) {
        assertEquals(label + " length", expected.length, actual.length);
        for (int i = 0; i < expected.length; ++i) {
            assertEquals(label + "[" + i + "]", expected[i], actual[i], 1.0e-12);
        }
    }
}

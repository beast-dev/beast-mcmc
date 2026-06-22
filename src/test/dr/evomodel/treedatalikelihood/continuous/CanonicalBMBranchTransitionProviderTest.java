package test.dr.evomodel.treedatalikelihood.continuous;

import dr.evolution.tree.SimpleNode;
import dr.evolution.tree.SimpleTree;
import dr.evomodel.branchratemodel.StrictClockBranchRates;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.continuous.ContinuousRateTransformation;
import dr.evomodel.treedatalikelihood.continuous.canonical.adapter.HomogeneousCanonicalBMBranchTransitionProvider;
import dr.evomodel.treedatalikelihood.continuous.canonical.math.GaussianFormConverter;
import dr.evomodel.treedatalikelihood.continuous.canonical.message.CanonicalGaussianTransition;
import dr.inference.model.MatrixParameter;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

public final class CanonicalBMBranchTransitionProviderTest extends TestCase {

    private static final double TOL = 1.0e-11;

    public void testHomogeneousProviderMatchesMomentConverter() {
        final int dim = 2;
        final TreeModel tree = createThreeTipTree();
        final MultivariateDiffusionModel diffusion = createDiffusion();
        final StrictClockBranchRates rates =
                new StrictClockBranchRates(new Parameter.Default(new double[] {2.0}));
        final ContinuousRateTransformation rateTransformation =
                new ContinuousRateTransformation.Default(tree, false, false);
        final HomogeneousCanonicalBMBranchTransitionProvider provider =
                new HomogeneousCanonicalBMBranchTransitionProvider(
                        tree,
                        diffusion,
                        rates,
                        rateTransformation);

        final int childNodeIndex = findNodeWithBranchLength(tree, 1.0);
        final double effectiveTime = provider.getEffectiveBranchLength(childNodeIndex);
        assertEquals("strict-clock effective branch time", 2.0, effectiveTime, TOL);
        assertNull("BM provider should not expose per-branch transition objects",
                provider.getCanonicalTransitionView(childNodeIndex));

        final CanonicalGaussianTransition direct = new CanonicalGaussianTransition(dim);
        provider.fillCanonicalTransition(childNodeIndex, direct);

        final double[] flatTransition = {
                1.0, 0.0,
                0.0, 1.0
        };
        final double[] offset = new double[dim];
        final double[] traitCovariance = new double[dim * dim];
        provider.fillTraitCovarianceFlat(traitCovariance);
        final double[] branchCovariance = traitCovariance.clone();
        for (int i = 0; i < branchCovariance.length; i++) {
            branchCovariance[i] *= effectiveTime;
        }

        final CanonicalGaussianTransition expected = new CanonicalGaussianTransition(dim);
        final GaussianFormConverter.Workspace workspace = new GaussianFormConverter.Workspace();
        workspace.ensureDim(dim);
        GaussianFormConverter.fillTransitionFromMoments(
                flatTransition,
                offset,
                branchCovariance,
                dim,
                workspace,
                expected);

        assertArrayClose("precisionXX", expected.precisionXX, direct.precisionXX);
        assertArrayClose("precisionXY", expected.precisionXY, direct.precisionXY);
        assertArrayClose("precisionYX", expected.precisionYX, direct.precisionYX);
        assertArrayClose("precisionYY", expected.precisionYY, direct.precisionYY);
        assertArrayClose("informationX", expected.informationX, direct.informationX);
        assertArrayClose("informationY", expected.informationY, direct.informationY);
        assertEquals("log normalizer", expected.logNormalizer, direct.logNormalizer, TOL);
    }

    public void testMomentProviderUsesCallerOwnedBuffers() {
        final int dim = 2;
        final TreeModel tree = createThreeTipTree();
        final HomogeneousCanonicalBMBranchTransitionProvider provider =
                new HomogeneousCanonicalBMBranchTransitionProvider(
                        tree,
                        createDiffusion(),
                        null,
                        new ContinuousRateTransformation.Default(tree, false, false));

        final double branchLength = 0.25;
        final double[] transition = {-1.0, -1.0, -1.0, -1.0};
        final double[] offset = {-1.0, -1.0};
        final double[] covariance = new double[dim * dim];
        provider.fillTransitionMomentsFlat(branchLength, transition, offset, covariance);

        assertArrayClose("identity transition",
                new double[] {1.0, 0.0, 0.0, 1.0},
                transition);
        assertArrayClose("zero offset", new double[] {0.0, 0.0}, offset);
        assertArrayClose("scaled covariance",
                new double[] {0.14285714285714285, -0.03571428571428571,
                        -0.03571428571428571, 0.07142857142857142},
                covariance);
    }

    private static MultivariateDiffusionModel createDiffusion() {
        final MatrixParameter precision = new MatrixParameter("bm.precision");
        precision.addParameter(new Parameter.Default(new double[] {2.0, 1.0}));
        precision.addParameter(new Parameter.Default(new double[] {1.0, 4.0}));
        return new MultivariateDiffusionModel(precision);
    }

    private static TreeModel createThreeTipTree() {
        final SimpleNode tip0 = new SimpleNode();
        final SimpleNode tip1 = new SimpleNode();
        final SimpleNode tip2 = new SimpleNode();
        final SimpleNode internal = new SimpleNode();
        internal.setHeight(0.5);
        internal.addChild(tip0);
        internal.addChild(tip1);

        final SimpleNode root = new SimpleNode();
        root.setHeight(1.0);
        root.addChild(internal);
        root.addChild(tip2);

        return new DefaultTreeModel(new SimpleTree(root));
    }

    private static int findNodeWithBranchLength(final TreeModel tree, final double branchLength) {
        for (int i = 0; i < tree.getNodeCount(); i++) {
            if (!tree.isRoot(tree.getNode(i))
                    && Math.abs(tree.getBranchLength(tree.getNode(i)) - branchLength) < TOL) {
                return i;
            }
        }
        throw new IllegalArgumentException("No node found with branch length " + branchLength);
    }

    private static void assertArrayClose(final String label, final double[] expected, final double[] actual) {
        assertEquals(label + " length", expected.length, actual.length);
        for (int i = 0; i < expected.length; i++) {
            assertEquals(label + "[" + i + "]", expected[i], actual[i], TOL);
        }
    }
}

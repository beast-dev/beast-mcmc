package test.dr.evomodel.continuous;

import dr.app.beast.BeastParser;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evomodel.continuous.GeneralOUTreeSimulator;
import dr.evomodel.continuous.MultivariateDiffusionModel;
import dr.evomodel.continuous.MultivariateElasticModel;
import dr.xml.XMLObject;
import dr.inference.model.MatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.inference.timeseries.gaussian.OUProcessModel;
import dr.inference.timeseries.representation.GaussianBranchTransitionKernel;
import dr.math.MathUtils;
import junit.framework.TestCase;
import org.ejml.data.DenseMatrix64F;
import org.ejml.ops.CommonOps;

import java.io.StringReader;
import java.util.Map;

public class GeneralOUTreeSimulatorTest extends TestCase {

    public void testTipMomentsMatchExactOUTransition() throws Exception {
        final Tree tree = new NewickImporter("(A:0.7,B:0.5);").importTree(null);

        final double[][] drift = {
                {0.8, 0.2},
                {-0.1, 1.1}
        };
        final double[][] precision = {
                {2.0, 0.3},
                {0.3, 1.5}
        };

        final MultivariateElasticModel elasticModel =
                new MultivariateElasticModel(matrixParameter("drift", drift));
        final MultivariateDiffusionModel diffusionModel =
                new MultivariateDiffusionModel(matrixParameter("precision", precision));
        final Parameter optimum = new Parameter.Default("optimum", new double[]{1.0, -0.5});
        final Parameter rootMean = new Parameter.Default("rootMean", new double[]{0.2, -0.4});
        final MatrixParameter rootCovariance = matrixParameter("rootCovariance", new double[][]{
                {0.0, 0.0},
                {0.0, 0.0}
        });

        final GeneralOUTreeSimulator simulator =
                new GeneralOUTreeSimulator(tree, elasticModel, diffusionModel, optimum, rootMean, rootCovariance);

        final MatrixParameter diffusionCovariance = matrixParameter("diffusionCovariance", invert(precision));
        final MatrixParameter zeroInitialCovariance = matrixParameter("initialCovariance", new double[][]{
                {0.0, 0.0},
                {0.0, 0.0}
        });
        final OUProcessModel processModel = new OUProcessModel(
                "test.ou",
                2,
                elasticModel.getStrengthOfSelectionMatrixParameter(),
                diffusionCovariance,
                optimum,
                zeroInitialCovariance);
        final GaussianBranchTransitionKernel kernel =
                processModel.getRepresentation(GaussianBranchTransitionKernel.class);

        final double[][] transitionMatrix = new double[2][2];
        final double[] transitionOffset = new double[2];
        final double[][] transitionCovariance = new double[2][2];
        kernel.fillTransitionMatrix(0.7, transitionMatrix);
        kernel.fillTransitionOffset(0.7, transitionOffset);
        kernel.fillTransitionCovariance(0.7, transitionCovariance);

        final double[] expectedMean = multiply(transitionMatrix, rootMean.getParameterValues());
        addInPlace(expectedMean, transitionOffset);

        final NodeRef tipA = findExternalNode(tree, "A");
        final int replicates = 6000;
        final double[] empiricalMean = new double[2];
        final double[][] empiricalSecondMoment = new double[2][2];

        MathUtils.setSeed(17L);
        for (int i = 0; i < replicates; ++i) {
            final double[][] nodeTraits = simulator.simulateNodeTraits();
            final double[] sample = nodeTraits[tipA.getNumber()];
            empiricalMean[0] += sample[0];
            empiricalMean[1] += sample[1];
            empiricalSecondMoment[0][0] += sample[0] * sample[0];
            empiricalSecondMoment[0][1] += sample[0] * sample[1];
            empiricalSecondMoment[1][0] += sample[1] * sample[0];
            empiricalSecondMoment[1][1] += sample[1] * sample[1];
        }

        for (int i = 0; i < 2; ++i) {
            empiricalMean[i] /= replicates;
        }

        final double[][] empiricalCovariance = new double[2][2];
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 2; ++j) {
                empiricalCovariance[i][j] =
                        empiricalSecondMoment[i][j] / replicates - empiricalMean[i] * empiricalMean[j];
            }
        }

        assertEquals(expectedMean[0], empiricalMean[0], 0.03);
        assertEquals(expectedMean[1], empiricalMean[1], 0.03);
        assertEquals(transitionCovariance[0][0], empiricalCovariance[0][0], 0.04);
        assertEquals(transitionCovariance[0][1], empiricalCovariance[0][1], 0.04);
        assertEquals(transitionCovariance[1][0], empiricalCovariance[1][0], 0.04);
        assertEquals(transitionCovariance[1][1], empiricalCovariance[1][1], 0.04);
    }

    public void testAnnotatedTreeCarriesNodeAndTaxonTraits() throws Exception {
        final Tree tree = new NewickImporter("(A:0.7,B:0.5);").importTree(null);

        final double[][] drift = {
                {0.8, 0.0},
                {0.0, 1.1}
        };
        final double[][] precision = {
                {2.0, 0.0},
                {0.0, 1.5}
        };

        final GeneralOUTreeSimulator simulator = new GeneralOUTreeSimulator(
                tree,
                new MultivariateElasticModel(matrixParameter("drift", drift)),
                new MultivariateDiffusionModel(matrixParameter("precision", precision)),
                new Parameter.Default("optimum", new double[]{1.0, -0.5}),
                new Parameter.Default("rootMean", new double[]{0.2, -0.4}),
                matrixParameter("rootCovariance", new double[][]{{0.0, 0.0}, {0.0, 0.0}}));

        MathUtils.setSeed(9L);
        final Tree annotated = simulator.simulateAnnotatedTree("ouTrait", true);

        final NodeRef root = annotated.getRoot();
        final Object rootTrait = annotated.getNodeAttribute(root, "ouTrait");
        assertTrue(rootTrait instanceof double[]);
        assertEquals(2, ((double[]) rootTrait).length);

        final NodeRef tipA = findExternalNode(annotated, "A");
        final Object taxonTrait = annotated.getNodeTaxon(tipA).getAttribute("ouTrait");
        assertTrue(taxonTrait instanceof String);
        assertTrue(((String) taxonTrait).contains(" "));
    }

    public void testXmlParserSimulatesAnnotatedTree() throws Exception {
        final String previousParsers = System.getProperty("parsers");
        System.setProperty("parsers", "development");
        try {
            final BeastParser parser = new BeastParser(new String[0], null, false, true, true, null);
            final String xml =
                    "<beast>" +
                            "  <taxa id=\"taxa\">" +
                            "    <taxon id=\"A\"/>" +
                            "    <taxon id=\"B\"/>" +
                            "  </taxa>" +
                            "  <generalOuTreeSimulator id=\"simTree\" traitName=\"ouTrait\" clone=\"true\" selectionChart=\"dense\">" +
                            "    <newick id=\"startTree\" usingDates=\"false\" usingHeights=\"false\">(A:0.7,B:0.5);</newick>" +
                            "    <multivariateDiffusionModel id=\"diffusion\">" +
                            "      <precisionMatrix>" +
                            "        <matrixParameter id=\"precision\">" +
                            "          <parameter value=\"2.0 0.3\"/>" +
                            "          <parameter value=\"0.3 1.5\"/>" +
                            "        </matrixParameter>" +
                            "      </precisionMatrix>" +
                            "    </multivariateDiffusionModel>" +
                            "    <strengthOfSelectionMatrix>" +
                            "      <matrixParameter id=\"selection\">" +
                            "        <parameter value=\"0.8 -0.1\"/>" +
                            "        <parameter value=\"0.2 1.1\"/>" +
                            "      </matrixParameter>" +
                            "    </strengthOfSelectionMatrix>" +
                            "    <optimalTraits>" +
                            "      <parameter id=\"optimum\" value=\"1.0 -0.5\"/>" +
                            "    </optimalTraits>" +
                            "    <rootMean>" +
                            "      <parameter id=\"rootMean\" value=\"0.2 -0.4\"/>" +
                            "    </rootMean>" +
                            "    <rootCovariance>" +
                            "      <matrixParameter id=\"rootCovariance\">" +
                            "        <parameter value=\"0.0 0.0\"/>" +
                            "        <parameter value=\"0.0 0.0\"/>" +
                            "      </matrixParameter>" +
                            "    </rootCovariance>" +
                            "  </generalOuTreeSimulator>" +
                            "</beast>";

            final Map<String, XMLObject> objects = parser.parse(new StringReader(xml), true);
            final Tree simulated = (Tree) objects.get("simTree").getNativeObject();
            final NodeRef tipA = findExternalNode(simulated, "A");

            final Object nodeTrait = simulated.getNodeAttribute(tipA, "ouTrait");
            assertTrue(nodeTrait instanceof double[]);
            assertEquals(2, ((double[]) nodeTrait).length);

            final Object taxonTrait = simulated.getNodeTaxon(tipA).getAttribute("ouTrait");
            assertTrue(taxonTrait instanceof String);
        } finally {
            if (previousParsers == null) {
                System.clearProperty("parsers");
            } else {
                System.setProperty("parsers", previousParsers);
            }
        }
    }

    public void testXmlParserRejectsDenseSelectionWithoutExplicitChart() throws Exception {
        final String previousParsers = System.getProperty("parsers");
        System.setProperty("parsers", "development");
        try {
            final BeastParser parser = new BeastParser(new String[0], null, false, true, true, null);
            final String xml =
                    "<beast>" +
                            "  <taxa id=\"taxa\">" +
                            "    <taxon id=\"A\"/>" +
                            "    <taxon id=\"B\"/>" +
                            "  </taxa>" +
                            "  <generalOuTreeSimulator id=\"simTree\" traitName=\"ouTrait\" clone=\"true\">" +
                            "    <newick id=\"startTree\" usingDates=\"false\" usingHeights=\"false\">(A:0.7,B:0.5);</newick>" +
                            "    <multivariateDiffusionModel id=\"diffusion\">" +
                            "      <precisionMatrix>" +
                            "        <matrixParameter id=\"precision\">" +
                            "          <parameter value=\"2.0 0.3\"/>" +
                            "          <parameter value=\"0.3 1.5\"/>" +
                            "        </matrixParameter>" +
                            "      </precisionMatrix>" +
                            "    </multivariateDiffusionModel>" +
                            "    <strengthOfSelectionMatrix>" +
                            "      <matrixParameter id=\"selection\">" +
                            "        <parameter value=\"0.8 -0.1\"/>" +
                            "        <parameter value=\"0.2 1.1\"/>" +
                            "      </matrixParameter>" +
                            "    </strengthOfSelectionMatrix>" +
                            "    <optimalTraits>" +
                            "      <parameter id=\"optimum\" value=\"1.0 -0.5\"/>" +
                            "    </optimalTraits>" +
                            "    <rootMean>" +
                            "      <parameter id=\"rootMean\" value=\"0.2 -0.4\"/>" +
                            "    </rootMean>" +
                            "    <rootCovariance>" +
                            "      <matrixParameter id=\"rootCovariance\">" +
                            "        <parameter value=\"0.0 0.0\"/>" +
                            "        <parameter value=\"0.0 0.0\"/>" +
                            "      </matrixParameter>" +
                            "    </rootCovariance>" +
                            "  </generalOuTreeSimulator>" +
                            "</beast>";

            try {
                parser.parse(new StringReader(xml), true);
                fail("Expected XMLParseException for dense selection matrix without explicit selectionChart");
            } catch (Exception expected) {
                assertTrue(expected.getMessage().contains("selectionChart=\"dense\""));
            }
        } finally {
            if (previousParsers == null) {
                System.clearProperty("parsers");
            } else {
                System.setProperty("parsers", previousParsers);
            }
        }
    }

    private static MatrixParameter matrixParameter(final String name, final double[][] values) {
        final Parameter[] columns = new Parameter[values[0].length];
        for (int j = 0; j < values[0].length; ++j) {
            final double[] column = new double[values.length];
            for (int i = 0; i < values.length; ++i) {
                column[i] = values[i][j];
            }
            columns[j] = new Parameter.Default(name + ".col" + j, column);
        }
        return new MatrixParameter(name, columns);
    }

    private static double[] multiply(final double[][] matrix, final double[] vector) {
        final double[] out = new double[matrix.length];
        for (int i = 0; i < matrix.length; ++i) {
            double value = 0.0;
            for (int j = 0; j < vector.length; ++j) {
                value += matrix[i][j] * vector[j];
            }
            out[i] = value;
        }
        return out;
    }

    private static void addInPlace(final double[] target, final double[] increment) {
        for (int i = 0; i < target.length; ++i) {
            target[i] += increment[i];
        }
    }

    private static NodeRef findExternalNode(final Tree tree, final String taxonId) {
        for (int i = 0; i < tree.getExternalNodeCount(); ++i) {
            final NodeRef node = tree.getExternalNode(i);
            if (taxonId.equals(tree.getNodeTaxon(node).getId())) {
                return node;
            }
        }
        throw new IllegalArgumentException("Taxon not found: " + taxonId);
    }

    private static double[][] invert(final double[][] matrix) {
        final DenseMatrix64F dense = new DenseMatrix64F(matrix.length, matrix.length);
        for (int i = 0; i < matrix.length; ++i) {
            for (int j = 0; j < matrix.length; ++j) {
                dense.set(i, j, matrix[i][j]);
            }
        }
        CommonOps.invert(dense);
        final double[][] out = new double[matrix.length][matrix.length];
        for (int i = 0; i < matrix.length; ++i) {
            for (int j = 0; j < matrix.length; ++j) {
                out[i][j] = dense.get(i, j);
            }
        }
        return out;
    }
}

/*
 * MascotGradientTest.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package test.dr.evomodel.coalescent;

import dr.evolution.alignment.Patterns;
import dr.evolution.datatype.GeneralDataType;
import dr.evolution.io.NewickImporter;
import dr.evomodel.coalescent.basta.ConstantPopulationSizeModel;
import dr.evomodel.coalescent.basta.StructuredCoalescentLikelihoodGradient;
import dr.evomodel.coalescent.mascot.MascotLikelihood;
import dr.evomodel.substmodel.ComplexSubstitutionModel;
import dr.evomodel.substmodel.ComplexSubstitutionModelAtStationarity;
import dr.evomodel.substmodel.ComplexSubstitutionModelGradientSupport;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SVSComplexSubstitutionModel;
import dr.evomodel.substmodel.SubstitutionModel;
import dr.evomodel.tree.DefaultTreeModel;
import dr.evomodel.tree.TreeModel;
import dr.inference.hmc.CompactGradient;
import dr.inference.model.DuplicatedParameter;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

/** MASCOT adjoint gradients against central finite differences, including at high coalescence pressure. */
public class MascotGradientTest extends TestCase {

    private static final double RELATIVE_STEP = 1.0e-4;

    // State of tests/TestXML/testMascotHmcGradientCheck.xml.
    private static final String MODERATE_TREE = "(((A:1.0,B:1.0):1.0,(C:1.0,D:1.0):1.0):1.0,(E:2.0,F:2.0):1.0);";

    // State of tests/TestXML/testMascotHighPressureLogging.xml: the (E,F) branch integrates a coalescent
    // pressure of about 20, which amplified roundoff off the probability simplex by about e^20 before M0.
    private static final String HIGH_PRESSURE_TREE =
            "(((A:2.8877280032784474,B:2.8877280032784474):1.5358558603633852," +
                    "(C:3.3883756235523608,D:3.3883756235523608):1.0352082400894718):1.4805952691315536," +
                    "(E:0.5493739204766552,F:0.5493739204766552):5.354805212296731);";
    private static final double[] HIGH_PRESSURE_RATES = {0.022840096758276696, 0.9991617072287613,
            0.6039584418543463, 0.006329302680772018, 0.6532680506509443, 0.004136765472939712};
    private static final double[] HIGH_PRESSURE_POP_SIZES = {0.15729183324413404, 0.20647504338270894,
            0.03191953596439705};
    private static final double HIGH_PRESSURE_CONVERGED_LOG_LIKELIHOOD = -73.73126497724;

    public void testModerateGradientsMatchFiniteDifferences() throws Exception {
        GeneralDataType dataType = dataType(2);
        Parameter migrationRates = new Parameter.Default("theta", new double[]{0.2, 0.3});
        Parameter popSizes = new Parameter.Default("popSizes", new double[]{1.5, 1.8});
        MascotLikelihood likelihood = mascot(MODERATE_TREE, new int[]{0, 1, 0, 1, 0, 1}, dataType,
                migrationRates, null, popSizes, 0.005);

        assertGradientsMatch(likelihood, migrationRates, popSizes, 1.0e-6);
    }

    public void testHighPressureGradientsMatchFiniteDifferences() throws Exception {
        GeneralDataType dataType = dataType(3);
        Parameter rates = new Parameter.Default("location.rates", HIGH_PRESSURE_RATES);
        Parameter popSizes = new Parameter.Default("popSizes", HIGH_PRESSURE_POP_SIZES);
        MascotLikelihood likelihood = highPressure(dataType, rates, popSizes, 0.01);

        assertGradientsMatch(likelihood, rates, popSizes, 1.0e-5);
    }

    /**
     * A rate whose BSSVS indicator is 0 does not enter the likelihood, so its gradient must be exactly 0. Index 3
     * (1 to 0) is where the rates-parameter order and the row-major migration order differ.
     */
    public void testBssvsGradientsMatchFiniteDifferences() throws Exception {
        GeneralDataType dataType = dataType(3);
        Parameter rates = new Parameter.Default("location.rates", HIGH_PRESSURE_RATES);
        Parameter indicators = new Parameter.Default("location.indicators", new double[]{1.0, 0.0, 1.0, 0.0, 1.0, 1.0});
        Parameter popSizes = new Parameter.Default("popSizes", HIGH_PRESSURE_POP_SIZES);
        FrequencyModel frequencies = new FrequencyModel(dataType, new double[]{0.2, 0.3, 0.5});
        SVSComplexSubstitutionModel migrationModel =
                new SVSComplexSubstitutionModel("location.model", dataType, frequencies, rates, indicators);
        migrationModel.setNormalization(false);
        migrationModel.setScaleRatesByFrequencies(true);
        MascotLikelihood likelihood = mascot(HIGH_PRESSURE_TREE, new int[]{0, 1, 2, 0, 1, 2}, dataType, null,
                migrationModel, popSizes, 0.01);

        double[] migrationGradient = likelihood.getMigrationGradientLogDensity();
        assertEquals(0.0, migrationGradient[1], 0.0);
        assertEquals(0.0, migrationGradient[3], 0.0);
        assertGradientsMatch(likelihood, rates, popSizes, 1.0e-5);
    }

    /**
     * A normalized model divides Q by s = sum_i pi_i sum_{j != i} Q_ij, so every rate also acts through s. Unequal
     * frequencies make pi_i matter, with and without frequency scaling of the rates.
     */
    public void testNormalizedGradientsMatchFiniteDifferences() throws Exception {
        for (boolean scaleByFrequencies : new boolean[]{false, true}) {
            GeneralDataType dataType = dataType(3);
            Parameter rates = new Parameter.Default("location.rates", HIGH_PRESSURE_RATES);
            Parameter popSizes = new Parameter.Default("popSizes", HIGH_PRESSURE_POP_SIZES);
            FrequencyModel frequencies = new FrequencyModel(dataType, new double[]{0.2, 0.3, 0.5});
            ComplexSubstitutionModel migrationModel =
                    new ComplexSubstitutionModel("location.model", dataType, frequencies, rates);
            migrationModel.setNormalization(true);
            migrationModel.setScaleRatesByFrequencies(scaleByFrequencies);
            assertNull(ComplexSubstitutionModelGradientSupport.getCompatibilityError(migrationModel, "MASCOT", true));
            assertNotNull(ComplexSubstitutionModelGradientSupport.getCompatibilityError(migrationModel, "BASTA"));
            MascotLikelihood likelihood = mascot(HIGH_PRESSURE_TREE, new int[]{0, 1, 2, 0, 1, 2}, dataType, null,
                    migrationModel, popSizes, 0.01);

            assertGradientsMatch(likelihood, rates, popSizes, 1.0e-5);
        }
    }

    /**
     * The XML parser builds anonymous subclasses that only choose the eigen system; those are supported. A model that
     * normalizes with Q's own stationary distribution is not, because its pi depends on the rates.
     */
    public void testNormalizationSupportCoversParserModelsOnly() {
        GeneralDataType dataType = dataType(3);
        Parameter rates = new Parameter.Default("location.rates", HIGH_PRESSURE_RATES);
        FrequencyModel frequencies = new FrequencyModel(dataType, new double[]{0.2, 0.3, 0.5});
        ComplexSubstitutionModel parserLike = new ComplexSubstitutionModel("parserLike", dataType, frequencies, rates) {
        };
        parserLike.setNormalization(true);
        assertNull(ComplexSubstitutionModelGradientSupport.getCompatibilityError(parserLike, "MASCOT", true));

        ComplexSubstitutionModelAtStationarity stationary =
                new ComplexSubstitutionModelAtStationarity("stationary", dataType, rates);
        stationary.setNormalization(true);
        assertNotNull(ComplexSubstitutionModelGradientSupport.getCompatibilityError(stationary, "MASCOT", true));
    }

    /** Normalization with BSSVS: an inactive rate enters neither Q nor s, so its gradient is still exactly 0. */
    public void testNormalizedBssvsGradientsMatchFiniteDifferences() throws Exception {
        GeneralDataType dataType = dataType(3);
        Parameter rates = new Parameter.Default("location.rates", HIGH_PRESSURE_RATES);
        Parameter indicators = new Parameter.Default("location.indicators", new double[]{1.0, 0.0, 1.0, 0.0, 1.0, 1.0});
        Parameter popSizes = new Parameter.Default("popSizes", HIGH_PRESSURE_POP_SIZES);
        FrequencyModel frequencies = new FrequencyModel(dataType, new double[]{0.2, 0.3, 0.5});
        SVSComplexSubstitutionModel migrationModel =
                new SVSComplexSubstitutionModel("location.model", dataType, frequencies, rates, indicators);
        migrationModel.setNormalization(true);
        migrationModel.setScaleRatesByFrequencies(true);
        MascotLikelihood likelihood = mascot(HIGH_PRESSURE_TREE, new int[]{0, 1, 2, 0, 1, 2}, dataType, null,
                migrationModel, popSizes, 0.01);

        double[] migrationGradient = likelihood.getMigrationGradientLogDensity();
        assertEquals(0.0, migrationGradient[1], 0.0);
        assertEquals(0.0, migrationGradient[3], 0.0);
        assertGradientsMatch(likelihood, rates, popSizes, 1.0e-5);
    }

    /** One Ne shared by all demes through a DuplicatedParameter: its gradient is the sum over the copies. */
    public void testCompactedDuplicatedPopSizeGradientMatchesFiniteDifferences() throws Exception {
        GeneralDataType dataType = dataType(3);
        Parameter popSize = new Parameter.Default("popSizes", 0.2);
        popSize.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, 0.0, 1));
        DuplicatedParameter popSizes = new DuplicatedParameter(popSize);
        popSizes.addDuplicationParameter(new Parameter.Default(3.0));
        MascotLikelihood likelihood = highPressure(dataType,
                new Parameter.Default("location.rates", HIGH_PRESSURE_RATES), popSizes, 0.01);

        CompactGradient gradient = new CompactGradient(new StructuredCoalescentLikelihoodGradient(likelihood, null,
                StructuredCoalescentLikelihoodGradient.WrtParameter.POPULATION_SIZE));
        double[] perCopy = likelihood.getPopSizeGradientLogDensity();
        double[] compacted = gradient.getGradientLogDensity();

        assertEquals(3, perCopy.length);
        assertEquals(1, compacted.length);
        assertEquals(perCopy[0] + perCopy[1] + perCopy[2], compacted[0], 0.0);
        assertMatchesFiniteDifferences("popSize", compacted, likelihood, popSize, 1.0e-5);
    }

    public void testHighPressureLogLikelihoodAndScores() throws Exception {
        GeneralDataType dataType = dataType(3);
        MascotLikelihood likelihood = highPressure(dataType,
                new Parameter.Default("location.rates", HIGH_PRESSURE_RATES),
                new Parameter.Default("popSizes", HIGH_PRESSURE_POP_SIZES), 0.001);

        assertEquals(HIGH_PRESSURE_CONVERGED_LOG_LIKELIHOOD, likelihood.getLogLikelihood(), 1.0e-7);

        double[] scores = likelihood.getAncestralStateScores();
        TreeModel tree = likelihood.getTreeModel();
        for (int i = 0; i < tree.getInternalNodeCount(); i++) {
            int offset = tree.getInternalNode(i).getNumber() * 3;
            double sum = 0.0;
            for (int s = 0; s < 3; s++) {
                assertTrue("non-finite score at internal node " + i, Double.isFinite(scores[offset + s]));
                sum += scores[offset + s];
            }
            assertEquals("score row sum at internal node " + i, 1.0, sum, 1.0e-12);
        }
    }

    private static void assertGradientsMatch(MascotLikelihood likelihood, Parameter migrationRates,
                                             Parameter popSizes, double tolerance) {
        double[] migrationGradient = likelihood.getMigrationGradientLogDensity();
        double[] popSizeGradient = likelihood.getPopSizeGradientLogDensity();
        assertMatchesFiniteDifferences("migration", migrationGradient, likelihood, migrationRates, tolerance);
        assertMatchesFiniteDifferences("popSize", popSizeGradient, likelihood, popSizes, tolerance);
    }

    private static void assertMatchesFiniteDifferences(String label, double[] analytic, MascotLikelihood likelihood,
                                                       Parameter parameter, double tolerance) {
        assertEquals(parameter.getDimension(), analytic.length);
        for (int i = 0; i < analytic.length; i++) {
            double x = parameter.getParameterValue(i);
            double xUp = x * (1.0 + RELATIVE_STEP);
            double xDown = x * (1.0 - RELATIVE_STEP);
            parameter.setParameterValue(i, xUp);
            double up = likelihood.getLogLikelihood();
            parameter.setParameterValue(i, xDown);
            double down = likelihood.getLogLikelihood();
            parameter.setParameterValue(i, x);

            double numeric = (up - down) / (xUp - xDown);
            if (numeric == 0.0) {
                assertEquals(label + "[" + i + "]", 0.0, analytic[i], 0.0);
                continue;
            }
            double relativeError = Math.abs(analytic[i] - numeric) / Math.abs(numeric);
            assertTrue(label + "[" + i + "]: analytic " + analytic[i] + ", numeric " + numeric +
                    ", relative error " + relativeError, relativeError < tolerance);
        }
    }

    private static MascotLikelihood highPressure(GeneralDataType dataType, Parameter rates, Parameter popSizes,
                                                 double maxStep) throws Exception {
        FrequencyModel frequencies = new FrequencyModel(dataType, new double[]{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0});
        ComplexSubstitutionModel migrationModel =
                new ComplexSubstitutionModel("location.model", dataType, frequencies, rates);
        migrationModel.setNormalization(false);
        migrationModel.setScaleRatesByFrequencies(false);
        return mascot(HIGH_PRESSURE_TREE, new int[]{0, 1, 2, 0, 1, 2}, dataType, null, migrationModel, popSizes,
                maxStep);
    }

    /** Tip states are indexed by taxon name: A is tipStates[0], B is tipStates[1], and so on. */
    private static MascotLikelihood mascot(String newick, int[] tipStates, GeneralDataType dataType,
                                           Parameter migrationRates, SubstitutionModel migrationModel,
                                           Parameter popSizes, double maxStep) throws Exception {
        int stateCount = dataType.getStateCount();
        TreeModel tree = new DefaultTreeModel("treeModel", new NewickImporter(newick).importTree(null));
        Patterns tipPatterns = new Patterns(dataType, tree);
        int[] pattern = new int[tree.getTaxonCount()];
        for (int i = 0; i < pattern.length; i++) {
            pattern[i] = tipStates[tree.getTaxonId(i).charAt(0) - 'A'];
        }
        tipPatterns.addPattern(pattern);

        return new MascotLikelihood(null, tree, tipPatterns, migrationRates,
                migrationModel == null ? null : new SubstitutionModel[]{migrationModel}, null,
                new ConstantPopulationSizeModel("constant", popSizes, stateCount, -1), stateCount, null, maxStep,
                false, MascotLikelihood.DEFAULT_ANCESTRAL_STATE_TAG_NAME, MascotLikelihood.DEFAULT_MODE_STATE_TAG_NAME);
    }

    private static GeneralDataType dataType(int stateCount) {
        String[] codes = new String[stateCount];
        for (int i = 0; i < stateCount; i++) {
            codes[i] = Integer.toString(i);
        }
        return new GeneralDataType(codes);
    }
}

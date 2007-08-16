package dr.evomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.evolution.datatype.TwoStateCovarion;
import dr.inference.model.Parameter;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * TwoStateCovarionModel Tester.
 *
 * @author Alexei Drummond
 */
public class TwoStateCovarionModelTest extends TestCase {
    public TwoStateCovarionModelTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();

        frequencies = new Parameter.Default(new double[]{0.25, 0.25, 0.25, 0.25});
        alpha = new Parameter.Default(0.1);
        switchingRate = new Parameter.Default(0.5);

        FrequencyModel freqModel = new FrequencyModel(TwoStateCovarion.INSTANCE, frequencies);
        model = new TwoStateCovarionModel(TwoStateCovarion.INSTANCE, freqModel, alpha, switchingRate);
        dataType = model.dataType;
    }

    public void testTransitionProbabilities() {

        // with alpha == 1, the transition probability should be the same as binary jukes cantor
        alpha.setParameterValue(0, 1.0);
        switchingRate.setParameterValue(0, 1.0);

        model.setupMatrix();

        double[] matrix = new double[16];

        double[] pi = model.getFrequencyModel().getFrequencies();

        for (double distance = 0.01; distance <= 1.005; distance += 0.01) {
            model.getTransitionProbabilities(distance, matrix);

            double pChange =
                    (matrix[1] + matrix[3]) * pi[0] +
                            (matrix[4] + matrix[6]) * pi[1] +
                            (matrix[9] + matrix[11]) * pi[2] +
                            (matrix[12] + matrix[14]) * pi[3];

            // analytical result for the probability of a mismatch in binary jukes cantor model
            double jc = 0.5 * (1 - Math.exp(-2.0 * distance));

            assertEquals(pChange, jc, 1e-14);
        }

        // with alpha < 1.0, the probability of difference should be smaller than binary jukes cantor
        alpha.setParameterValue(0, 0.0);

        model.setupMatrix();

        for (double distance = 0.01; distance <= 1.005; distance += 0.01) {
            model.getTransitionProbabilities(distance, matrix);

            double pChange =
                    (matrix[1] + matrix[3]) * pi[0] +
                            (matrix[4] + matrix[6]) * pi[1] +
                            (matrix[9] + matrix[11]) * pi[2] +
                            (matrix[12] + matrix[14]) * pi[3];

            // analytical result for the probability of a mismatch in binary jukes cantor model
            double jc = 0.5 * (1 - Math.exp(-2.0 * distance));

            System.out.println(distance + "\t" + jc + "\t" + pChange);

            // this is an extremely weak test
            // I need to find some analytical results for this probability!
            assertTrue(pChange < jc);
        }
    }

    public void testSetupRelativeRates() throws Exception {

        model.setupMatrix();

        assertEquals(alpha.getParameterValue(0), model.relativeRates[0], 1e-8);
        assertEquals(switchingRate.getParameterValue(0), model.relativeRates[1], 1e-8);
        assertEquals(0.0, model.relativeRates[2], 1e-8);
        assertEquals(0.0, model.relativeRates[3], 1e-8);
        assertEquals(switchingRate.getParameterValue(0), model.relativeRates[4], 1e-8);
        assertEquals(1.0, model.relativeRates[5], 1e-8);
    }

    public void testNormalize() {

        model.setupMatrix();

        double[] pi = model.freqModel.getFrequencies();

        int stateCount = dataType.getStateCount();

        double totalRate = 0.0;
        for (int i = 0; i < stateCount; i++) {
            for (int j = 0; j < stateCount; j++) {
                int diff = Math.abs(i - j);
                if (diff != 2 && diff != 0) {
                    totalRate += model.q[i][j] * pi[i];
                }
            }
        }

        System.out.println(SubstitutionModelUtils.toString(model.q, dataType, 2));

        assertEquals(1.0, totalRate, 1e-8);
    }

    public static Test suite() {
        return new TestSuite(TwoStateCovarionModelTest.class);
    }

    TwoStateCovarionModel model;
    DataType dataType;
    Parameter frequencies;
    Parameter switchingRate;
    Parameter alpha;
}

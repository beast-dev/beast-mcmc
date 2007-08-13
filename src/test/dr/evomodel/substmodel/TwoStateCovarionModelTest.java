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

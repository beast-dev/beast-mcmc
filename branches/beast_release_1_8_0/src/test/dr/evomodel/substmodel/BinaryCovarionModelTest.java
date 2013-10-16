package test.dr.evomodel.substmodel;

import dr.evomodel.substmodel.*;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.TwoStateCovarion;
import dr.evolution.datatype.TwoStates;
import dr.inference.model.Parameter;
import dr.math.matrixAlgebra.IllegalDimension;
import dr.math.matrixAlgebra.Matrix;
import dr.math.matrixAlgebra.Vector;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * BinaryCovarionModel Tester.
 *
 * @author Alexei Drummond
 * @version 1.0
 * @since <pre>08/26/2007</pre>
 */
public class BinaryCovarionModelTest extends TestCase {
    public BinaryCovarionModelTest(String name) {
        super(name);
    }

    @Override
	public void setUp() throws Exception {
        super.setUp();

        frequencies = new Parameter.Default(new double[]{0.5, 0.5});
        hiddenFrequencies = new Parameter.Default(new double[]{0.5, 0.5});
        alpha = new Parameter.Default(0.0);
        switchingRate = new Parameter.Default(1.0);

        model = new BinaryCovarionModel(TwoStateCovarion.INSTANCE, frequencies, hiddenFrequencies, alpha, switchingRate);
        dataType = model.getDataType();
    }

    @Override
	public void tearDown() throws Exception {
        super.tearDown();
    }

    /**
     * Tests that pi*Q = 0
     */
    public void testEquilibriumDistribution() {

        alpha.setParameterValue(0, 0.1);
        switchingRate.setParameterValue(0, 1.0);

        model.setupMatrix();

        double[] pi = model.getFrequencyModel().getFrequencies();

        try {
            Matrix m = new Matrix(model.getQ());
            Vector p = new Vector(pi);
            Vector y = m.product(p);

            assertEquals(0.0, y.norm(), 1e-14);
        } catch (IllegalDimension illegalDimension) {
        }
    }

    public void testTransitionProbabilitiesAgainstEqualBaseFreqsEqualRates() {

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
    }

    public void testTransitionProbabilitiesUnequalBaseFreqsEqualRates() {

        // with alpha == 1, the transition probability should be the same as binary jukes cantor
        alpha.setParameterValue(0, 1.0);
        switchingRate.setParameterValue(0, 1.0);
        frequencies.setParameterValue(0, 0.25);
        frequencies.setParameterValue(1, 0.75);

        FrequencyModel freqModel = new FrequencyModel(TwoStates.INSTANCE, frequencies);

        GeneralSubstitutionModel modelToCompare = new GeneralSubstitutionModel(TwoStates.INSTANCE, freqModel, null, 0);

        model.setupMatrix();

        double[] matrix = new double[16];
        double[] m = new double[4];

        double[] pi = model.getFrequencyModel().getFrequencies();

        for (double distance = 0.01; distance <= 1.005; distance += 0.01) {
            model.getTransitionProbabilities(distance, matrix);

            modelToCompare.getTransitionProbabilities(distance, m);

            double pChange =
                    (matrix[1] + matrix[3]) * pi[0] +
                            (matrix[4] + matrix[6]) * pi[1] +
                            (matrix[9] + matrix[11]) * pi[2] +
                            (matrix[12] + matrix[14]) * pi[3];

            double pChange2 = m[1] * frequencies.getParameterValue(0) +
                    m[2] * frequencies.getParameterValue(1);

            System.out.println(distance + "\t" + pChange2 + "\t" + pChange);

            assertEquals(pChange2, pChange, 1e-14);
        }
    }

    public void testTransitionProbabilitiesUnequalBaseFreqsUnequalRateFreqsEqualRates() {

        // with alpha == 1, the transition probability should be the same as binary jukes cantor
        alpha.setParameterValue(0, 1.0);
        switchingRate.setParameterValue(0, 1.0);
        frequencies.setParameterValue(0, 0.25);
        frequencies.setParameterValue(1, 0.75);
        hiddenFrequencies.setParameterValue(0, 0.1);
        hiddenFrequencies.setParameterValue(1, 0.9);

        FrequencyModel freqModel = new FrequencyModel(TwoStates.INSTANCE, frequencies);

        GeneralSubstitutionModel modelToCompare = new GeneralSubstitutionModel(TwoStates.INSTANCE, freqModel, null, 0);

        model.setupMatrix();

        double[] matrix = new double[16];
        double[] m = new double[4];

        double[] pi = model.getFrequencyModel().getFrequencies();

        for (double distance = 0.01; distance <= 1.005; distance += 0.01) {
            model.getTransitionProbabilities(distance, matrix);

            modelToCompare.getTransitionProbabilities(distance, m);

            double pChange =
                    (matrix[1] + matrix[3]) * pi[0] +
                            (matrix[4] + matrix[6]) * pi[1] +
                            (matrix[9] + matrix[11]) * pi[2] +
                            (matrix[12] + matrix[14]) * pi[3];

            double pChange2 = m[1] * frequencies.getParameterValue(0) +
                    m[2] * frequencies.getParameterValue(1);

            //System.out.println(distance + "\t" + pChange2 + "\t" + pChange);

            assertEquals(pChange2, pChange, 1e-14);
        }
    }


    public void testTransitionProbabilitiesAgainstMatLab() {
        // test againt Matlab results for alpha = 0.5 and switching rate = 1.0
        // and visible state base frequencies = 0.25, 0.75

        alpha.setParameterValue(0, 0.5);
        switchingRate.setParameterValue(0, 1.0);
        frequencies.setParameterValue(0, 0.25);
        frequencies.setParameterValue(1, 0.75);

        model.setupMatrix();
        double[] matrix = new double[16];

        double[] pi = model.getFrequencyModel().getFrequencies();
        System.out.println(pi[0] + " " + pi[1] + " " + pi[2] + " " + pi[3]);

        System.out.println(SubstitutionModelUtils.toString(model.getQ(), dataType, 2));

        int index = 0;
        for (double distance = 0.01; distance <= 1.005; distance += 0.01) {
            model.getTransitionProbabilities(distance, matrix);


            double pChange =
                    (matrix[1] + matrix[3]) * pi[0] +
                            (matrix[4] + matrix[6]) * pi[1] +
                            (matrix[9] + matrix[11]) * pi[2] +
                            (matrix[12] + matrix[14]) * pi[3];


            double pChangeIndependent = TwoStateCovarionModelTest.matLabPChange[index];

            //System.out.println(distance + "\t" + pChange + "\t" + pChangeIndependent);
            assertEquals(pChange, pChangeIndependent, 1e-14);

            index += 1;
        }
    }

    public static Test suite() {
        return new TestSuite(BinaryCovarionModelTest.class);
    }

    BinaryCovarionModel model;
    DataType dataType;
    Parameter frequencies;
    Parameter hiddenFrequencies;
    Parameter switchingRate;
    Parameter alpha;
}

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
        alpha = new Parameter.Default(0.0);
        switchingRate = new Parameter.Default(1.0);

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

        // test againt Matlab results for alpha = 0.0 and switching rate = 1.0

        alpha.setParameterValue(0, 0.0);

        model.setupMatrix();

        int index = 0;
        for (double distance = 0.01; distance <= 1.005; distance += 0.01) {
            model.getTransitionProbabilities(distance, matrix);


            double pChange =
                    (matrix[1] + matrix[3]) * pi[0] +
                            (matrix[4] + matrix[6]) * pi[1] +
                            (matrix[9] + matrix[11]) * pi[2] +
                            (matrix[12] + matrix[14]) * pi[3];

            //System.out.println(distance + "\t" + pChange + "\t");

            double pChangeIndependent = matLabPChange[index];

            assertEquals(pChange, pChangeIndependent, 1e-12);

            index += 1;
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


    /*
      The following Matlab code was written by David Bryant to generate these probabilities:

      format long;
      Q = [-1 1; 1 -1]; %Rate matrix
      G = [-1 1; 1 -1]; %Transition matrix for rate classes
      D = diag([0 1]);
      I = eye(2,2);
      R = kron(D,Q) + kron(G,I);

      pi = diag([0.25 0.25 0.25 0.25]);
      C = [0 1 0 1; 1 0 1 0; 0 1 0 1; 1 0 1 0]

      rate = sum(sum(C.*(pi*R)))

      %Normalise
      R = R/rate;


      for i=1:1:100
          t = i*0.01;
          P = expm(R*t);

          pchange(i) = sum(sum(C.*(pi*P)));
      end

      pchange'

     */
    static final double[] matLabPChange = {
            0.00980393422979,
            0.01923096169941,
            0.02830281301764,
            0.03703982111886,
            0.04546101308968,
            0.05358419593960,
            0.06142603671519,
            0.06900213733121,
            0.07632710446737,
            0.08341461485617,
            0.09027747626621,
            0.09692768446489,
            0.10337647642583,
            0.10963438002890,
            0.11571126048426,
            0.12161636369662,
            0.12735835677156,
            0.13294536585258,
            0.13838501146502,
            0.14368444153122,
            0.14885036221082,
            0.15388906670942,
            0.15880646218995,
            0.16360809491163,
            0.16829917371373,
            0.17288459195312,
            0.17736894799786,
            0.18175656437174,
            0.18605150563920,
            0.19025759511336,
            0.19437843046503,
            0.19841739830503,
            0.20237768780772,
            0.20626230343884,
            0.21007407684676,
            0.21381567797243,
            0.21748962542940,
            0.22109829620214,
            0.22464393470764,
            0.22812866126218,
            0.23155447999260,
            0.23492328622860,
            0.23823687341039,
            0.24149693954353,
            0.24470509323091,
            0.24786285930967,
            0.25097168411915,
            0.25403294042419,
            0.25704793201638,
            0.26001789801471,
            0.26294401688518,
            0.26582741019804,
            0.26866914613991,
            0.27147024279689,
            0.27423167122378,
            0.27695435831350,
            0.27963918947981,
            0.28228701116573,
            0.28489863318896,
            0.28747483093525,
            0.29001634740951,
            0.29252389515418,
            0.29499815804348,
            0.29743979296179,
            0.29984943137373,
            0.30222768079306,
            0.30457512615706,
            0.30689233111257,
            0.30917983921961,
            0.31143817507783,
            0.31366784538091,
            0.31586933990373,
            0.31804313242664,
            0.32018968160088,
            0.32230943175922,
            0.32440281367517,
            0.32647024527434,
            0.32851213230089,
            0.33052886894220,
            0.33252083841431,
            0.33448841351085,
            0.33643195711777,
            0.33835182269605,
            0.34024835473464,
            0.34212188917540,
            0.34397275381195,
            0.34580126866418,
            0.34760774632983,
            0.34939249231488,
            0.35115580534388,
            0.35289797765168,
            0.35461929525773,
            0.35632003822405,
            0.35800048089794,
            0.35966089214044,
            0.36130153554147,
            0.36292266962238,
            0.36452454802695,
            0.36610741970139,
            0.36767152906409
    };

}

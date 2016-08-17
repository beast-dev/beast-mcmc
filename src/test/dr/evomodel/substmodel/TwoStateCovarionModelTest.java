package test.dr.evomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.evolution.datatype.TwoStateCovarion;
import dr.oldevomodel.substmodel.FrequencyModel;
import dr.oldevomodel.substmodel.SubstitutionModelUtils;
import dr.oldevomodel.substmodel.TwoStateCovarionModel;
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
        dataType = model.getDataType();
    }

    public void testTransitionProbabilities() {

        // with alpha == 1, the transition probability should be the same as binary jukes cantor
        alpha.setParameterValue(0, 1.0);
        switchingRate.setParameterValue(0, 1.0);

        model.setupMatrix();

        double[] matrix = new double[16];

        double[] pi = model.getFrequencyModel().getFrequencies();

        for (double distance = 0.01; distance <= 1; distance += 0.01) {
            model.getTransitionProbabilities(distance, matrix);

            double pChange =
                    (matrix[1] + matrix[3]) * pi[0] +
                            (matrix[4] + matrix[6]) * pi[1] +
                            (matrix[9] + matrix[11]) * pi[2] +
                            (matrix[12] + matrix[14]) * pi[3];

            // analytical result for the probability of a mismatch in binary jukes cantor model
            double jc = 0.5 * (1 - Math.exp(-2.0 * distance));

//            System.err.println("Testing d=" + distance);
            assertEquals(pChange, jc, 1e-14);
        }
    }
    /*
    public void testCompareToScilabCode() {

        // test against Scilab results for alpha = 0.0 and switching rate = 1.0, visible state freq = {0.25, 0.75}

        frequencies = new Parameter.Default(new double[]{0.125, 0.125, 0.375, 0.375});

        FrequencyModel freqModel = new FrequencyModel(TwoStateCovarion.INSTANCE, frequencies);
        model = new TwoStateCovarionModel(TwoStateCovarion.INSTANCE, freqModel, alpha, switchingRate);
        dataType = model.getDataType();

        alpha.setParameterValue(0, 0.5);
        switchingRate.setParameterValue(0, 1.0);

        model.setupMatrix();

        double[] matrix = new double[16];

        double[] pi = model.getFrequencyModel().getFrequencies();

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

            System.err.println("Testing against scilab d=" + distance);
            assertEquals(pChange, pChangeIndependent, 1e-14);

            index += 1;
        }
    } */

    public void testSetupRelativeRates() throws Exception {

        model.setupMatrix();

        assertEquals(alpha.getParameterValue(0), model.getRelativeRates()[0], 1e-8);
        assertEquals(switchingRate.getParameterValue(0), model.getRelativeRates()[1], 1e-8);
        assertEquals(0.0, model.getRelativeRates()[2], 1e-8);
        assertEquals(0.0, model.getRelativeRates()[3], 1e-8);
        assertEquals(switchingRate.getParameterValue(0), model.getRelativeRates()[4], 1e-8);
        assertEquals(1.0, model.getRelativeRates()[5], 1e-8);
    }

    public void testNormalize() {

        model.setupMatrix();

        double[] pi = model.getFrequencyModel().getFrequencies();

        int stateCount = dataType.getStateCount();

        double totalRate = 0.0;
        for (int i = 0; i < stateCount; i++) {
            for (int j = 0; j < stateCount; j++) {
                int diff = Math.abs(i - j);
                if (diff != 2 && diff != 0) {
                    totalRate += model.getQ()[i][j] * pi[i];
                }
            }
        }

        System.out.println(SubstitutionModelUtils.toString(model.getQ(), dataType, 2));

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
      // The following Scilab code was written by Alexei Drummond
      // (adapted from Matlab code by David Bryant)

      Q = [-1,1;1,-1];

      // frequencies of visible states
      pi = diag([0.25,0.75]);

      Q = Q*pi;

      G = [-1,1;1,-1];
      D = diag([0.5,1]);
      I = eye(2,2);
      R = kron(D,Q) + kron(G, I);

      // frequencies of hidden states
      f = diag([0.5, 0.5]);

      // frequencies of big matrix
      pif = kron(f, pi);

      C = [0,1,0,1;1,0,1,0;0,1,0,1;1,0,1,0]

      rate = sum(sum(C.*(pif*R)))

      Rn = R/rate;

      for i = 1:1:100
        t = i*0.01;
        P = expm(Rn*t);
        pchange(i) = sum(sum(C .*(pif*P)));
      end;


     */
    static final double[] matLabPChange = {
            0.009853754858257556,
            0.01942241145387763,
            0.028716586161758886,
            0.03774630552222203,
            0.046521051427055954,
            0.055049802312796506,
            0.06334107073053451,
            0.07140293762693993,
            0.07924308363985522,
            0.08686881768339709,
            0.09428710307179183,
            0.1015045814078466,
            0.10852759444084098,
            0.11536220407949097,
            0.1220142107282943,
            0.12848917009986108,
            0.134792408641597,
            0.14092903770220686,
            0.14690396655180843,
            0.15272191435884597,
            0.15838742121740867,
            0.16390485830985116,
            0.1692784372817475,
            0.17451221889905671,
            0.17961012105092147,
            0.18457592615563362,
            0.18941328802201546,
            0.19412573821362156,
            0.1987166919588106,
            0.20318945364578556,
            0.20754722193809755,
            0.21179309454286796,
            0.21593007266103278,
            0.2199610651462345,
            0.2238888923965756,
            0.22771629000123267,
            0.231445912161957,
            0.23508033490766292,
            0.23862205911867845,
            0.2420735133757424,
            0.24543705664747917,
            0.24871498082886684,
            0.25190951314210097,
            0.2550228184102454,
            0.25805700121315506,
            0.2610141079343201,
            0.26389612870652945,
            0.26670499926357183,
            0.2694426027045621,
            0.2721107711769275,
            0.27471128748357165,
            0.27724588661926514,
            0.27971625724089505,
            0.2821240430758165,
            0.2844708442722035,
            0.28675821869497176,
            0.28898768317056994,
            0.2911607146836553,
            0.2932787515284428,
            0.29534319441730017,
            0.2973554075489424,
            0.2993167196384335,
            0.30122842491099766,
            0.3030917840615263,
            0.3049080251815066,
            0.30667834465498195,
            0.3084039080250403,
            0.31008585083221957,
            0.3117252794261185,
            0.3133232717514247,
            0.314880878109485,
            0.3163991218964633,
            0.3178790003190806,
            0.31932148508884917,
            0.3207275230956765,
            0.3220980370616438,
            0.3234339261757309,
            0.3247360667102066,
            0.32600531261936083,
            0.32724249612123,
            0.3284484282629134,
            0.32962389947006454,
            0.33076968008109947,
            0.33188652086664117,
            0.3329751535346926,
            0.3340362912220076,
            0.3350706289721107,
            0.3360788442003818,
            0.337061597146628,
            0.33801953131551693,
            0.33895327390525487,
            0.3398634362248641,
            0.34075061410039437,
            0.3416153882704134,
            0.34245832477107263,
            0.34327997531106624,
            0.3440808776367703,
            0.34486155588783984,
            0.34562252094354134,
            0.3463642707600772
    };

}

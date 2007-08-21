package dr.evomodel.substmodel;

import dr.evolution.datatype.HiddenNucleotides;
import dr.inference.model.Parameter;
import junit.framework.Test;
import junit.framework.TestCase;
import junit.framework.TestSuite;

/**
 * CovarionHKY Tester.
 *
 * @author Alexei Drummond
 * @version 1.0
 * @since <pre>08/20/2007</pre>
 */
public class CovarionHKYTest extends TestCase {
    public CovarionHKYTest(String name) {
        super(name);
    }

    public void setUp() throws Exception {
        super.setUp();

        frequencies = new Parameter.Default(new double[]
                {0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125});
        alpha = new Parameter.Default(0.0);
        switchingRate = new Parameter.Default(1.0);
        kappa = new Parameter.Default(2.0);

        dataType = new HiddenNucleotides(2);

        FrequencyModel freqModel = new FrequencyModel(dataType, frequencies);
        model = new CovarionHKY(dataType, kappa, alpha, switchingRate, freqModel);
    }

    public void tearDown() throws Exception {
        super.tearDown();
    }

    public void testTransitionProbabilities() {

        alpha.setParameterValue(0, 0.0);
        switchingRate.setParameterValue(0, 1.0);
        kappa.setParameterValue(0, 2.0);

        model.setupMatrix();

        System.out.println(SubstitutionModelUtils.toString(model.q, dataType, 2));

        double[] matrix = new double[64];
        double[] pi = model.getFrequencyModel().getFrequencies();

        int index = 0;
        for (double distance = 0.01; distance <= 1.005; distance += 0.01) {
            model.getTransitionProbabilities(distance, matrix);

            double pChange = 0.0;
            double pNoOrHiddenChange = 0.0;
            int k = 0;
            for (int i = 0; i < 8; i++) {
                for (int j = 0; j < 8; j++) {

                    if ((i % 4) != (j % 4)) {
                        pChange += matrix[k] * pi[i];
                    } else {
                        pNoOrHiddenChange += matrix[k] * pi[i];
                    }
                    k += 1;
                }
            }

            double totalP = pChange + pNoOrHiddenChange;

            assertEquals(1.0, totalP, 1e-10);

            //System.out.print(distance + "\t" + "\t" + pChange + "\t");
            //if (index < 100) System.out.print(matLabPChange[index]);
            //System.out.println();

            assertEquals(matLabPChange[index], pChange, 1e-14);

            index += 1;
        }
    }

    public void testGetRelativeDNARates() throws Exception {
        //TODO: Test goes here...
    }

    public void testSetGetKappa() throws Exception {
        //TODO: Test goes here...
    }

    public static Test suite() {
        return new TestSuite(CovarionHKYTest.class);
    }

    CovarionHKY model;
    HiddenNucleotides dataType;
    Parameter frequencies;
    Parameter kappa;
    Parameter switchingRate;
    Parameter alpha;

    /*
     %The following Matlab code provided by David Bryant was used to calculate
     %the probabilities for a covarion K2P model with hidden rates of 0 and 1,
     %kappa of 2 and a switching rate of 1

     format long;
     clear all;

     r=4;
     k=2;

     Q = [-2-k 1 k 1; 1 -2-k 1 k; k 1 -2-k 1; 1 k 1 -2-k];

     piQ = diag([0.25 0.25 0.25 0.25]);

     %Next step necessary for unequal base cases. Commented out
     %to match Alexei's code.
     %Q = Q*piQ;

     G = [-1 1; 1 -1]; %Transition matrix for rate classes
     piG = diag([1/2 1/2]);
     D = diag([0 1]);
     I = eye(r,r);
     R = kron(D,Q) + kron(G,I);

     pi = kron(piG,piQ);

     %Matrix C is a 0 one matrix. Picks out where bases change
     C = kron(ones(2,2),ones(r,r)-eye(r));

     %Normalise
     rate = sum(sum(C.*(pi*R)));
     R = R/rate;


     for i=1:1:100
         t = i*0.01;
         P = expm(R*t);
         pchange(i) = sum(sum(C.*(pi*P)));
     end

     pchange'

    */
    static final double[] matLabPChange = {
            0.00986400785088,
            0.01946196037798,
            0.02880252532240,
            0.03789407770533,
            0.04674470983328,
            0.05536224095885,
            0.06375422660885,
            0.07192796759140,
            0.07989051869303,
            0.08764869707661,
            0.09520909039030,
            0.10257806459771,
            0.10976177153876,
            0.11676615623066,
            0.12359696391788,
            0.13025974687993,
            0.13675987100520,
            0.14310252213901,
            0.14929271221364,
            0.15533528516793,
            0.16123492266363,
            0.16699614960557,
            0.17262333947249,
            0.17812071946489,
            0.18349237547643,
            0.18874225689480,
            0.19387418123805,
            0.19889183863200,
            0.20379879613422,
            0.20859850190993,
            0.21329428926481,
            0.21788938053981,
            0.22238689087256,
            0.22678983183007,
            0.23110111491715,
            0.23532355496481,
            0.23945987340274,
            0.24351270142001,
            0.24748458301758,
            0.25137797795665,
            0.25519526460616,
            0.25893874269311,
            0.26261063595893,
            0.26621309472519,
            0.26974819837174,
            0.27321795773026,
            0.27662431739624,
            0.27996915796208,
            0.28325429817402,
            0.28648149701564,
            0.28965245572033,
            0.29276881971518,
            0.29583218049872,
            0.29884407745474,
            0.30180599960433,
            0.30471938729832,
            0.30758563385225,
            0.31040608712556,
            0.31318205104731,
            0.31591478708987,
            0.31860551569264,
            0.32125541763742,
            0.32386563537701,
            0.32643727431884,
            0.32897140406488,
            0.33146905960964,
            0.33393124249748,
            0.33635892194069,
            0.33875303589965,
            0.34111449212649,
            0.34344416917334,
            0.34574291736640,
            0.34801155974717,
            0.35025089298173,
            0.35246168823930,
            0.35464469204110,
            0.35680062708056,
            0.35893019301571,
            0.36103406723493,
            0.36311290559680,
            0.36516734314491,
            0.36719799479871,
            0.36920545602088,
            0.37119030346234,
            0.37315309558549,
            0.37509437326636,
            0.37701466037663,
            0.37891446434595,
            0.38079427670539,
            0.38265457361258,
            0.38449581635920,
            0.38631845186142,
            0.38812291313386,
            0.38990961974765,
            0.39167897827306,
            0.39343138270730,
            0.39516721488804,
            0.39688684489292,
            0.39859063142581,
            0.40027892219004
    };
}

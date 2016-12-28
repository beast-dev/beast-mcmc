package test.dr.evomodel.substmodel;

import dr.evolution.datatype.OldHiddenNucleotides;
import dr.evomodel.substmodel.*;
import dr.inference.model.Parameter;
import dr.oldevomodel.substmodel.CovarionHKY;
import dr.oldevomodel.substmodel.FrequencyModel;
import dr.oldevomodel.substmodel.SubstitutionModelUtils;
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

        alpha = new Parameter.Default(0.0);
        switchingRate = new Parameter.Default(1.0);
        kappa = new Parameter.Default(2.0);

        dataType = new OldHiddenNucleotides(2);

    }

    public void tearDown() throws Exception {
        super.tearDown();
    }

    public void testK2PTransitionProbabilities() {

        double[] frequencies = new double[]
                {0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125, 0.125};

        transitionProbabilitiesTester(frequencies, matLabPChange_K2P);
    }

    public void testHKYTransitionProbabilities() {

        double[] baseFrequencies = new double[]{0.15, 0.2, 0.3, 0.35};

        double[] frequencies = new double[8];
        for (int i = 0; i < frequencies.length; i++) {
            frequencies[i] = baseFrequencies[i % 4] / 2.0;
        }

        transitionProbabilitiesTester(frequencies, matLabPChange_HKY);
    }


    public void transitionProbabilitiesTester(double[] freqs, double[] expectedP) {

        Parameter frequencies = new Parameter.Default(freqs);
        FrequencyModel freqModel = new FrequencyModel(dataType, frequencies);
        model = new CovarionHKY(dataType, kappa, alpha, switchingRate, freqModel);

        alpha.setParameterValue(0, 0.0);
        switchingRate.setParameterValue(0, 1.0);
        kappa.setParameterValue(0, 2.0);

        model.setupMatrix();

        System.out.println(SubstitutionModelUtils.toString(model.getQ(), dataType, 2));

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

            System.out.print(distance + "\t" + "\t" + pChange + "\t");
            if (index < 100) System.out.print(expectedP[index]);
            System.out.println();

            //assertEquals(expectedP[index], pChange, 1e-14);

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
    OldHiddenNucleotides dataType;
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
    static final double[] matLabPChange_K2P = {
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

    static final double[] matLabPChange_HKY = {
            0.00985954695768,
            0.01944735733283,
            0.02877658973544,
            0.03785964647441,
            0.04670821854965,
            0.05533332794691,
            0.06374536739747,
            0.07195413775516,
            0.07996888313383,
            0.08779832394009,
            0.09545068792763,
            0.10293373939203,
            0.11025480661789,
            0.11742080768318,
            0.12443827471949,
            0.13131337672096,
            0.13805194098895,
            0.14465947329452,
            0.15114117683550,
            0.15750197006065,
            0.16374650342878,
            0.16987917516689,
            0.17590414608721,
            0.18182535351970,
            0.18764652441305,
            0.19337118765401,
            0.19900268565191,
            0.20454418523250,
            0.20999868788237,
            0.21536903938304,
            0.22065793887109,
            0.22586794735891,
            0.23100149574820,
            0.23606089236674,
            0.24104833005692,
            0.24596589284278,
            0.25081556220095,
            0.25559922295905,
            0.26031866884393,
            0.26497560770062,
            0.26957166640174,
            0.27410839546591,
            0.27858727340240,
            0.28300971079857,
            0.28737705416541,
            0.29169058955549,
            0.29595154596716,
            0.30016109854757,
            0.30432037160662,
            0.30843044145304,
            0.31249233906332,
            0.31650705259341,
            0.32047552974251,
            0.32439867997787,
            0.32827737662883,
            0.33211245885791,
            0.33590473351622,
            0.33965497689022,
            0.34336393634610,
            0.34703233187806,
            0.35066085756610,
            0.35425018294872,
            0.35780095431559,
            0.36131379592506,
            0.36478931115072,
            0.36822808356150,
            0.37163067793909,
            0.37499764123646,
            0.37832950348097,
            0.38162677862531,
            0.38488996534947,
            0.38811954781650,
            0.39131599638496,
            0.39447976828045,
            0.39761130822876,
            0.40071104905293,
            0.40377941223623,
            0.40681680845316,
            0.40982363807038,
            0.41280029161931,
            0.41574715024200,
            0.41866458611203,
            0.42155296283168,
            0.42441263580698,
            0.42724395260180,
            0.43004725327227,
            0.43282287068266,
            0.43557113080388,
            0.43829235299547,
            0.44098685027222,
            0.44365492955622,
            0.44629689191518,
            0.44891303278793,
            0.45150364219766,
            0.45406900495387,
            0.45660940084344,
            0.45912510481158,
            0.46161638713329,
            0.46408351357576,
            0.46652674555230
    };
}

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

            System.out.print(distance + "\t" + "\t" + pChange + "\t");
            if (index < 100) System.out.print(matLabPChange[index]);
            System.out.println();

            //assertEquals(matLabPChange[index], pChange, 1e-2);

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
    The following Matlab generated probabilities were provided by David Bryant
    */
    static final double[] matLabPChange = {
            0.00986467734329,
            0.01946717653895,
            0.02881967209327,
            0.03793366886586,
            0.04682004016820,
            0.05548906367886,
            0.06395045530081,
            0.07221340107909,
            0.08028658729021,
            0.08817822880831,
            0.09589609584712,
            0.10344753917083,
            0.11083951386199,
            0.11807860172905,
            0.12517103243190,
            0.13212270339876,
            0.13893919860391,
            0.14562580627166,
            0.15218753556809,
            0.15862913233866,
            0.16495509394649,
            0.17116968326280,
            0.17727694185824,
            0.18328070244081,
            0.18918460058374,
            0.19499208578385,
            0.20070643188890,
            0.20633074693013,
            0.21186798239388,
            0.21732094196467,
            0.22269228976989,
            0.22798455815467,
            0.23320015501386,
            0.23834137070647,
            0.24341038457651,
            0.24840927110271,
            0.25334000569846,
            0.25820447018181,
            0.26300445793469,
            0.26774167876885,
            0.27241776351551,
            0.27703426835435,
            0.28159267889695,
            0.28609441403850,
            0.29054082959121,
            0.29493322171183,
            0.29927283013496,
            0.30356084122345,
            0.30779839084614,
            0.31198656709290,
            0.31612641283628,
            0.32021892814848,
            0.32426507258187,
            0.32826576732097,
            0.33222189721310,
            0.33613431268466,
            0.34000383154964,
            0.34383124071634,
            0.34761729779830,
            0.35136273263474,
            0.35506824872573,
            0.35873452458701,
            0.36236221502887,
            0.36595195236364,
            0.36950434754558,
            0.37301999124733,
            0.37649945487625,
            0.37994329153421,
            0.38335203692407,
            0.38672621020576,
            0.39006631480492,
            0.39337283917683,
            0.39664625752796,
            0.39988703049786,
            0.40309560580335,
            0.40627241884734,
            0.40941789329415,
            0.41253244161336,
            0.41561646559387,
            0.41867035682989,
            0.42169449718048,
            0.42468925920406,
            0.42765500656937,
            0.43059209444424,
            0.43350086986327,
            0.43638167207586,
            0.43923483287539,
            0.44206067691092,
            0.44485952198215,
            0.44763167931877,
            0.45037745384488,
            0.45309714442957,
            0.45579104412416,
            0.45845944038711,
            0.46110261529708,
            0.46372084575493,
            0.46631440367526,
            0.46888355616801,
            0.47142856571076,
            0.47394969031224
    };
}

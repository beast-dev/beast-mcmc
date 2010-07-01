package test.dr.evomodel.substmodel;
import junit.framework.TestCase;
import dr.evolution.datatype.Microsatellite;
import dr.evomodel.substmodel.OnePhaseModel;
import dr.evomodel.substmodel.AsymmetricQuadraticModel;
import dr.evomodel.substmodel.LinearBiasModel;
import dr.inference.model.Parameter;

/**
 * @author Chieh-Hsi Wu
 *
 * Tests the LinearBiasModel of microsatellites.
 */
public class LinearBiasTest extends TestCase {
    interface Instance {

        public OnePhaseModel getSubModel();

        public double getBiasLinearParam();

        public double getBiasConstantParam();

        double getDistance();

        double[] getExpectedPi();

        double[] getExpectedResult();

        public boolean isLogistics();

    }

    Instance test0 = new Instance() {

        public OnePhaseModel getSubModel(){
            return new AsymmetricQuadraticModel(new Microsatellite(1,4),null,
                    new Parameter.Default(5.0), new Parameter.Default(3.0),new Parameter.Default(0.0),
                    new Parameter.Default(5.0), new Parameter.Default(3.0),new Parameter.Default(0.0), false);
        }
        public double getBiasConstantParam(){
            return 0.6;
        }


        public double getBiasLinearParam(){
            return -0.3;
        }


        public double getDistance() {
            return 0.1;
        }

        public double[] getExpectedPi() {
            return new double[]{
                    0.605108055009823,   0.324165029469548,   0.070726915520629,   0
            };
        }

        public double[] getExpectedResult() {
            return new double[]{
                    0.946658398354944,   0.052177678283089,   0.001163923361967,                   0,
                    0.097398332795100,   0.863963320210341,   0.038638346994560,                   0,
                    0.009958010985717,   0.177092423725065,   0.812949565289218,                   0,
                    0.000867275762735,   0.023191450073206,   0.212503945856679,   0.763437328307380
            };
        }

        public boolean isLogistics(){
            return false;
        }
    };

    Instance test1 = new Instance() {

        public OnePhaseModel getSubModel(){
            return new AsymmetricQuadraticModel(new Microsatellite(1,4), null,
                    new Parameter.Default(2.0), new Parameter.Default(5.0), new Parameter.Default(1.0),
                    new Parameter.Default(2.0), new Parameter.Default(5.0), new Parameter.Default(1.0), false);
        }
        public double getBiasConstantParam(){
            return 0.7;
        }


        public double getBiasLinearParam(){
            return -0.1;
        }


        public double getDistance() {
            return 0.76;
        }

        public double[] getExpectedPi() {
            return new double[]{
                    0.545073375262055,   0.238469601677149,   0.143081761006289,   0.073375262054507
            };
        }

        public double[] getExpectedResult() {
            return new double[]{
                    0.863562128338329,   0.108051762152360,   0.022697835971058,   0.005688273538253,
                    0.246975456348250,   0.483877169454682,   0.197140800250442,   0.072006573946625,
                    0.086467946556413,   0.328568000417404,   0.374198643716438,   0.210765409309745,
                    0.042255746284165,   0.234021365326532,   0.410992548154004,   0.312730340235300
            };
        }

        public boolean isLogistics(){
            return false;
        }

    };



    Instance test2 = new Instance() {

        public OnePhaseModel getSubModel(){
            return new AsymmetricQuadraticModel(new Microsatellite(1,6), null);
        }
        public double getBiasConstantParam(){
            return -0.1;
        }


        public double getBiasLinearParam(){
            return 0.2;
        }


        public double getDistance() {
            return 0.135;
        }

        public double[] getExpectedPi() {
            return new double[]{
                     0.0596248692859803, 0.0596248692859803, 0.0735548466843399, 0.111916502631505, 0.209948475478618, 0.485330436633581
            };
        }

        public double[] getExpectedResult() {
            return new double[]{
                9.06886087194020e-01, 8.80044177328101e-02, 4.90144978250249e-03, 0.000201113041761387, 6.73323616634035e-06, 1.99012740290826e-07,
                8.80044177328099e-02, 8.14017231936112e-01, 9.20028255435521e-02, 0.005709760680862604, 2.56214890258364e-04, 9.54921640436405e-06,
                3.97320252528957e-03, 7.45791296461646e-02, 8.13688411677705e-01, 0.100651260606737944, 6.76913614468138e-03, 3.38859399421337e-04,
                1.07145403446126e-04, 3.04194400509115e-03, 6.61509953263142e-02, 0.813368537652168544, 1.09045253167639e-01, 8.28612444534153e-03,
                1.91222310785769e-06, 7.27644213941150e-05, 2.37154744835561e-03, 0.058128373331922951, 8.13110664106114e-01, 1.26314738469105e-01,
                2.44495455773264e-08, 1.17316108146806e-06, 5.13562498673616e-05, 0.001910768413216536, 5.46423318430593e-02, 9.43394345883230e-01
            };
        }

        public boolean isLogistics(){
            return true;
        }

    };

    Instance[] all = {test0, test1, test2};
    public void testLinearBiasModel() {

        for (Instance test : all) {
            OnePhaseModel subModel = test.getSubModel();
            Microsatellite microsat = (Microsatellite)subModel.getDataType();
            Parameter biasLinear = new Parameter.Default(1, test.getBiasLinearParam());
            Parameter biasConstant = new Parameter.Default(1, test.getBiasConstantParam());
            LinearBiasModel lbm = new LinearBiasModel(
                    microsat,
                    null,
                    subModel,
                    biasConstant,
                    biasLinear,
                    test.isLogistics(),
                    false,
                    false);

            lbm.computeStationaryDistribution();

            double[] statDist = lbm.getStationaryDistribution();
            final double[] expectedStatDist = test.getExpectedPi();
            for (int k = 0; k < statDist.length; ++k) {
                assertEquals(statDist[k], expectedStatDist[k], 1e-10);
            }
            int stateCount = microsat.getStateCount();
            double[] mat = new double[stateCount*stateCount];
            lbm.getTransitionProbabilities(test.getDistance(), mat);
            final double[] result = test.getExpectedResult();

            int k;
            for (k = 0; k < mat.length; ++k) {
                assertEquals(result[k], mat[k], 5e-9);
                //System.out.print(" " + (mat[k] - result[k]));
            }

            k = 0;
            for(int i = 0; i < microsat.getStateCount(); i ++){
                for(int j = 0; j < microsat.getStateCount(); j ++){
                    assertEquals(result[k++], lbm.getOneTransitionProbabilityEntry(test.getDistance(), i , j), 5e-9);

                }
            }
            
            for(int j = 0; j < microsat.getStateCount();j ++){
                double[] colTransitionProb = lbm.getColTransitionProbabilities(test.getDistance(), j);
                for(int i =0 ; i < microsat.getStateCount(); i++){
                    assertEquals(result[i*microsat.getStateCount()+j], colTransitionProb[i], 5e-9);
                }
            }

            for(int i = 0; i < microsat.getStateCount();i ++){
                double[] rowTransitionProb = lbm.getRowTransitionProbabilities(test.getDistance(), i);
                for(int j =0 ; j < microsat.getStateCount(); j++){
                    assertEquals(result[i*microsat.getStateCount()+j], rowTransitionProb[j], 5e-9);
                }
            }
        }
    }

}

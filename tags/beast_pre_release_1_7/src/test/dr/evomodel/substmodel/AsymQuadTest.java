package test.dr.evomodel.substmodel;

import junit.framework.TestCase;
import dr.inference.model.Parameter;
import dr.evomodel.substmodel.AsymmetricQuadraticModel;
import dr.evolution.datatype.Microsatellite;

/**
 * @author Chieh-Hsi Wu
 *
 * Tests AsymmetricQuadraticModel exponentiation.
 */
public class AsymQuadTest extends TestCase {

    interface Instance {
        public Microsatellite getDataType();

        public double getExpanConst();

        public double getExpanLin();

        public double getExpanQuad();

        public double getContractConst();

        public double getContractLin();

        public double getContractQuad();

        double getDistance();

        double[] getExpectedPi();

        public double[] getExpectedResult();
    }

    Instance test0 = new Instance() {
        public Microsatellite getDataType(){
            return new Microsatellite(1,4);
        }
        public double getExpanConst(){
            return 1.0;
        }

        public double getExpanLin(){
            return 5.0;
        }

        public double getExpanQuad(){
            return 0.0;
        }

        public double getContractConst(){
            return 1.0;
        }

        public double getContractLin(){
            return 5.0;
        }

        public double getContractQuad(){
            return 0.0;
        }
                                                  
        public double getDistance() {
            return 0.1;
        }

        public double[] getExpectedPi() {
            return new double[]{
                    0.757532281205165, 0.126255380200861, 0.068866571018651, 0.047345767575323
            };
        }

        public double[] getExpectedResult() {
            return new double[]{
                    0.979555040783480,   0.019216583311851,   0.001139116232520,   0.000089259672149,
                    0.115299499871107,   0.780702902910835,   0.092806213742576,   0.011191383475483,
                    0.012530278557716,   0.170144725194722,   0.654730453041978,   0.162594543205584,
                    0.001428154754382,   0.029843689267955,   0.236501153753577,   0.732227002224086
            };
        }
    };

    Instance test1 = new Instance() {
        public Microsatellite getDataType(){
            return new Microsatellite(1,4);
        }
        public double getExpanConst(){
            return 1.0;
        }

        public double getExpanLin(){
            return 5.0;
        }

        public double getExpanQuad(){
            return 0.0;
        }

        public double getContractConst(){
            return 2.0;
        }

        public double getContractLin(){
            return 3.0;
        }

        public double getContractQuad(){
            return 0.0;
        }

        public double getDistance() {
            return 0.2;
        }

        public double[] getExpectedPi() {
            return new double[]{
                    0.666666666666667, 0.133333333333333, 0.100000000000000, 0.100000000000000
            };
        }

        public double[] getExpectedResult() {
            return new double[]{
                    0.965025560544615,   0.031394424214122,   0.003139726429413,   0.000440288811850,
                    0.156972121070610,   0.676199129838700,   0.136694646485644,   0.030134102605046,
                    0.020931509529421,   0.182259528647526,   0.546569017275914,   0.250239944547139,
                    0.002935258745666,   0.040178803473394,   0.250239944547139,   0.706645993233800
            };
        }


    };

    Instance test2 = new Instance() {
        public Microsatellite getDataType(){
            return new Microsatellite(1,4);
        }
        public double getExpanConst(){
            return 1.0;
        }

        public double getExpanLin(){
            return 5.0;
        }

        public double getExpanQuad(){
            return 3.0;
        }

        public double getContractConst(){
            return 2.0;
        }

        public double getContractLin(){
            return 3.0;
        }

        public double getContractQuad(){
            return 5.0;
        }

        public double getDistance() {
            return 0.3;
        }

        public double[] getExpectedPi() {
            return new double[]{
                    0.873099838521076, 0.087309983852108, 0.028063923381035, 0.011526254245782
            };
        }

        public double[] getExpectedResult() {
            return new double[]{
                    0.951679358560076,   0.039846739718488,   0.006618862899417,   0.001855038822019,
                    0.398467397184877,   0.419766204810971,   0.131559686528619,   0.050206711475533,
                    0.205920179092971,   0.409296802533481,   0.257041454224293,   0.127741564149256,
                    0.140516950382903,   0.380309775814669,   0.311022938798188,   0.168150335004240
            };
        }

    };

    Instance[] all = {test0, test1, test2};
    public void testAsymmetricQuadraticModel() {

        for (Instance test : all) {
            Parameter expanConst = new Parameter.Default(1,test.getExpanConst());
            Parameter expanLin = new Parameter.Default(1, test.getExpanLin());
            Parameter expanQuad = new Parameter.Default(1, test.getExpanQuad());
            Parameter contractConst = new Parameter.Default(1,test.getContractConst());
            Parameter contractLin = new Parameter.Default(1, test.getContractLin());
            Parameter contractQuad = new Parameter.Default(1, test.getContractQuad());

            Microsatellite microsat = test.getDataType();
            AsymmetricQuadraticModel aqm = new AsymmetricQuadraticModel(microsat,null,
                    expanConst, expanLin, expanQuad, contractConst, contractLin, contractQuad, false);

            aqm.computeStationaryDistribution();

            double[] statDist = aqm.getStationaryDistribution();
            final double[] expectedStatDist = test.getExpectedPi();
            for (int k = 0; k < statDist.length; ++k) {
                assertEquals(statDist[k], expectedStatDist[k], 1e-10);
            }

            double[] mat = new double[4*4];
            aqm.getTransitionProbabilities(test.getDistance(), mat);
            final double[] result = test.getExpectedResult();

            int k;
            for (k = 0; k < mat.length; ++k) {
                assertEquals(result[k], mat[k], 1e-10);
                // System.out.print(" " + (mat[k] - result[k]));
            }

            k = 0;
            for(int i = 0; i < microsat.getStateCount(); i ++){
                for(int j = 0; j < microsat.getStateCount(); j ++){
                    assertEquals(result[k++], aqm.getOneTransitionProbabilityEntry(test.getDistance(), i , j), 1e-10);

                }
            }

            for(int j = 0; j < microsat.getStateCount();j ++){
                double[] colTransitionProb = aqm.getColTransitionProbabilities(test.getDistance(), j);
                for(int i =0 ; i < microsat.getStateCount(); i++){
                    assertEquals(result[i*microsat.getStateCount()+j], colTransitionProb[i], 1e-10);
                }
            }

            for(int i = 0; i < microsat.getStateCount();i ++){
                double[] rowTransitionProb = aqm.getRowTransitionProbabilities(test.getDistance(), i);
                for(int j =0 ; j < microsat.getStateCount(); j++){
                    assertEquals(result[i*microsat.getStateCount()+j], rowTransitionProb[j], 1e-10);
                }
            }
        }
    }
}

package test.dr.evomodel.substmodel;

import junit.framework.TestCase;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.NtdBMA;
import dr.evolution.datatype.Nucleotides;

/**
 * @author Chieh-Hsi Wu
 *
 * JUnit test for NtdBMA model
 */
public class testNtdBMA extends TestCase {
    interface Instance {
        double[] getPi();

        double getLogKappa();

        double getLogTN();

        double getLogAC();

        double getLogAT();
        double getLogGC();
        double getLogGT();
        Variable<Integer> getModelChoose();

        double getDistance();

        double[] getExpectedResult();
    }

    /*
     * Results obtained by running the following scilab code,
     *
     * k = 5 ; piQ = diag([.2, .3, .25, .25]) ; d = 0.1 ;
     * % Q matrix with zeroed diagonal
     * XQ = [0 1 k 1; 1 0 1 k; k 1 0 1; 1 k 1 0];
     *
     * xx = XQ * piQ ;
     *
     * % fill diagonal and normalize by total substitution rate
     * q0 = (xx + diag(-sum(xx,2))) / sum(piQ * sum(xx,2)) ;
     * expm(q0 * d)
     */

    //A HKY model
    Instance test0 = new Instance() {
        public double[] getPi() {
            return new double[]{0.25, 0.25, 0.25, 0.25};
        }

        public double getLogKappa() {
            return Math.log(2);
        }

        public double getLogTN(){
            return Math.log(1.2);
        }

        public double getLogAC(){
            return Math.log(0.5);
        }

        public double getLogAT(){
            return Math.log(0.5);
        }

        public double getLogGC(){
            return Math.log(0.5);
        }
        public double getLogGT(){
            return Math.log(0.5);
        }

        public Variable<Integer> getModelChoose(){
            return new Variable.I(new int[]{0, 0});
        }

        public double getDistance() {
            return 0.1;
        }

        public double[] getExpectedResult() {
            return new double[]{
                    0.906563342722, 0.023790645491, 0.045855366296, 0.023790645491,
                    0.023790645491, 0.906563342722, 0.023790645491, 0.045855366296,
                    0.045855366296, 0.023790645491, 0.906563342722, 0.023790645491,
                    0.023790645491, 0.045855366296, 0.023790645491, 0.906563342722
            };
        }
    };

    //A TN93 model
    Instance test1 = new Instance() {
        public double[] getPi() {
            return new double[]{0.1, 0.2, 0.3, 0.4};
        }

        public double getLogKappa() {
            return Math.log(3)-0.5*Math.log(1.5);
        }

        public double getLogTN(){
            return -0.5*Math.log(1.5);
        }

        public double getLogAC(){
            return Math.log(0.5);
        }

        public double getLogAT(){
            return Math.log(0.5);
        }

        public double getLogGC(){
            return Math.log(0.5);
        }
        public double getLogGT(){
            return Math.log(0.5);
        }

        public Variable<Integer> getModelChoose(){
            return new Variable.I(new int[]{1, 0});
        }


        public double getDistance() {
            return 0.1;
        }

        public double[] getExpectedResult() {
            return new double[]{
                    0.895550254199242, 0.017687039418335, 0.051388627545752, 0.035374078836670,
                    0.008843519709168, 0.865344657365451, 0.026530559127503, 0.099281263797879,
                    0.017129542515251, 0.017687039418335, 0.929809339229744, 0.035374078836670,
                    0.008843519709168, 0.049640631898940, 0.026530559127503, 0.914985289264390
            };
        }
    };

    //GTR example
    Instance test2 = new Instance() {
        public double[] getPi() {
            return new double[]{0.20, 0.30, 0.25, 0.25};
        }

        public double getLogKappa() {
            return Math.log(3)-0.5*Math.log(1.5);
        }

        public double getLogTN(){
            return -0.5*Math.log(1.5);
        }

        public double getLogAC(){
            return Math.log(1.2);
        }

        public double getLogAT(){
            return Math.log(0.6);
        }

        public double getLogGC(){
            return Math.log(0.5);
        }
        public double getLogGT(){
            return Math.log(0.8);
        }

        public Variable<Integer> getModelChoose(){
            return new Variable.I(new int[]{1, 1});
        }

        public double getDistance() {
            return 0.1;
        }

        public double[] getExpectedResult() {
            return new double[]{
                    0.9078362845301, 0.0325116185198, 0.0449673267333, 0.0146847702168,
                    0.0216744123465, 0.9006273487178, 0.0122790622489, 0.0654191766868,
                    0.0359738613867, 0.0147348746987, 0.9308468616493, 0.0184444022654,
                    0.0117478161734, 0.0785030120241, 0.0184444022654, 0.8913047695370
            };
        }
    };

    Instance test3 = new Instance() {
        public double[] getPi() {
            return new double[]{0.25, 0.25, 0.25, 0.25};
        }

        public double getLogKappa() {
            return Math.log(2);
        }

        public double getLogTN(){
            return Math.log(1.2);
        }

        public double getLogAC(){
            return Math.log(0.5);
        }

        public double getLogAT(){
            return Math.log(0.5);
        }

        public double getLogGC(){
            return Math.log(0.5);
        }
        public double getLogGT(){
            return Math.log(0.5);
        }

        public Variable<Integer> getModelChoose(){
            return new Variable.I(new int[]{0, 0});
        }

        public double getDistance() {
            return 1.8;
        }

        public double[] getExpectedResult() {
            return new double[]{
                    0.324927478425, 0.208675277945, 0.257721965686, 0.208675277945,
                    0.208675277945, 0.324927478425, 0.208675277945, 0.257721965686,
                    0.257721965686, 0.208675277945, 0.324927478425, 0.208675277945,
                    0.208675277945, 0.257721965686, 0.208675277945, 0.324927478425
            };
        }
    };

    Instance test4 = new Instance() {
        public double[] getPi() {
            return new double[]{0.1, 0.2, 0.3, 0.4};
        }

        public double getLogKappa() {
            return Math.log(3)-0.5*Math.log(1.5);
        }

        public double getLogTN(){
            return -0.5*Math.log(1.5);
        }

        public double getLogAC(){
            return Math.log(0.5);
        }

        public double getLogAT(){
            return Math.log(0.5);
        }

        public double getLogGC(){
            return Math.log(0.5);
        }
        public double getLogGT(){
            return Math.log(0.5);
        }

        public Variable<Integer> getModelChoose(){
            return new Variable.I(new int[]{1, 0});
        }


        public double getDistance() {
            return 2.5;
        }

        public double[] getExpectedResult() {
            return new double[]{
                    0.144168843021, 0.180243104854, 0.315101842417, 0.360486209708,
                    0.090121552427, 0.217265980316, 0.270364657281, 0.422247809976,
                    0.105033947472, 0.180243104854, 0.354236737965, 0.360486209708,
                    0.090121552427, 0.211123904988, 0.270364657281, 0.428389885304
            };
        }
    };

    //GTR example
    Instance test5 = new Instance() {
        public double[] getPi() {
            return new double[]{0.20, 0.30, 0.25, 0.25};
        }

        public double getLogKappa() {
            return Math.log(3)-0.5*Math.log(1.5);
        }

        public double getLogTN(){
            return -0.5*Math.log(1.5);
        }

        public double getLogAC(){
            return Math.log(1.2);
        }

        public double getLogAT(){
            return Math.log(0.6);
        }

        public double getLogGC(){
            return Math.log(0.5);
        }
        public double getLogGT(){
            return Math.log(0.8);
        }

        public Variable<Integer> getModelChoose(){
            return new Variable.I(new int[]{1, 1});
        }

        public double getDistance() {
            return 2.5;
        }

        public double[] getExpectedResult() {
            return new double[]{
                    0.246055801088, 0.266163561908, 0.273492078437, 0.214288558567,
                    0.177442374606, 0.341245862774, 0.202456669039, 0.278855093581,
                    0.218793662750, 0.242948002847, 0.333259562500, 0.204998771904,
                    0.171430846853, 0.334626112297, 0.204998771904, 0.288944268946
            };
        }
    };

    Instance[] all = {test0,test1,test2,test3,test4,test5};

    public void testNtdBMA() {
        for (Instance test : all) {
            Parameter logKappa = new Parameter.Default(1, test.getLogKappa());
            Parameter logTN = new Parameter.Default(1, test.getLogTN());
            Parameter logAC = new Parameter.Default(1, test.getLogAC());
            Parameter logAT = new Parameter.Default(1, test.getLogAT());
            Parameter logGC = new Parameter.Default(1, test.getLogGC());
            Parameter logGT = new Parameter.Default(1, test.getLogGT());
            Variable<Integer> modelChoose = test.getModelChoose();

            double[] pi = test.getPi();

            Parameter freqs = new Parameter.Default(pi);
            FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
            NtdBMA ntdBMA = new NtdBMA(
                    logKappa,
                    logTN,
                    logAC,
                    logAT,
                    logGC,
                    logGT,
                    modelChoose,
                    f
            );

            double distance = test.getDistance();

            double[] mat = new double[4 * 4];
            ntdBMA.getTransitionProbabilities(distance, mat);
            final double[] result = test.getExpectedResult();

            for (int k = 0; k < mat.length; ++k) {
                assertEquals(mat[k], result[k], 5e-10);
                // System.out.print(" " + (mat[k] - result[k]));
            }
        }
    }
}

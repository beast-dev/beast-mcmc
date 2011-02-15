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

    Instance test1 = new Instance() {
        public double[] getPi() {
            return new double[]{0.50, 0.20, 0.2, 0.1};
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
                    0.928287993055, 0.021032136637, 0.040163801989, 0.010516068319,
                    0.052580341593, 0.906092679369, 0.021032136637, 0.020294842401,
                    0.100409504972, 0.021032136637, 0.868042290072, 0.010516068319,
                    0.052580341593, 0.040589684802, 0.021032136637, 0.885797836968
            };
        }
    };

    Instance test2 = new Instance() {
        public double[] getPi() {
            return new double[]{0.20, 0.30, 0.25, 0.25};
        }

        public double getLogKappa() {
            return Math.log(5);
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
                    0.904026219693, 0.016708646875, 0.065341261036, 0.013923872396,
                    0.011139097917, 0.910170587813, 0.013923872396, 0.064766441875,
                    0.052273008829, 0.016708646875, 0.917094471901, 0.013923872396,
                    0.011139097917, 0.077719730250, 0.013923872396, 0.897217299437
            };
        }
    };

    Instance[] all = {test2, test1, test0};

    public void testHKY() {
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
                assertEquals(mat[k], result[k], 1e-10);
                // System.out.print(" " + (mat[k] - result[k]));
            }
        }
    }
}

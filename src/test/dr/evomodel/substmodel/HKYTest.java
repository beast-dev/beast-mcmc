package test.dr.evomodel.substmodel;

import dr.evolution.datatype.Nucleotides;
import dr.oldevomodel.substmodel.FrequencyModel;
import dr.oldevomodel.substmodel.HKY;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

/**
 * Test HKY matrix exponentiation
 *
 * @author Joseph Heled
 *         Date: 7/11/2007
 */
public class HKYTest extends TestCase {

    interface Instance {
        double[] getPi();

        double getKappa();

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

        public double getKappa() {
            return 2;
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

        public double getKappa() {
            return 2;
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

        public double getKappa() {
            return 5;
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
            Parameter kappa = new Parameter.Default(1, test.getKappa());
            double[] pi = test.getPi();

            Parameter freqs = new Parameter.Default(pi);
            FrequencyModel f = new FrequencyModel(Nucleotides.INSTANCE, freqs);
            HKY hky = new HKY(kappa, f);

            double distance = test.getDistance();

            double[] mat = new double[4 * 4];
            hky.getTransitionProbabilities(distance, mat);
            final double[] result = test.getExpectedResult();

            for (int k = 0; k < mat.length; ++k) {
                assertEquals(mat[k], result[k], 1e-10);
                // System.out.print(" " + (mat[k] - result[k]));
            }
        }
    }
}

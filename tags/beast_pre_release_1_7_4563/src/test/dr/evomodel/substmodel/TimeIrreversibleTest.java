package test.dr.evomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.evomodel.substmodel.ComplexSubstitutionModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SVSComplexSubstitutionModel;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Random;

/**
 *
 */
public class TimeIrreversibleTest extends TestCase {
    private static final double time = 0.01;

    private static NumberFormat formatter = new DecimalFormat("###0.000000");

    private static ArrayList<Double> ratioSummary = new ArrayList<Double> ();

    static class Original {
        public double[] getRates() {
            return new double[]{
                    38.505, 23.573, 2.708, 3.35, 11.641, 0.189, 0.127, 0.511, 0.272, 1.788E-2,
                    0.214, 1.322E-2, 3.015E-2, 0.449, 0.177, 0.305, 1.517E-2, 3.924E-2, 0.18, 0.14,
                    1.273E-2, 0.265, 1.422E-2, 1.474E-2, 0.911, 0.17, 0.217, 4.078, 0.206, 3.309E-2,
                    0.657, 1.874E-2, 3.141E-2, 0.403, 2.003E-2, 0.582, 0.732, 0.106, 9.147E-2, 0.248,
                    1.516E-2, 0.524, 0.1, 1.986, 0.819, 0.146, 7.519E-2, 1.35, 0.166, 0.204,
                    1.753E-2, 0.59, 0.691, 3.308, 0.377, 1.785E-2
            }; // 56
        }

        public double[] getIndicators() {
            return new double[]{
                    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                    1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
                    1.0, 1.0, 1.0, 1.0, 1.0, 1.0
            }; // 56
        }

        public double[] getFrequencies() {
//            return new double[]{0.141, 7.525E-2, 0.117, 0.283, 7.743E-2, 0.156, 1.555E-2, 0.135};  // sum to 1.00023
            return new double[]{0.141, 0.075, 0.117, 0.283, 0.077, 0.156, 0.016, 0.135};
        }

        public DataType getDataType() {
            return new DataType() {
                public String getDescription() {
                    return null;
                }

                public int getType() {
                    return 0;   // NUCLEOTIDES = 0;
                }

                @Override
                public char[] getValidChars() {
                    return null;
                }

                public int getStateCount() {
                    return getFrequencies().length;  // = frequency
                }
            };
        }

        public String toString() {
            return "original data test";
        }
    }

    class Test extends Original {
        private final double x;

        public Test(double x) {
            this.x = x;
        }

        public double[] getRates(int id) {
            double[] originalRates = super.getRates();
            System.out.println("original rates:");
            printRateMatrix(originalRates, getDataType().getStateCount());
            double[] newRates = new double[originalRates.length];

            double[] uniform = new double[originalRates.length];

            for (int r = 0; r < originalRates.length; r++) {
                if (r == id) {
                    uniform[r] = (new Random()).nextDouble() * ((1 / x) - x) + x;
                    newRates[r] = originalRates[r] * uniform[r];
                } else {
                    newRates[r] = originalRates[r];
                }
            }
            System.out.println("random ratio:");
            printRateMatrix(uniform, getDataType().getStateCount());
            System.out.println("new rates:");
            printRateMatrix(newRates, getDataType().getStateCount());
            return newRates;
        }

        public String toString() {
            return "test using random number : " + x;
        }
    }

    public void tests() {
        Original originalTest = new Original();
        double[] csm_orig = testComplexSubstitutionModel(originalTest, originalTest.getRates());
        double[] svs_orig = testSVSComplexSubstitutionModel(originalTest, originalTest.getRates());

        Test test = new Test(0.8);
        for (int r = 0; r < test.getRates().length; r++) {
            System.out.println("==================== changing index = " + r + " (start from 0) ====================");

            double[] newRate = test.getRates(r);
            double[] csm_test = testComplexSubstitutionModel(test, newRate);
            reportMatrix(csm_orig, csm_test);

            double[] svs_test = testSVSComplexSubstitutionModel(test, newRate);
            reportMatrix(svs_orig, svs_test);
        }

        System.out.println("==================== Biggest Ratio Summary ====================\n");
        int i = 1;
        double bigget = 0;
        int biggetId = 0;

        for (Double r : ratioSummary) {
            if (i % 2 != 0) {
                System.out.print(i/2 + "   ");
            }
            System.out.print(formatter.format(r) + ",  ");
            if (bigget < r) {
                bigget = r;
                biggetId = i;
            }
            if (i % 2 == 0) {
                System.out.println("");
            }

            i++;
        }

        System.out.println("bigget = " + formatter.format(bigget) + ", where index is " + biggetId/2);
    }

    private double[] testComplexSubstitutionModel(Original test, double[] rates) {
        System.out.println("\n*** Complex Substitution Model Test: " + test + " ***");

        Parameter ratesP = new Parameter.Default(rates);

        DataType dataType = test.getDataType();

        FrequencyModel freqModel = new FrequencyModel(dataType, new Parameter.Default(test.getFrequencies()));

        ComplexSubstitutionModel substModel = new ComplexSubstitutionModel("Complex Substitution Model Test",
                dataType,
                freqModel,
                ratesP);
        double logL = substModel.getLogLikelihood();

        System.out.println("Prior = " + logL);

        double[] finiteTimeProbs = null;
        if (!Double.isInfinite(logL)) {
            finiteTimeProbs = new double[substModel.getDataType().getStateCount() * substModel.getDataType().getStateCount()];
            substModel.getTransitionProbabilities(time, finiteTimeProbs);
            System.out.println("Probs = ");
            printRateMatrix(finiteTimeProbs, substModel.getDataType().getStateCount());
        }
//            assertEquals(1, 1, 1e-10);
        return finiteTimeProbs;
    }

    private double[] testSVSComplexSubstitutionModel(Original test, double[] rates) {
        System.out.println("\n*** SVS Complex Substitution Model Test: " + test + " ***");

        double[] indicators = test.getIndicators();

        Parameter ratesP = new Parameter.Default(rates);
        Parameter indicatorsP = new Parameter.Default(indicators);

        DataType dataType = test.getDataType();

        FrequencyModel freqModel = new FrequencyModel(dataType, new Parameter.Default(test.getFrequencies()));

        SVSComplexSubstitutionModel substModel = new SVSComplexSubstitutionModel("SVS Complex Substitution Model Test",
                dataType,
                freqModel,
                ratesP, indicatorsP);

        double logL = substModel.getLogLikelihood();

        System.out.println("Prior = " + logL);

        double[] finiteTimeProbs = null;
        if (!Double.isInfinite(logL)) {
            finiteTimeProbs = new double[substModel.getDataType().getStateCount() * substModel.getDataType().getStateCount()];
            substModel.getTransitionProbabilities(time, finiteTimeProbs);
            System.out.println("Probs = ");
            printRateMatrix(finiteTimeProbs, substModel.getDataType().getStateCount());
        }
//            assertEquals(1, 1, 1e-10);
        return finiteTimeProbs;
    }

    public static void printRateMatrix(double[] m, int a) {
        int id = 0;
        for (int i = 0; i < a; i++) {
            if (i == 0) {
                System.out.print("/  ");
            } else if (i == a - 1) {
                System.out.print("\\  ");
            } else {
                System.out.print("|  ");
            }

            for (int j = 0; j < a; j++) {
                if (i == j) {
                    System.out.print("null");
                    for (int n = 0; n < formatter.getMaximumFractionDigits(); n++) {
                        System.out.print(" ");
                    }
                } else {
                    System.out.print(formatter.format(m[id]) + "  ");
                    id++;                                
                }
            }

            if (i == 0) {
                System.out.print("\\");
            } else if (i == a - 1) {
                System.out.print("/");
            } else {
                System.out.print("|");
            }

            System.out.println();
        }
        System.out.println("\n");
    }

    public static void reportMatrix(double[] orig, double[] test) {
        double bigRatio = 0;
        double ratio;
        int index = -1;

        if (orig.length != test.length)
            System.err.println("Error : 2 matrix should have same length ! " + orig.length + " " + test.length);

        for (int i = 0; i < orig.length; i++) {
            ratio = Math.abs(orig[i] / test[i]);
            if (bigRatio < ratio) {
                bigRatio = ratio;
                index = i;
            }
        }

        ratioSummary.add(bigRatio);

        System.out.println("Biggest Ratio = " + formatter.format(bigRatio) + ", between " + formatter.format(orig[index])
                + " and " + formatter.format(test[index]));
        System.out.println("index = " + index + " (start from 0)");
        System.out.println("\n");
    }
}

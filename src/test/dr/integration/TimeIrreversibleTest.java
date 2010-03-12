package test.dr.integration;

import dr.evolution.datatype.DataType;
import dr.evomodel.substmodel.ComplexSubstitutionModel;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.SVSComplexSubstitutionModel;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.util.Random;

/**
 *
 */
public class TimeIrreversibleTest extends TestCase {
    private static NumberFormat formatter = new DecimalFormat("###0.000");

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

        public Test (double x) {
            this.x = x;
        }

        public double[] getRates() {
            double[] originalRates = super.getRates();
            System.out.println("original rates:"); 
            printArrayMatrix(originalRates, getDataType().getStateCount(), getDataType().getStateCount());
            double[] newRates = new double[originalRates.length];

            double[] uniform = new double[originalRates.length];

            for (int r = 0; r < originalRates.length; r++) {
                uniform[r] = (new Random()).nextDouble() * ((1/x) - x) + x;
                newRates[r] = originalRates[r] * uniform[r];
            }
            System.out.println("random ratio:");
            printArrayMatrix(uniform, getDataType().getStateCount(), getDataType().getStateCount());
            System.out.println("new rates:");
            printArrayMatrix(newRates, getDataType().getStateCount(), getDataType().getStateCount());
            return newRates;
        }

        public String toString() {
            return "test using random number : " + x;
        }
    }

    public void tests() {
        Original originalTest = new Original();
        testComplexSubstitutionModel(originalTest);
        testSVSComplexSubstitutionModel(originalTest);

        Test test = new Test(0.8);
        testComplexSubstitutionModel(test);
        testSVSComplexSubstitutionModel(test);
    }

    private void testComplexSubstitutionModel(Original test) {
        System.out.println("Complex Substitution Model Test: " + test);
                
        double[] rates = test.getRates();
        Parameter ratesP = new Parameter.Default(rates);

        DataType dataType = test.getDataType();

        FrequencyModel freqModel = new FrequencyModel(dataType, new Parameter.Default(test.getFrequencies()));

        ComplexSubstitutionModel substModel = new ComplexSubstitutionModel("Complex Substitution Model Test",
                dataType,
                freqModel,
                ratesP);
        double logL = substModel.getLogLikelihood();

        System.out.println("Prior = " + logL);

        if (!Double.isInfinite(logL)) {
            double[] finiteTimeProbs = new double[substModel.getDataType().getStateCount() * substModel.getDataType().getStateCount()];
            double time = 1.0;
            substModel.getTransitionProbabilities(time, finiteTimeProbs);
            System.out.println("Probs = ");
            printArrayMatrix(finiteTimeProbs, substModel.getDataType().getStateCount(), substModel.getDataType().getStateCount());
        }
//            assertEquals(1, 1, 1e-10);

    }

    private void testSVSComplexSubstitutionModel(Original test) {
        System.out.println("SVS Complex Substitution Model Test: " + test);
                
        double[] rates = test.getRates();
        double[] indicators = test.getIndicators();

        Parameter ratesP = new Parameter.Default(rates);
        Parameter indicatorsP = new Parameter.Default(indicators);

        DataType dataType = test.getDataType();

        FrequencyModel freqModel = new FrequencyModel(dataType, new Parameter.Default(test.getFrequencies()));

        SVSComplexSubstitutionModel substModel = new SVSComplexSubstitutionModel("SVS Complex Substitution Model Test",
                dataType,
                freqModel,
                ratesP, indicatorsP, true);

        double logL = substModel.getLogLikelihood();

        System.out.println("Prior = " + logL);

        if (!Double.isInfinite(logL)) {
            double[] finiteTimeProbs = new double[substModel.getDataType().getStateCount() * substModel.getDataType().getStateCount()];
            double time = 1.0;
            substModel.getTransitionProbabilities(time, finiteTimeProbs);
            System.out.println("Probs = ");
            printArrayMatrix(finiteTimeProbs, substModel.getDataType().getStateCount(), substModel.getDataType().getStateCount());
        }
//            assertEquals(1, 1, 1e-10);
    }

    public static void printArrayMatrix(double[] m, int a, int b) {
        for (int i = 0; i < a; i++) {
           if (i == 0) {
               System.out.print("/  ");
           } else if (i == a - 1) {
               System.out.print("\\  ");
           } else {
               System.out.print("|  ");
           }

           for (int j = 0; j < b; j++) {
               System.out.print(formatter.format(m[i+j]) + "  ");
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
}

package test.dr.evomodel.substmodel;

import junit.framework.TestCase;
import dr.app.beagle.evomodel.substmodel.*;
import dr.evolution.datatype.TwoStates;
import dr.math.matrixAlgebra.Vector;
import dr.inference.model.Parameter;

import java.util.List;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * @author Marc A. Suchard
 */

public class ProductChainSubstitutionModelTest extends TestCase {

    private void setUpTwoStates() {
        FrequencyModel freqModel0 = new FrequencyModel(TwoStates.INSTANCE,
                new double[]{1.0 / 3.0, 2.0 / 3.0});
        FrequencyModel freqModel1 = new FrequencyModel(TwoStates.INSTANCE,
                new double[]{1.0 / 4.0, 3.0 / 4.0});

        GeneralSubstitutionModel substModel0 = new GeneralSubstitutionModel("model0",
                TwoStates.INSTANCE, freqModel0, new Parameter.Default(new double[]{1}), 1);
        GeneralSubstitutionModel substModel1 = new GeneralSubstitutionModel("model1",
                TwoStates.INSTANCE, freqModel1, new Parameter.Default(new double[]{1}), 1);

        baseModels = new ArrayList<SubstitutionModel>();
        baseModels.add(substModel0);
        baseModels.add(substModel1);

        productChainModel = new ProductChainSubstitutionModel("productChain", baseModels);
        stateCount = 2 * 2;

        /*
            model0 = as.eigen.two.state(1.5, 0.75)
            model1 = as.eigen.two.state(2.0, 2.0 / 3.0)
            pc = ind.two.eigen(model0, model1)
            pc$rate.matrix
            matexp(pc, 0.5)
        */

        markovJumpsInfinitesimalResult = new double[]{
                -3.5000000, 2.000000, 1.5000000, 0.000000,
                0.6666667, -2.166667, 0.0000000, 1.500000,
                0.7500000, 0.000000, -2.7500000, 2.000000,
                0.0000000, 0.750000, 0.6666667, -1.416667
        };
        markovJumpsEigenValues = new double[]{
                -4.916667, -2.666667,  -2.250000, 0.000000
        };

        markovJumpsProbs = new double[]{
                0.24613009, 0.3036382, 0.20156776, 0.2486639,
                0.10121274, 0.4485556, 0.08288798, 0.3673437,
                0.10078388, 0.1243320, 0.34691397, 0.4279702,
                0.04144399, 0.1836719, 0.14265673, 0.6322274
        };
    }

    private void assertEquals(double[] a, double[] b, double accuracy) {
        assertEquals(a.length, b.length);
        for (int i = 0; i < a.length; i++) {
            assertEquals(a[i], b[i], accuracy);
        }
    }

    public void testTwoStateProductChain() {
        setUpTwoStates();
        loop("TwoState");
    }

    private void loop(String name) {
        System.out.println("Running " + name + "...");
        int m = 0;
        for (SubstitutionModel substModel : baseModels) {
            double[] out = new double[
                    substModel.getDataType().getStateCount() *
                            substModel.getDataType().getStateCount()];
            substModel.getInfinitesimalMatrix(out);
            System.out.println("R" + m + " = " + new Vector(out));
            m++;
        }

        double[] rates = new double[stateCount * stateCount];
        productChainModel.getInfinitesimalMatrix(rates);
        System.out.println("Product Chain Rate Matrix = " + new Vector(rates));
        assertEquals(rates, markovJumpsInfinitesimalResult, accuracy);

        EigenDecomposition eigen = productChainModel.getEigenDecomposition();
        double[] eval = eigen.getEigenValues();
        double[] sortedEval = new double[eval.length];
        System.arraycopy(eval, 0, sortedEval, 0, stateCount);
        Arrays.sort(sortedEval);
        System.out.println("Eigenvalues = " + new Vector(eigen.getEigenValues()));        
        assertEquals(sortedEval, markovJumpsEigenValues, accuracy);
  
        double[] testProbs = new double[stateCount * stateCount];
        productChainModel.getTransitionProbabilities(0.0, testProbs);
        System.out.println("Finite time (0) probabilities = ");
        printSquareMatrix(testProbs, stateCount);
        double[] trueProbs = new double[stateCount * stateCount];
        for (int i = 0; i < stateCount; i++) {
            trueProbs[i * stateCount + i] = 1.0;
        }
        assertEquals(testProbs, trueProbs, accuracy);

        productChainModel.getTransitionProbabilities(0.5, testProbs);
        System.out.println("Finite time (0.5) probabilities = ");
        printSquareMatrix(testProbs, stateCount);
        assertEquals(testProbs, markovJumpsProbs, accuracy);
    }

    private void printSquareMatrix(double[] A, int dim) {
        double[] row = new double[dim];
        for (int i = 0; i < dim; i++) {
            System.arraycopy(A, i * dim, row, 0, dim);
            System.out.println(new Vector(row));
        }
    }

    List<SubstitutionModel> baseModels;
    int stateCount;

    ProductChainSubstitutionModel productChainModel;
    double[] markovJumpsInfinitesimalResult;
    double[] markovJumpsEigenValues;
    double[] markovJumpsProbs;

    private static final double accuracy = 1E-5;
}

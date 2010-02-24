package test.dr.evomodel.substmodel;

import junit.framework.TestCase;
import dr.app.beagle.evomodel.substmodel.*;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.datatype.DataType;
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

    public void setUp() throws Exception {
        super.setUp();
//        setUpNucleotides();
        setUpTwoStates();
    }

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

        List<DataType> dataTypes = new ArrayList<DataType>();
        dataTypes.add(TwoStates.INSTANCE);
        dataTypes.add(TwoStates.INSTANCE);

        productChainModel = new ProductChainSubstitutionModel("productChain", dataTypes, baseModels);
        stateCount = 2 * 2;

        /*
            model0 = as.eigen.two.state(1.5, 0.75)
            model1 = as.eigen.two.state(2.0, 2.0 / 3.0)
            pc = ind.two.eigen(model0, model1)
            pc$rate.matrix
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
    }

    private void assertEquals(double[] a, double[] b, double accuracy) {
        assertEquals(a.length, b.length);
        for (int i = 0; i < a.length; i++) {
            assertEquals(a[i], b[i], accuracy);
        }
    }

    private void setUpNucleotides() {

        FrequencyModel freqModel0 = new FrequencyModel(Nucleotides.INSTANCE,
                new double[]{0.25, 0.25, 0.25, 0.25});
        FrequencyModel freqModel1 = new FrequencyModel(Nucleotides.INSTANCE,
                new double[]{0.25, 0.25, 0.25, 0.25});
        FrequencyModel freqModel2 = new FrequencyModel(Nucleotides.INSTANCE,
                new double[]{0.25, 0.25, 0.25, 0.25});

        HKY baseModel0 = new HKY(2.0, freqModel0);
        HKY baseModel1 = new HKY(2.0, freqModel1);
        HKY baseModel2 = new HKY(2.0, freqModel2);

        baseModels = new ArrayList<SubstitutionModel>(3);
        baseModels.add(baseModel0);
        baseModels.add(baseModel1);
        baseModels.add(baseModel2);

        List<DataType> dataTypes = new ArrayList<DataType>(3);
        dataTypes.add(Nucleotides.INSTANCE);
        dataTypes.add(Nucleotides.INSTANCE);
        dataTypes.add(Nucleotides.INSTANCE);

        productChainModel = new ProductChainSubstitutionModel("productChain", dataTypes, baseModels);
        stateCount = 4 * 4 * 4;
    }

    public void testInfinitesimalRateMatrices() {
        System.out.println("Running testInfinitesimalRateMatrices...");
        int i = 0;
        for (SubstitutionModel substModel : baseModels) {
            double[] out = new double[
                    substModel.getDataType().getStateCount() *
                            substModel.getDataType().getStateCount()];
            substModel.getInfinitesimalMatrix(out);
            System.out.println("R" + i + " = " + new Vector(out));
            i++;
        }

        double[] out = new double[stateCount * stateCount];
        productChainModel.getInfinitesimalMatrix(out);
        System.out.println("Product Chain Rate Matrix = " + new Vector(out));
        assertEquals(out, markovJumpsInfinitesimalResult, accuracy);

        EigenDecomposition eigen = productChainModel.getEigenDecomposition();
        double[] eval = eigen.getEigenValues();
        Arrays.sort(eval);
        System.out.println("Eigenvalues = " + new Vector(eigen.getEigenValues()));        
        assertEquals(eigen.getEigenValues(), markovJumpsEigenValues, accuracy);
    }

    List<SubstitutionModel> baseModels;
    int stateCount;

    ProductChainSubstitutionModel productChainModel;
    double[] markovJumpsInfinitesimalResult;
    double[] markovJumpsEigenValues;

    private static final double accuracy = 1E-5;
}

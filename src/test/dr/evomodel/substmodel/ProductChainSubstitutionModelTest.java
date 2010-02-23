package test.dr.evomodel.substmodel;

import junit.framework.TestCase;
import dr.app.beagle.evomodel.substmodel.*;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.datatype.DataType;
import dr.math.matrixAlgebra.Vector;

import java.util.List;
import java.util.ArrayList;

/**
 * @author Marc A. Suchard
 */

public class ProductChainSubstitutionModelTest extends TestCase {

    public void setUp() throws Exception {
        super.setUp();

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
//        EigenDecomposition eigen =
                productChainModel.getEigenDecomposition();

        double[] rateMatrix = new double[stateCount * stateCount];
        productChainModel.getInfinitesimalMatrix(rateMatrix);
        System.err.println("Rate Matrix = " + new Vector(rateMatrix));

    }

    List<SubstitutionModel> baseModels;
    int stateCount;

    ProductChainSubstitutionModel productChainModel;
}

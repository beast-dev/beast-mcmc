package dr.app.beagle.evomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.math.KroneckerOperation;
import dr.math.matrixAlgebra.Vector;

import java.util.List;
import java.util.Arrays;

/**
 * @author Marc A. Suchard
 * @author Vladimir Minin
 *         <p/>
 *         A class for implementing a kronecker sum of CTMC models in BEAST using BEAGLE
 *         This work is supported by NSF grant 0856099
 *         <p/>
 *         O'Brien JD, Minin VN and Suchard MA (2009) Learning to count: robust estimates for labeled distances between
 *         molecular sequences. Molecular Biology and Evolution, 26, 801-814
 */
public class ProductChainSubstitutionModel extends BaseSubstitutionModel {

    public ProductChainSubstitutionModel(String name, List<DataType> dataTypes,
                                         List<SubstitutionModel> baseModels) {

        super(name);

        this.baseModels = baseModels;
        numBaseModel = baseModels.size();

        if (numBaseModel == 0) {
            throw new RuntimeException("May not construct ProductChainSubstitutionModel with 0 base models");
        }

        if (numBaseModel != dataTypes.size()) {
            throw new RuntimeException("Each SubstitutionModel must have a DataType in ProductChainSubstitutionModel");
        }

        stateSizes = new int[numBaseModel];
        stateCount = 1;
        for (int i = 0; i < numBaseModel; i++) {
            DataType dataType = dataTypes.get(i);
            stateSizes[i] = dataType.getStateCount();
            stateCount *= dataType.getStateCount();
        }

        updateMatrix = true;
    }

    public EigenDecomposition getEigenDecomposition() {
        synchronized (this) {
            if (updateMatrix) {
                computeKroneckerSumsAndProducts();
            }
        }
        return eigenDecomposition;
    }

// Function 'ind.codon.eigen' from MarkovJumps-R
//  rate.mat = kronecker.sum(kronecker.sum(codon1.eigen$rate.matrix, codon2.eigen$rate.matrix),codon3.eigen$rate.matrix)
//
//  stat = codon1.eigen$stationary%x%codon2.eigen$stationary%x%codon3.eigen$stationary
//
//  ident.vec = rep(1,length(codon1.eigen$stationary))
//
//  eigen.val = (codon1.eigen$values%x%ident.vec + ident.vec%x%codon2.eigen$values)%x%
//    ident.vec + ident.vec%x%ident.vec%x%codon3.eigen$values
//
//  right.eigen.vec = (codon1.eigen$vectors%x%codon2.eigen$vectors)%x%codon3.eigen$vectors
//
//  left.eigen.vec = t((t(codon1.eigen$invvectors)%x%t(codon2.eigen$invvectors))%x%
//    t(codon3.eigen$invvectors))

    public void getInfinitesimalMatrix(double[] out) {
        getEigenDecomposition(); // Updates rate matrix if necessary
        System.arraycopy(rateMatrix, 0, out, 0, stateCount * stateCount);
    }

    private void computeKroneckerSumsAndProducts() {

        int currentStateSize = stateSizes[0];
        double[] currentRate = new double[currentStateSize * currentStateSize];
        baseModels.get(0).getInfinitesimalMatrix(currentRate);
        EigenDecomposition currentED = baseModels.get(0).getEigenDecomposition();
        double[] currentEval = currentED.getEigenValues();

//        double[] ievc = oneEigenDecomposition.getInverseEigenVectors();
//        double[] evec = oneEigenDecomposition.getEigenVectors();

        for (int i = 1; i < numBaseModel; i++) {
            SubstitutionModel nextModel = baseModels.get(i);
            int nextStateSize = stateSizes[i];
            double[] nextRate = new double[nextStateSize * nextStateSize];
            nextModel.getInfinitesimalMatrix(nextRate);
            currentRate = KroneckerOperation.sum(currentRate, currentStateSize, nextRate, nextStateSize);

            EigenDecomposition nextED = nextModel.getEigenDecomposition();
            double[] nextEval = nextED.getEigenValues();
//            System.err.println("cEval = " + new Vector(currentEval));
//            System.err.println("nEval = " + new Vector(nextEval));
//            System.exit(-1);
            currentEval = KroneckerOperation.sum(currentEval, nextEval);
//            System.err.println("final = " + new Vector(currentEval));
//            System.exit(-1);

//            System.err.println("nextStateSize = " + nextStateSize);
//
//            System.err.println("\nCalling with:");
//            System.err.println("currentRate.length = " + currentRate.length);
//            System.err.println("currentSize = " + currentStateSize);
//            System.err.println("nextRate.length = "+ nextRate.length);
//            System.err.println("nextSize = " + nextStateSize);

            currentStateSize *= nextStateSize;
//            System.err.println("Current size = " + currentStateSize);
//            System.exit(-1);
        }

        rateMatrix = currentRate;
        
        eigenDecomposition = new EigenDecomposition(null, null, currentEval);
        updateMatrix = false;
    }

    public FrequencyModel getFrequencyModel() {
        throw new RuntimeException("KroneckerSumSubstitionModel does have a FrequencyModel");
    }

    protected void frequenciesChanged() {
        // Do nothing
    }

    protected void ratesChanged() {
        // Do nothing
    }

    protected void setupRelativeRates(double[] rates) {
        // Do nothing
    }

    private final int numBaseModel;
    private final List<SubstitutionModel> baseModels;
    private final int[] stateSizes;
    //    private final List<DataType> dataTypes;
    private double[] rateMatrix = null;
}

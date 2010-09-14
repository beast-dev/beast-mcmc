package dr.app.beagle.evomodel.substmodel;

import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneralDataType;
import dr.math.KroneckerOperation;
import dr.util.Citable;
import dr.util.Citation;
import dr.util.CommonCitations;
import dr.app.beagle.evomodel.sitemodel.SiteRateModel;
//import dr.math.matrixAlgebra.Vector;

import java.util.List;
import java.util.ArrayList;
import java.util.logging.Logger;

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
public class ProductChainSubstitutionModel extends BaseSubstitutionModel implements Citable {

    public ProductChainSubstitutionModel(String name, List<SubstitutionModel> baseModels) {
        this(name, baseModels, null);
    }

    public ProductChainSubstitutionModel(String name,
                                         List<SubstitutionModel> baseModels,
                                         List<SiteRateModel> rateModels) {
        super(name);

        this.baseModels = baseModels;
        this.rateModels = rateModels;
        numBaseModel = baseModels.size();

        if (numBaseModel == 0) {
            throw new RuntimeException("May not construct ProductChainSubstitutionModel with 0 base models");
        }

        if (rateModels != null) {
            for(SiteRateModel rateModel : rateModels) {
                if (rateModel.getCategoryCount() > 1) {
                    throw new RuntimeException("ProductChainSubstitutionModels with multiple categories not yet implemented");
                }
            }
        }

        List<FrequencyModel> freqModels = new ArrayList<FrequencyModel>();
        stateSizes = new int[numBaseModel];
        stateCount = 1;
        for (int i = 0; i < numBaseModel; i++) {
            freqModels.add(baseModels.get(i).getFrequencyModel());
            DataType dataType = baseModels.get(i).getDataType();
            stateSizes[i] = dataType.getStateCount();
            stateCount *= dataType.getStateCount();
        }

        pcFreqModel = new ProductChainFrequencyModel("pc",freqModels);

        String[] codeStrings = getCharacterStrings();

        dataType = new GeneralDataType(codeStrings);

        updateMatrix = true;

        Logger.getLogger("dr.app.beagle").info("\tConstructing a product chain substition model,  please cite:\n"
                + Citable.Utils.getCitationString(this));
    }

    public List<Citation> getCitations() {
        List<Citation> citations = new ArrayList<Citation>();
        citations.add(
                CommonCitations.OBRIEN_2009
        );
        return citations;
    }

    public EigenDecomposition getEigenDecomposition() {
        synchronized (this) {
            if (updateMatrix) {
                computeKroneckerSumsAndProducts();
            }
        }
        return eigenDecomposition;
    }

    private String[] getCharacterStrings() {
        String[] strings = null;
        for (int i = numBaseModel - 1; i >= 0; i--) {
            strings = recursivelyAppendCharacterStates(baseModels.get(i).getDataType(), strings);
        }

        return strings;
    }

    private String[] recursivelyAppendCharacterStates(DataType dataType, String[] inSubStates) {

        String[] subStates = inSubStates;
        if (subStates == null) {
            subStates = new String[]{""};
        }

        final int previousStateCount = subStates.length;
        final int inStateCount = dataType.getStateCount();
        String[] states = new String[previousStateCount * inStateCount];

        for (int i = 0; i < inStateCount; i++) {
            String code = dataType.getCode(i);
            for (int j = 0; j < previousStateCount; j++) {
                states[i * previousStateCount + j] = code + subStates[j];
            }
        }
        return states;
    }

    public void getInfinitesimalMatrix(double[] out) {
        getEigenDecomposition(); // Updates rate matrix if necessary
        System.arraycopy(rateMatrix, 0, out, 0, stateCount * stateCount);
    }

    protected double[] scaleForProductChain(double[] in, int model) {
        if (rateModels == null) {
            return in;
        }
        final double scalar = rateModels.get(model).getRateForCategory(0);
        if (scalar == 1.0) {
            return in;
        }
        final int len = in.length;
        double[] out = new double[len];
        for (int i = 0; i < len; i++) {
            out[i] = scalar * in[i];
        }
        return out;
    }    

    private void computeKroneckerSumsAndProducts() {

        int currentStateSize = stateSizes[0];
        double[] currentRate = new double[currentStateSize * currentStateSize];
        baseModels.get(0).getInfinitesimalMatrix(currentRate);
        currentRate = scaleForProductChain(currentRate, 0);
        
        EigenDecomposition currentED = baseModels.get(0).getEigenDecomposition();
        double[] currentEval = scaleForProductChain(currentED.getEigenValues(), 0);
        double[] currentEvec = currentED.getEigenVectors();
        double[] currentIevcT = transpose(currentED.getInverseEigenVectors(), currentStateSize);

        for (int i = 1; i < numBaseModel; i++) {
            SubstitutionModel nextModel = baseModels.get(i);
            int nextStateSize = stateSizes[i];
            double[] nextRate = new double[nextStateSize * nextStateSize];
            nextModel.getInfinitesimalMatrix(nextRate);
            nextRate = scaleForProductChain(nextRate, i);
            currentRate = KroneckerOperation.sum(currentRate, currentStateSize, nextRate, nextStateSize);

            EigenDecomposition nextED = nextModel.getEigenDecomposition();
            double[] nextEval = scaleForProductChain(nextED.getEigenValues(), i);
            double[] nextEvec = nextED.getEigenVectors();
            double[] nextIevcT = transpose(nextED.getInverseEigenVectors(), nextStateSize);

            currentEval = KroneckerOperation.sum(currentEval, nextEval);

            currentEvec = KroneckerOperation.product(
                    currentEvec, currentStateSize, currentStateSize,
                    nextEvec, nextStateSize, nextStateSize);

            currentIevcT = KroneckerOperation.product(
                    currentIevcT, currentStateSize, currentStateSize,
                    nextIevcT, nextStateSize, nextStateSize);
            currentStateSize *= nextStateSize;

        }

        rateMatrix = currentRate;

        eigenDecomposition = new EigenDecomposition(
                currentEvec,
                transpose(currentIevcT, currentStateSize),
                currentEval);
        updateMatrix = false;
    }

//   private static void printSquareMatrix(double[] A, int dim) {
//        double[] row = new double[dim];
//        for (int i = 0; i < dim; i++) {
//            System.arraycopy(A, i * dim, row, 0, dim);
//            System.err.println(new Vector(row));
//        }
//    }

    // transposes a square matrix
    private static double[] transpose(double[] mat, int dim) {
        double[] out = new double[dim * dim];
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++) {
                out[j * dim + i] = mat[i * dim + j];
            }
        }
        return out;
    }

    public FrequencyModel getFrequencyModel() {
        return pcFreqModel;
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

    protected final int numBaseModel;
    protected final List<SubstitutionModel> baseModels;
    protected final List<SiteRateModel> rateModels;
    protected final int[] stateSizes;
    protected final ProductChainFrequencyModel pcFreqModel;
    protected double[] rateMatrix = null;
}

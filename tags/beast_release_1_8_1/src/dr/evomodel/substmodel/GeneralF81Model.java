package dr.evomodel.substmodel;

import dr.evolution.datatype.GeneralDataType;

/**
 * @author Alexei Drummond
 */
public class GeneralF81Model extends AbstractSubstitutionModel {

    /**
     * @param freqModel
     */
    public GeneralF81Model(FrequencyModel freqModel) {
        super("generalF81", freqModel.getDataType(), freqModel);
        setupMatrix();
    }

    protected void frequenciesChanged() {
        // frequencyModel changed
    }

    @Override
    protected void ratesChanged() {
        //To change body of implemented methods use File | Settings | File Templates.
    }

    /**
     * get the complete transition probability matrix for the given distance
     *
     * @param distance the expected number of substitutions
     * @param matrix   an array to store the matrix
     */
    public void getTransitionProbabilities(double distance, double[] matrix) {

        int stateSize = freqModel.getFrequencyCount();

        double[] pi = freqModel.getFrequencies();

        double beta = 1;
        for (double p : pi) {
            beta -= p * p;
        }
        beta = 1.0 / beta;

        int c = 0;
        for (int i = 0; i < stateSize; i++) {
            for (int j = 0; j < stateSize; j++) {

                if (i == j) {
                    matrix[c] = pi[i] + (1 - pi[i]) * Math.exp(-beta * distance);
                } else {
                    matrix[c] = pi[j] * (1 - Math.exp(-beta * distance));
                }
                c += 1;
            }
        }
    }

    protected void setupRelativeRates() {
    }

    // *****************************************************************
    // Interface Model
    // *****************************************************************

    public void storeState() {
        super.storeState();
    }

    /**
     * Restore the stored state
     */
    public void restoreState() {
        super.restoreState();
    }

    // **************************************************************
    // XHTMLable IMPLEMENTATION
    // **************************************************************

    public String toXHTML() {
        StringBuffer buffer = new StringBuffer();

        buffer.append("<em>General F81 Model</em>");

        return buffer.toString();
    }

    public static void main(String[] args) {

        GeneralDataType general = new GeneralDataType(new String[]{"0", "1", "2"});

        FrequencyModel freqModel = new FrequencyModel(general, new double[]{0.2, 0.3, 0.5});

        GeneralF81Model f81 = new GeneralF81Model(freqModel);

        int S = general.getStateCount();

        double[] P = new double[S * S];

        f81.getTransitionProbabilities(0.01, P);

        int c = 0;
        for (int i = 0; i < S; i++) {
            System.out.print(P[c]);
            double rowSum = P[c];
            c += 1;
            for (int j = 1; j < S; j++) {
                System.out.print(", " + P[c]);
                rowSum += P[c];
                c += 1;
            }
            System.out.println(" : " + rowSum);
        }

    }
}


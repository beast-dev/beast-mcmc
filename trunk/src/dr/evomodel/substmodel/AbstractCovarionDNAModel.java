package dr.evomodel.substmodel;

import dr.evolution.datatype.HiddenNucleotides;
import dr.inference.model.Parameter;

/**
 * A general time reversible model of nucleotide evolution with
 * covarion hidden rate categories.
 *
 * @author Alexei Drummond
 * @version $Id$
 */
abstract public class AbstractCovarionDNAModel extends AbstractSubstitutionModel {

    public static final String HIDDEN_CLASS_RATES = "hiddenClassRates";
    public static final String SWITCHING_RATES = "switchingRates";
    public static final String FREQUENCIES = "frequencies";

    /**
     * @param name             the name of the covarion substitution model
     * @param dataType         the data type
     * @param freqModel        the equlibrium frequencies
     * @param hiddenClassRates the relative rates of the hidden categories
     *                         (first hidden category has rate 1.0 so this parameter
     *                         has dimension one less than number of hidden categories.
     *                         each hidden category.
     * @param switchingRates   rate of switching between hidden categories
     */
    public AbstractCovarionDNAModel(String name,
                                    HiddenNucleotides dataType,
                                    Parameter hiddenClassRates,
                                    Parameter switchingRates,
                                    FrequencyModel freqModel) {

        super(name, dataType, freqModel);

        hiddenClassCount = dataType.getHiddenClassCount();

        this.hiddenClassRates = hiddenClassRates;
        this.switchingRates = switchingRates;

        assert hiddenClassRates.getDimension() == hiddenClassCount - 1;

        int hiddenClassCount = getHiddenClassCount();

        int switchingClassCount = hiddenClassCount * (hiddenClassCount - 1) / 2;

        if (switchingRates.getDimension() != switchingClassCount) {
            throw new IllegalArgumentException("switching rate parameter must have " +
                    switchingClassCount + " rates for " + hiddenClassCount + " classes");
        }
        addParameter(switchingRates);
        addParameter(hiddenClassRates);
    }

    /**
     * @return the relative rates of A<->C, A<->G, A<->T, C<->G, C<->T and G<->T substitutions
     */
    abstract double[] getRelativeDNARates();

    /**
     * @return the number of hidden classes in this covarion model.
     */
    public final int getHiddenClassCount() {
        return hiddenClassCount;
    }

    public void frequenciesChanged() {
        // DO NOTHING
    }

    public void ratesChanged() {
        setupRelativeRates();
    }

    protected void setupRelativeRates() {

        double[] phi = switchingRates.getParameterValues();
        double[] rr = getRelativeDNARates();

        int x = 0;
        int y = 0;
        for (int i = 0; i < hiddenClassCount; i++) {
            double hiddenRate = hiddenClassRates.getParameterValue(i);
            for (int j = i; j < hiddenClassCount; j++) {

                if (i == j) {
                    // within hidden rate class
                    for (int k = 0; k < rr.length; i++) {
                        relativeRates[x] = rr[k] * hiddenRate;
                        x += 1;
                    }
                } else {
                    // between two hidden rate classes
                    relativeRates[x] = phi[y];
                    relativeRates[x + 1] = 0;
                    relativeRates[x + 2] = 0;
                    relativeRates[x + 3] = 0;
                    relativeRates[x + 4] = phi[y];
                    relativeRates[x + 5] = 0;
                    relativeRates[x + 6] = 0;
                    relativeRates[x + 7] = 0;
                    relativeRates[x + 8] = phi[y];
                    relativeRates[x + 9] = 0;
                    relativeRates[x + 10] = 0;
                    relativeRates[x + 11] = 0;
                    relativeRates[x + 12] = phi[y];
                    relativeRates[x + 13] = 0;
                    relativeRates[x + 14] = 0;
                    relativeRates[x + 15] = 0;
                    x += 16;
                    y += 1;
                }

            }
        }
    }

    /**
     * Normalize rate matrix to one expected substitution per unit time
     *
     * @param matrix the matrix to normalize to one expected substitution
     * @param pi     the equilibrium distribution of states
     */
    void normalize(double[][] matrix, double[] pi) {
        double subst = 0.0;
        int dimension = pi.length;

        for (int i = 0; i < dimension; i++) {
            subst += -matrix[i][i] * pi[i];
        }

        // normalize, including switches
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                matrix[i][j] = matrix[i][j] / subst;
            }
        }

        double switchingProportion = 0.0;
        for (int i = 0; i < hiddenClassCount; i++) {
            for (int j = i + 1; j < hiddenClassCount; j++) {
                for (int l = 0; l < 4; l++) {
                    switchingProportion += matrix[i * 4 + l][j * 4 + l] * pi[j * 4 + l];
                    switchingProportion += matrix[j * 4 + l][i * 4 + l] * pi[i * 4 + l];
                }
            }
        }

        // normalize, removing switches
        for (int i = 0; i < dimension; i++) {
            for (int j = 0; j < dimension; j++) {
                matrix[i][j] = matrix[i][j] / (1.0 - switchingProportion);
            }
        }

        System.out.println("proportion of changes that are switching = " + switchingProportion);
    }

    Parameter switchingRates;
    Parameter hiddenClassRates;

    private int hiddenClassCount;
}

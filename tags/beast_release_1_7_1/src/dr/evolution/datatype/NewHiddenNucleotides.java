package dr.evolution.datatype;

/**
 * @author Marc A. Suchard
 */
public class NewHiddenNucleotides extends Nucleotides implements HiddenDataType {

    public static final String DESCRIPTION = "hiddenNucleotide";

    public static final NewHiddenNucleotides NUCLEOTIDE_HIDDEN_1 = new NewHiddenNucleotides(1);
    public static final NewHiddenNucleotides NUCLEOTIDE_HIDDEN_2 = new NewHiddenNucleotides(2);

    /**
     * Private constructor - DEFAULT_INSTANCE provides the only instance
     */
    private NewHiddenNucleotides(int hiddenClassCount) {
        super();
        this.hiddenClassCount = hiddenClassCount;
    }

    /**
     * returns an array containing the non-ambiguous states that this state represents.
     */
    public boolean[] getStateSet(int state) {

        boolean[] stateSet = new boolean[stateCount * hiddenClassCount];

        if (!isAmbiguousState(state)) {
            for (int h = 0; h < hiddenClassCount; h++)
                stateSet[h * stateCount + state] = true;
        } else {
            for (int i = 0; i < stateCount; i++) {
                stateSet[i] = true;
            }
        }

        return stateSet;
    }

    public int getStateCount() {
        return stateCount * hiddenClassCount;
    }

    private int hiddenClassCount;

    public int getHiddenClassCount() {
        return hiddenClassCount;
    }
}
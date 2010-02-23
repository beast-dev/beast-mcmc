package dr.inference.markovjumps;

/**
 * @author Marc A. Suchard
 * @author Vladimir Minin
 *         <p/>
 *         An enum for different types of codon labelings in BEAST using BEAGLE
 *         This work is supported by NSF grant 0856099
 *         <p/>
 *         O'Brien JD, Minin VN and Suchard MA (2009) Learning to count: robust estimates for labeled distances between
 *         molecular sequences. Molecular Biology and Evolution, 26, 801-814
 */

public enum CodonLabeling {
    SYN("S"), // synonymous mutations
    NON_SYN("N"); // non-synonymous mutations

    CodonLabeling(String text) {
        this.text = text;
    }

    public String getText() {
        return text;
    }

    public double[][] getRegisterMatrix() {
        return null;
    }

    private final String text;

    public static CodonLabeling parseFromString(String text) {
        for (CodonLabeling scheme : CodonLabeling.values()) {
            if (scheme.getText().compareToIgnoreCase(text) == 0)
                return scheme;
        }
        return null;
    }
}

package dr.inference.markovjumps;

/**
 * @author Marc Suchard
 * @author Vladimir Minin
 *         <p/>
 *         A base class for implementing Markov chain-induced counting processes (markovjumps) in BEAST using BEAGLE
 *         This work is supported by NSF grant 0856099
 *         <p/>
 *         Minin VN and Suchard MA (2008) Counting labeled transitions in continous-time Markov models of evolution.
 *         Journal of Mathematical Biology, 56, 391-412.
 */

public enum MarkovJumpsType {

    HISTORY("history"),
    COUNTS("counts"),
    REWARDS("rewards");

    MarkovJumpsType(String text) {
        this.text = text;
    }

    public String getText() {
        return text;
    }

    private final String text;

    public static MarkovJumpsType parseFromString(String text) {
        for (MarkovJumpsType scheme : MarkovJumpsType.values()) {
            if (scheme.getText().compareToIgnoreCase(text) == 0)
                return scheme;
        }
        return null;
    }
}

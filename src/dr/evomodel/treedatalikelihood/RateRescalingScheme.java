package dr.evomodel.treedatalikelihood;

/**
 * @author Marc A. Suchard
 */
public enum RateRescalingScheme {

    NONE("none"),          // no scaling
    TREE_LENGTH("length"), // rescale to one-unit per tree length (in time)
    TREE_HEIGHT("height"); // rescale to one-unit per tree height (in time)

    RateRescalingScheme(String text) {
        this.text = text;
    }

    public String getText() {
        return text;
    }

    private final String text;

    public static RateRescalingScheme parseFromString(String text) {
        for (RateRescalingScheme scheme : RateRescalingScheme.values()) {
            if (scheme.getText().compareToIgnoreCase(text) == 0)
                return scheme;
        }
        return NONE;
    }

    @Override
    public String toString() {
        return text;
    }
}

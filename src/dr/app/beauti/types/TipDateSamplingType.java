package dr.app.beauti.types;

/**
 * @author Alexei Drummond
 */
public enum TipDateSamplingType {

    NO_SAMPLING("Off"),
    SAMPLE_INDIVIDUALLY("Sampling with individual priors"),
    SAMPLE_JOINT("Sampling with joint priors");

    TipDateSamplingType(String displayName) {
        this.displayName = displayName;
    }

    public String toString() {
        return displayName;
    }

    private final String displayName;
}
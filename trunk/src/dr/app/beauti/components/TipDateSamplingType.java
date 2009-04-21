package dr.app.beauti.components;

/**
 * @author Alexei Drummond
 */
public enum TipDateSamplingType {

    NO_SAMPLING("Off"),
    SAMPLE_SET("Sampling a set of tip dates"),
    SAMPLE_ALL("Sampling all tip dates");

    TipDateSamplingType(String displayName) {
        this.displayName = displayName;
    }

    public String toString() {
        return displayName;
    }

    private final String displayName;
}
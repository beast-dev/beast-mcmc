package dr.app.beauti.options;

/**
 * @author Alexei Drummond
 */
public enum FrequencyPolicy {

    ESTIMATED("Estimated"), EMPIRICAL("Empirical"), ALLEQUAL("All equal");

    FrequencyPolicy(String displayName) {
        this.displayName = displayName;
    }

    public String toString() {
        return displayName;
    }

    private String displayName;
}

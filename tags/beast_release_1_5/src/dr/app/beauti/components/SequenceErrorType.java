package dr.app.beauti.components;

/**
 * @author Alexei Drummond
 */
public enum SequenceErrorType {

    NO_ERROR("Off"),
    AGE_TRANSITIONS("Age-dependent: Transistions only"),
    AGE_ALL("Age-dependent: All substitutions"),
    BASE_TRANSITIONS("Age-independent: Transistions only"),
    BASE_ALL("Age-independent: All substitutions");

    SequenceErrorType(String displayName) {
        this.displayName = displayName;
    }

    public String toString() {
        return displayName;
    }

    private final String displayName;
}
package dr.app.beauti.types;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public enum SequenceErrorType {

    NO_ERROR("Off"),
    AGE_TRANSITIONS("Age-dependent: Transitions only"),
    AGE_ALL("Age-dependent: All substitutions"),
    BASE_TRANSITIONS("Age-independent: Transitions only"),
    BASE_ALL("Age-independent: All substitutions"),
    HYPERMUTATION_HA3G("RT Hypermutation: hA3G"),
    HYPERMUTATION_HA3F("RT Hypermutation: hA3F"),
    HYPERMUTATION_BOTH("RT Hypermutation: hA3G + hA3F"),
    HYPERMUTATION_ALL("RT Hypermutation: G->A");

    SequenceErrorType(String displayName) {
        this.displayName = displayName;
    }

    public String toString() {
        return displayName;
    }

    private final String displayName;
}
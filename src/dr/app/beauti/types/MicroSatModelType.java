package dr.app.beauti.types;

/**
 * @author Alexei Drummond
 * @author Walter Xie
 */
public enum MicroSatModelType {

    ASYM_QUAD_MODEL("Asymmetric Quadratic Model");

    MicroSatModelType (String displayName) {
        this.displayName = displayName;
    }

    public String toString() {
        return displayName;
    }

    private String displayName;
}

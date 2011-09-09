package dr.app.beauti.types;

/**
 * @author Andrew Rambaut
 */
public enum ContinuousSubstModelType {
    HOMOGENOUS("Homogenous Brownian model"),
    CAUCHY_RRW("Cauchy RRW model"),
    GAMMA_RRW("GAMMA RRW model");

    ContinuousSubstModelType(String name) {
        this.name = name;
    }

    public String toString() {
        return name;
    }

    private final String name;
}

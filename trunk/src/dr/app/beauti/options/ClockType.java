package dr.app.beauti.options;

/**
 * @author Alexei Drummond
 */
public enum ClockType {

    STRICT_CLOCK("Strict Clock"),
    UNCORRELATED_EXPONENTIAL("Relaxed Clock: Uncorrelated Exp"),
    UNCORRELATED_LOGNORMAL("Relaxed Clock: Uncorrelated Lognormal"),
    RANDOM_LOCAL_CLOCK("Random local clock model");

    ClockType(String name) {
        this.name = name;
    }

    public String toString() {
        return name;
    }

    private final String name;

}

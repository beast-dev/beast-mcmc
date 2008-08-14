package dr.app.beauti.options;

/**
 * @author Alexei Drummond
 */
public enum StartingTreeType {

    RANDOM("randomly generated"),
    USER("user-specified"),
    UPGMA("UPGMA generated");

    StartingTreeType(String name) {
        this.name = name;
    }

    public String getName() {
        return name;
    }

    public String toString() {
        return name;
    }

    private final String name;
}
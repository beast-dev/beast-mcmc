package dr.evolution.datatype;

/**
 * @author Alexei Drummond
 * @author Walter Xie
 */
public enum PloidyType {

    AUTOSOMAL_NUCLEAR("autosomal nuclear", 2.0),
    X("X", 1.5),
    Y("Y", 0.5),
    MITOCHONDRIAL("mitochondrial", 0.5);
    
    PloidyType(String name, double value) {
        this.name = name;
        this.value = value;
    }

    public String toString() {
        return name;
    }
    
    public double getValue() {
        return value;
    }

    private final String name;
    private final double value;
}

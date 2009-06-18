package dr.app.beauti.options;

/**
 * @author Alexei Drummond
 */
public enum NucModelType {

    JC, HKY, GTR;

    public final String getXMLName() {
        return name() + "Model";
    }
}

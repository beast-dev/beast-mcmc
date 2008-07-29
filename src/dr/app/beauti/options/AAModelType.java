package dr.app.beauti.options;

/**
 * @author Alexei Drummond
 */
public enum AAModelType {

    BLOSUM_62("Blosum62", "blosum62"),
    DAYHOFF("Dayhoff", "dayhoff"),
    JTT("JTT"),
    MT_REV_24("mtREV"),
    CP_REV_45("cpREV"),
    WAG("WAG");

    AAModelType(String displayName) {
        this(displayName, displayName);
    }

    AAModelType(String displayName, String xmlName) {
        this.displayName = displayName;
        this.xmlName = xmlName;
    }

    public String toString() {
        return displayName;
    }


    public String getXMLName() {
        return xmlName;
    }

    public static String[] xmlNames() {

        AAModelType[] values = values();

        String[] xmlNames = new String[values.length];
        for (int i = 0; i < values.length; i++) {
            xmlNames[i] = values[i].getXMLName();
        }
        return xmlNames;
    }

    String displayName, xmlName;
}

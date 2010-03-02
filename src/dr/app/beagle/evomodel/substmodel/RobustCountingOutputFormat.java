package dr.app.beagle.evomodel.substmodel;

/**
 * @author Marc A. Suchard
 * @author Vladimir Minin
 *         <p/>
 *         A class for enumerating different robust counting output formats in BEAST using BEAGLE
 *         This work is supported by NSF grant 0856099
 */

public enum RobustCountingOutputFormat {
    SUM_OVER_SITES("sumOverAllSites"),
    PER_SITE("perSite"),
    ARBITRARY_SITES("arbitrarySites");

    RobustCountingOutputFormat(String text) {
        this.text = text;
    }

    public String getText() {
        return text;
    }

    public static RobustCountingOutputFormat parseFromString(String text) {
        for (RobustCountingOutputFormat format : RobustCountingOutputFormat.values()) {
            if (format.getText().compareToIgnoreCase(text) == 0)
                return format;
        }
        return null;
    }

    private final String text;
}

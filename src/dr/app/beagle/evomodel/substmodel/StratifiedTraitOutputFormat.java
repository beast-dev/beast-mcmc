package dr.app.beagle.evomodel.substmodel;

/**
 * @author Marc A. Suchard
 * @author Vladimir Minin
 *         <p/>
 *         A class for enumerating different robust counting output formats in BEAST using BEAGLE
 *         This work is supported by NSF grant 0856099
 */

public enum StratifiedTraitOutputFormat {
    SUM_OVER_SITES("sumOverAllSites"),
    PER_SITE("perSite"),
    PER_SITE_WITH_UNCONDITIONED("perSiteWithUnconditioned"),
    ARBITRARY_SITES("arbitrarySites"),
    ARBITRARY_SITES_WITH_UNCONDITIONED("arbitrarySitesWithUnconditioned");

    StratifiedTraitOutputFormat(String text) {
        this.text = text;
    }

    public String getText() {
        return text;
    }

    public static StratifiedTraitOutputFormat parseFromString(String text) {
        for (StratifiedTraitOutputFormat format : StratifiedTraitOutputFormat.values()) {
            if (format.getText().compareToIgnoreCase(text) == 0)
                return format;
        }
        return null;
    }

    private final String text;
}

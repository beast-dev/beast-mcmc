package dr.app.beagle.evomodel.substmodel;

/**
 * @author Marc A. Suchard
 * @author Vladimir Minin
 *         <p/>
 *         A class for enumerating different robust counting output formats in BEAST using BEAGLE
 *         This work is supported by NSF grant 0856099
 */

public enum StratifiedTraitOutputFormat {
    SUM_OVER_SITES("sumOverAllSites", false),
    SUM_OVER_SITES_WITH_CONDITIONED("sumOverAllSitesWithUnconditioned", true),
    PER_SITE("perSite", false),
    PER_SITE_WITH_UNCONDITIONED("perSiteWithUnconditioned", true),
    ARBITRARY_SITES("arbitrarySites", false),
    ARBITRARY_SITES_WITH_UNCONDITIONED("arbitrarySitesWithUnconditioned", true);

    private StratifiedTraitOutputFormat(String text, boolean supportsUnconditioned) {
        this.text = text;
        this.supportsUnconditioned = supportsUnconditioned;
    }

    public String getText() {
        return text;
    }

    public boolean getSupportsUnconditioned() {
        return supportsUnconditioned;
    }

    public static StratifiedTraitOutputFormat parseFromString(String text) {
        for (StratifiedTraitOutputFormat format : StratifiedTraitOutputFormat.values()) {
            if (format.getText().compareToIgnoreCase(text) == 0)
                return format;
        }
        return null;
    }

    private final String text;
    private final boolean supportsUnconditioned;
}

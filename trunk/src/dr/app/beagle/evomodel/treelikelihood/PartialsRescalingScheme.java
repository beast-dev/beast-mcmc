package dr.app.beagle.evomodel.treelikelihood;

/**
 * @author Marc Suchard
 * @author Andrew Rambaut
 */
public enum PartialsRescalingScheme {

    DEFAULT("default"),
    NONE("none"),
    ALWAYS_RESCALE("alwaysRescale"),
    STATIC_RESCALING("staticRescaling"),
    DYNAMIC_RESCALING("dynamicRescaling");

    PartialsRescalingScheme(String text) {
        this.text = text;
    }

    public String getText() {
        return text;
    }

    private final String text;

    public static PartialsRescalingScheme parseFromString(String text) {
        for(PartialsRescalingScheme scheme : PartialsRescalingScheme.values()) {
            if (scheme.getText().compareToIgnoreCase(text) == 0)
                return scheme;
        }
        return null;
    }

}

/*
 * StochasticMapSampler.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent.mascot;

/** How {@link MascotLikelihood} draws node states and branch histories. */
public enum StochasticMapSampler {

    /**
     * Exact draws from Pi_h, the frozen-environment surrogate on MASCOT's RK4 grid: grid forward filtering and
     * backward sampling with per-step uniformized bridges. The default.
     */
    GRID("grid"),

    /**
     * The Doob h-transform of BEAST2's MappedMascot on cubic Hermite interpolants of MASCOT's lineage
     * probabilities, with exact thinning; approximate for the continuous surrogate, and inaccurate where RK4 steps
     * are stiff (MascotStochasticMap.getMaxRateTimesLength() well above one).
     */
    H_TRANSFORM("hTransform");

    private final String xmlName;

    StochasticMapSampler(String xmlName) {
        this.xmlName = xmlName;
    }

    /** The value of the {@code sampler} attribute of {@code <stochasticMapping>}. */
    public String getXmlName() {
        return xmlName;
    }

    /** Case-insensitive lookup by XML name. */
    public static StochasticMapSampler fromXmlName(String name) {
        for (StochasticMapSampler sampler : values()) {
            if (sampler.xmlName.equalsIgnoreCase(name)) {
                return sampler;
            }
        }
        StringBuilder known = new StringBuilder();
        for (StochasticMapSampler sampler : values()) {
            known.append(known.length() == 0 ? "" : ", ").append('"').append(sampler.xmlName).append('"');
        }
        throw new IllegalArgumentException("unknown stochastic-mapping sampler \"" + name + "\"; use " + known);
    }
}

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
    GRID,

    /**
     * The Doob h-transform of BEAST2's MappedMascot on cubic Hermite interpolants of MASCOT's lineage
     * probabilities, with exact thinning; approximate for the continuous surrogate.
     */
    H_TRANSFORM
}

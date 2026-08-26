/*
 * EpochBoundaries.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package dr.evomodel.coalescent;

import dr.inference.model.Parameter;

/**
 * Shared parsing/validation for "a Parameter of strictly-increasing
 * backward-time breakpoints" -- the concept behind both MASCOT's
 * {@code epochTimes} (feeding {@code GenericMascotLikelihoodDelegate}'s epoch-major {@code theta}
 * layout) and BASTA's {@code gridPoints} (feeding {@code
 * PiecewiseConstantPopulationSizeModel}'s per-segment population sizes).
 *
 * The two engines consume breakpoints differently -- MASCOT integrates
 * continuously through them, BASTA computes closed-form per-segment
 * size/integral pairs -- so this class only owns the one thing that actually
 * was duplicated: turning a {@code Parameter} into a validated, sorted
 * {@code double[]} of times, with a consistent error message. It does not
 * unify the two engines' numeric consumption of those times.
 */
public final class EpochBoundaries {

    private EpochBoundaries() {
    }

    /**
     * @return {@code times}' values, validated strictly increasing; an empty
     * array if {@code times} is {@code null}.
     */
    public static double[] validateSortedTimes(Parameter times, String elementName) {
        if (times == null) {
            return new double[0];
        }
        double[] values = times.getParameterValues();
        for (int i = 1; i < values.length; i++) {
            if (!(values[i] > values[i - 1])) {
                throw new IllegalArgumentException(elementName + " values must be strictly increasing " +
                        "in backward time, got " + values[i - 1] + " followed by " + values[i] +
                        " at index " + i);
            }
        }
        return values;
    }

    /**
     * The {@code [0.0, times[0], times[1], ..., +Infinity]} form MASCOT's
     * {@code GenericMascotLikelihoodDelegate} requires: {@code times.length} breakpoints produce
     * {@code times.length + 1} epochs, epoch 0 spanning {@code [0, times[0])}.
     * A {@code null} (or zero-dimension) {@code times} yields a single epoch
     * spanning {@code [0, +Infinity)}.
     */
    public static double[] withSentinels(Parameter times, String elementName) {
        double[] interior = validateSortedTimes(times, elementName);
        double[] boundaries = new double[interior.length + 2];
        boundaries[0] = 0.0;
        System.arraycopy(interior, 0, boundaries, 1, interior.length);
        boundaries[boundaries.length - 1] = Double.POSITIVE_INFINITY;
        return boundaries;
    }
}

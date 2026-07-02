/*
 * EmbeddedOrdinalParameter.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 *
 */

package dr.inference.hmc;

/**
 * Maps a continuous latent coordinate into ordinal states through interval cuts.
 *
 * State {@code k} corresponds to the interval
 * {@code (cuts[k], cuts[k + 1]]}.
 *
 * @author Filippo Monti (powered by OpenAI)
 */
public class EmbeddedOrdinalParameter {

    private final double[] cuts;

    public EmbeddedOrdinalParameter(double[] cuts) {
        if (cuts == null || cuts.length < 2) {
            throw new IllegalArgumentException("At least two cuts are required");
        }

        this.cuts = cuts.clone();
        for (int i = 1; i < cuts.length; i++) {
            if (!(this.cuts[i] > this.cuts[i - 1])) {
                throw new IllegalArgumentException("Cuts must be strictly increasing");
            }
        }
    }

    public int getStateCount() {
        return cuts.length - 1;
    }

    public int getCutCount() {
        return cuts.length;
    }

    public double getCut(int index) {
        if (index < 0 || index >= cuts.length) {
            throw new IllegalArgumentException("Invalid cut index " + index);
        }
        return cuts[index];
    }

    public int getBoundaryIndex(double value) {
        for (int i = 0; i < cuts.length; i++) {
            if (Double.doubleToLongBits(cuts[i]) == Double.doubleToLongBits(value)) {
                return i;
            }
        }
        return -1;
    }

    public double getNextBoundary(double value, double direction) {
        if (direction > 0.0) {
            for (int i = 0; i < cuts.length; i++) {
                if (cuts[i] > value) {
                    return cuts[i];
                }
            }
        } else if (direction < 0.0) {
            for (int i = cuts.length - 1; i >= 0; i--) {
                if (cuts[i] < value) {
                    return cuts[i];
                }
            }
        }
        return Double.NaN;
    }

    public int getStateIndex(double latentValue) {
        // Binary search for the first cut >= latentValue ("lower bound"), since cuts
        // is sorted; state k's interval is (cuts[k], cuts[k + 1]], except state 0's
        // interval is closed at the lower end ([cuts[0], cuts[1]]).
        final int stateCount = getStateCount();
        int lower = 0;
        int upper = cuts.length;
        while (lower < upper) {
            final int mid = (lower + upper) >>> 1;
            if (cuts[mid] < latentValue) {
                lower = mid + 1;
            } else {
                upper = mid;
            }
        }

        if (lower == 0 && cuts[0] == latentValue) {
            return 0;
        }
        if (lower >= 1 && lower <= stateCount) {
            return lower - 1;
        }

        throw new IllegalArgumentException("Value " + latentValue + " is outside the embedding intervals");
    }

    public double getLowerCut(int state) {
        checkState(state);
        return cuts[state];
    }

    public double getUpperCut(int state) {
        checkState(state);
        return cuts[state + 1];
    }

    private void checkState(int state) {
        if (state < 0 || state >= getStateCount()) {
            throw new IllegalArgumentException("Invalid state index " + state);
        }
    }
}

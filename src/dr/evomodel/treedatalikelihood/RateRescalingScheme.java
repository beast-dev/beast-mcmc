/*
 * RateRescalingScheme.java
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

package dr.evomodel.treedatalikelihood;

/**
 * @author Marc A. Suchard
 */
public enum RateRescalingScheme {

    NONE("none"),          // no scaling
    TREE_LENGTH("length"), // rescale to one-unit per tree length (in time)
    TREE_HEIGHT("height"); // rescale to one-unit per tree height (in time)

    RateRescalingScheme(String text) {
        this.text = text;
    }

    public String getText() {
        return text;
    }

    private final String text;

    public static RateRescalingScheme parseFromString(String text) {
        for (RateRescalingScheme scheme : RateRescalingScheme.values()) {
            if (scheme.getText().compareToIgnoreCase(text) == 0)
                return scheme;
        }
        return NONE;
    }

    @Override
    public String toString() {
        return text;
    }
}

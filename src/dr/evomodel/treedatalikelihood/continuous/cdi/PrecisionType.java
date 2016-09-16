/*
 * PrecisionType.java
 *
 * Copyright (c) 2002-2016 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 */

package dr.evomodel.treedatalikelihood.continuous.cdi;

/**
 * @author Marc A. Suchard
 */
public enum PrecisionType {
    SCALAR("proportional scaling per branch", 0),
    MIXED("mixed method", 1),
    FULL("full precision matrix per branch", 2);

    private final int power;
    private final String name;

    PrecisionType(String name, int power) {
        this.name = name;
        this.power = power;
    }

    public String toString() {
        return name;
    }

    public int getPower() {
        return power;
    }

    public int getPartialLength(final int numTraits, final int dimTrait) {
        int precision = 1;
        final int pow = getPower();
        for (int i = 0; i < pow; ++i) {
            precision *= numTraits;
        }
        return numTraits * (dimTrait + precision);
    }
}

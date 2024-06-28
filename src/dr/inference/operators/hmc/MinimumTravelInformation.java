/*
 * MinimumTravelInformation.java
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

package dr.inference.operators.hmc;

public class MinimumTravelInformation {

    final double time;
    final int[] index;
    final AbstractParticleOperator.Type type;

    MinimumTravelInformation(double minTime, int[] minIndex, AbstractParticleOperator.Type type) { //todo: merge multiple constructors
        this.time = minTime;
        this.index = minIndex;
        this.type = type;
    }

    MinimumTravelInformation(double minTime, int[] minIndex) {
        this.time = minTime;
        this.index = minIndex;
        this.type = AbstractParticleOperator.Type.NONE;
    }

    public MinimumTravelInformation(double minTime, int minIndex) {
        this.time = minTime;
        this.index = new int[]{minIndex};
        this.type = AbstractParticleOperator.Type.NONE;
    }

    MinimumTravelInformation(double minTime, int minIndex, AbstractParticleOperator.Type type) {
        this.time = minTime;
        this.index = new int[]{minIndex};
        this.type = type;
    }

    @SuppressWarnings("unused")
    public MinimumTravelInformation(double minTime, int[] minIndex, int ordinal) {
        this(minTime, minIndex, AbstractParticleOperator.Type.castFromInt(ordinal));
    }

    public boolean equals(Object obj) {

        if (obj instanceof MinimumTravelInformation) {
            MinimumTravelInformation rhs = (MinimumTravelInformation) obj;
//            return Math.abs(this.time - rhs.time) < 1E-5 && this.index == rhs.index;
            return this.time == rhs.time && this.index == rhs.index;
        }
        return false;
    }

    public String toString() {
        return "time = " + time + " @ " + index;
    }
}

/*
 * PriorScaleType.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.app.beauti.types;


/**
 * @author Alexei Drummond
 * @author Walter Xie
 */
public enum PriorScaleType {

    NONE,
    TIME_SCALE,
    LOG_TIME_SCALE,
    GROWTH_RATE_SCALE,
    BIRTH_RATE_SCALE,
    SUBSTITUTION_RATE_SCALE,
    LOG_STDEV_SCALE,
    SUBSTITUTION_PARAMETER_SCALE,
    T50_SCALE,
    ROOT_RATE_SCALE,
    LOG_VAR_SCALE,
    ORIGIN_SCALE,
    GAMMA_STAR_BEAST;


//    PriorScaleType(String name) {
//        this.name = name;
//    }
//
//    public String toString() {
//        return name;
//    }
//
//    private final String name;

}

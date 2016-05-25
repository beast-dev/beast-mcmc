/*
 * TreePriorType.java
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
 */
public enum TreePriorType {

    CONSTANT("Coalescent: Constant Size"),
    EXPONENTIAL("Coalescent: Exponential Growth"),
    LOGISTIC("Coalescent: Logistic Growth"),
    EXPANSION("Coalescent: Expansion Growth"),
    SKYGRID("Coalescent: Bayesian SkyGrid"),
    GMRF_SKYRIDE("Coalescent: GMRF Bayesian Skyride"),
    SKYLINE("Coalescent: Bayesian Skyline"),
    EXTENDED_SKYLINE("Coalescent: Extended Bayesian Skyline Plot"),
    YULE("Speciation: Yule Process"),
    YULE_CALIBRATION("Speciation: Calibrated Yule"),
    BIRTH_DEATH("Speciation: Birth-Death Process"),
    BIRTH_DEATH_INCOMPLETE_SAMPLING("Speciation: Birth-Death Incomplete Sampling"),
    BIRTH_DEATH_SERIAL_SAMPLING("Speciation: Birth-Death Serially Sampled"),
    BIRTH_DEATH_BASIC_REPRODUCTIVE_NUMBER("Epidemiology: Birth-Death Basic Reproductive Number"),
    SPECIES_YULE("Species Tree: Yule Process"),
    SPECIES_YULE_CALIBRATION("Species Tree: Calibrated Yule"),
    SPECIES_BIRTH_DEATH("Species Tree: Birth-Death Process");

    TreePriorType(String name) {
        this.name = name;
    }

    public String toString() {
        return name;
    }

    private final String name;
}

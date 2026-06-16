/*
 * TreePriorType.java
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

package dr.app.beauti.types;

/**
 * @author Alexei Drummond
 */
public enum TreePriorType {

    CONSTANT("Coalescent: Constant Size"),
    EXPONENTIAL("Coalescent: Exponential Growth"),
    LOGISTIC("Coalescent: Logistic Growth"),
    EXPANSION("Coalescent: Expansion Growth"),
    //    SKYLINE("Coalescent: Bayesian Skyline"), // @deprecated
    SKYGRID_HMC("Coalescent: Hamiltonian Monte Carlo SkyGrid"),
    SKYGRID("Coalescent: Bayesian SkyGrid"),
    GMRF_SKYRIDE("Coalescent: GMRF Bayesian Skyride"),
    EPISODIC_BIRTH_DEATH_SAMPLING("Epidemiology: Episodic Birth-Death Sampling"),
    BIRTH_DEATH_SERIAL_SAMPLING("Epidemiology: Birth-Death Serially Sampled"),
    //    BIRTH_DEATH_BASIC_REPRODUCTIVE_NUMBER("Epidemiology: Birth-Death Basic Reproductive Number"),  // @deprecated
    YULE("Speciation: Yule Process"),
    YULE_CALIBRATION("Speciation: Calibrated Yule"),
    BIRTH_DEATH("Speciation: Birth-Death Process"),
    BIRTH_DEATH_INCOMPLETE_SAMPLING("Speciation: Birth-Death Incomplete Sampling"),
    SET_BY_BIT("Set by backward-in-time model");

    TreePriorType(String name) {
        this.name = name;
    }

    public String toString() {
        return name;
    }

    private final String name;
}

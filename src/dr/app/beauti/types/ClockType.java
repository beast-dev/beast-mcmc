/*
 * ClockType.java
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
public enum ClockType {

    STRICT_CLOCK("Strict clock"),
    UNCORRELATED("Uncorrelated relaxed clock"),
    HMC_CLOCK("Hamiltonian Monte Carlo relaxed clock"),

    SHRINKAGE_LOCAL_CLOCK("Shrinkage local clock"),
    RANDOM_LOCAL_CLOCK("Classic random local clock"),
    FIXED_LOCAL_CLOCK("Fixed local clock"),
    MIXED_EFFECTS_CLOCK("Mixed effects clock"),
    AUTOCORRELATED("Autocorrelated relaxed clock");

    ClockType(String displayName) {
        this.displayName = displayName;
    }

    public String toString() {
        return displayName;
    }

    private final String displayName;

    //maintain unique parameter IDs to enable correct prior and transition kernel retrieval
    final public static String LOCAL_CLOCK = "localClock";
    final public static String UCED_MEAN = "uced.mean";
    final public static String UCLD_MEAN = "ucld.mean";
    final public static String UCLD_STDEV = "ucld.stdev";
    final public static String UCGD_MEAN = "ucgd.mean";
    final public static String UCGD_SHAPE = "ucgd.shape";
    final public static String SHRINKAGE_CLOCK_LOCATION = "branchRatesShrinkage.rate";
    final public static String HMC_CLOCK_LOCATION = "branchRatesHMC.rate";
    final public static String HMC_CLOCK_BRANCH_RATES = "branchRatesHMC.rates";
    final public static String HMCLN_SCALE = "branchRatesHMC.scale";
    final public static String ME_CLOCK_LOCATION = "branchRatesME.rate";
    final public static String ME_CLOCK_SCALE = "branchRatesME.scale";

    final public static String ACLD_MEAN = "acld.mean";
    final public static String ACLD_STDEV = "acld.stdev";

    final public static String HMC_CLOCK_RATES_DESCRIPTION = "HMC relaxed clock branch rates";
    final public static String HMC_CLOCK_LOCATION_SCALE_DESCRIPTION = "HMC relaxed clock location and scale";
}
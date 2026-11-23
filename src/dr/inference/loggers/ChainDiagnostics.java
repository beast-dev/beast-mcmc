/*
 * ChainDiagnostics.java
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

package dr.inference.loggers;

/**
 * Convergence diagnostics of the chains of a run, as shown by {@link ChainDiagnosticColumn}s.
 * A value that is not available is NaN (and its label null).
 *
 * @author Andrew Rambaut
 */
public interface ChainDiagnostics {

    /**
     * @return the smallest effective sample size (summed over the chains) of the logged values
     */
    double getMinESS();

    /**
     * @return the label of the logged value with the smallest effective sample size
     */
    String getMinESSLabel();

    /**
     * @return the largest split R-hat of the logged values
     */
    double getMaxRhat();

    /**
     * @return the label of the logged value with the largest split R-hat
     */
    String getMaxRhatLabel();

    /**
     * @return the average standard deviation of clade frequencies across the chains
     */
    double getASDSF();

    /**
     * @return the burn-in the chains share, in states
     */
    double getBurnin();

    /**
     * @return the proportion of the samples left out as burn-in when computing the diagnostics
     * (this is never NaN: when the chains share no burn-in yet it is the largest one tried)
     */
    double getBurninProportion();
}

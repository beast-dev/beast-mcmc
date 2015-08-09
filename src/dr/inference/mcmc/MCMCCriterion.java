/*
 * MCMCCriterion.java
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

package dr.inference.mcmc;

import dr.math.MathUtils;
import dr.inference.markovchain.Acceptor;

/**
 * This class encapsulates the acceptance criterion for an MCMC proposal.
 *
 * @author Alexei Drummond
 *
 * @version $Id: MCMCCriterion.java,v 1.12 2005/05/24 20:25:59 rambaut Exp $
 */
public class MCMCCriterion implements Acceptor {

    // this parameter is actually 1/T, when the temperature parameter is 0.0, then the distribution
    // is flat and will always accept (symmetric) proposals, i.e. hastings ratio of 0 in log space.
    // As this temperature parameter increases, the posterior gets more peaked.
    protected double temperature = 1.0;

    public MCMCCriterion() {

        temperature = 1.0;
    }

    public MCMCCriterion(double t) {
        temperature = t;
    }

    public double getAcceptanceValue(double oldScore, double hastingsRatio) {

        final double acceptanceValue =
                (MathUtils.randomLogDouble() + (oldScore * temperature) - hastingsRatio) / temperature;

        return acceptanceValue;
    }

    public boolean accept(double oldScore, double newScore, double hastingsRatio, double[] logr) {

        logr[0] = (newScore - oldScore) * temperature + hastingsRatio;

        // for coercedAcceptanceProbability
        if (logr[0] > 0) logr[0] = 0.0;

        final double v = MathUtils.randomLogDouble();
        final boolean accept = v < logr[0];

        return accept;
    }

    public double getTemperature() {
        return temperature;
    }

    public void setTemperature(double temperature) {
        this.temperature = temperature;
    }

}

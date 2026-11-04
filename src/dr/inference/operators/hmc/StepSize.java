/*
 * StepSize.java
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

/**
 * @author Marc A. Suchard
 */
class StepSize {

    final private Options options;
    final private double mu;

    private double stepSize;
    private double logStepSize;
    private double averageLogStepSize;
    private double h;

    StepSize(double initialStepSize) {
        this(initialStepSize, new Options());
    }

    private StepSize(double initialStepSize, Options options) {
        this.options = options;
        this.mu = Math.log(options.muFactor * initialStepSize);
        this.stepSize = initialStepSize;
        this.logStepSize = Math.log(stepSize);
        this.averageLogStepSize = 0;
        this.h = 0;
    }

    void update(long m, double cumAcceptProb, double numAcceptProbStates) {

        if (m <= options.adaptLength) {

            h = (1 - 1 / (m + options.t0)) * h + 1 / (m + options.t0) * (options.targetAcceptRate - (cumAcceptProb / numAcceptProbStates));
            logStepSize = mu - Math.sqrt(m) / options.gamma * h;
            averageLogStepSize = Math.pow(m, -options.kappa) * logStepSize +
                    (1 - Math.pow(m, -options.kappa)) * averageLogStepSize;
            stepSize = Math.exp(logStepSize);
        }
    }

    double getStepSize() {
        return stepSize;
    }

    static class Options { //TODO: these values might be adjusted for dual averaging.
        double kappa = 0.75;
        double t0 = 10.0;
        double gamma = 0.05;
        double targetAcceptRate = 0.85;
        double muFactor = 10.0;
        int adaptLength = 1000;
    }
}

/*
 * LikelihoodBasedDiscontinuousPotentialProvider.java
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

package dr.inference.hmc;

import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;

/**
 * Generic discontinuous potential provider that evaluates a likelihood after a
 * temporary one-coordinate move. The default path is deliberately eventful:
 * it fires the normal parameter change notifications while probing and while
 * restoring the original value, so listener-driven caches stay correct. The
 * quiet path is a faster debug/toy shortcut and must only be used when the
 * wrapped likelihood is self-contained with respect to the probed parameter.
 *
 * @author Filippo Monti (powered by OpenAI)
 */
public class LikelihoodBasedDiscontinuousPotentialProvider implements DiscontinuousPotentialProvider {

    private final Likelihood likelihood;
    private final Parameter parameter;
    private final boolean quiet;

    public LikelihoodBasedDiscontinuousPotentialProvider(Likelihood likelihood, Parameter parameter) {
        this(likelihood, parameter, false);
    }

    public LikelihoodBasedDiscontinuousPotentialProvider(Likelihood likelihood, Parameter parameter, boolean quiet) {
        this.likelihood = likelihood;
        this.parameter = parameter;
        this.quiet = quiet;
    }

    public Likelihood getLikelihood() {
        return likelihood;
    }

    public boolean isQuiet() {
        return quiet;
    }

    @Override
    public Parameter getParameter() {
        return parameter;
    }

    @Override
    public int getDimension() {
        return parameter.getDimension();
    }

    @Override
    public double getLogDensity() {
        return likelihood.getLogLikelihood();
    }

    @Override
    public double getLogDensityAfterSingleCoordinateMove(int index, double proposedValue) {
        final double currentValue = parameter.getParameterValue(index);

        try {
            if (quiet) {
                // quiet=true is a speed/debug option: it avoids BEAST's normal
                // parameter-change event cascade and only dirties this likelihood.
                // Use it only when no listener-maintained cache depends on the
                // temporary value being probed.
                parameter.setParameterValueQuietly(index, proposedValue);
            } else {
                parameter.setParameterValue(index, proposedValue);
            }
            likelihood.makeDirty();
            return likelihood.getLogLikelihood();
        } finally {
            if (quiet) {
                parameter.setParameterValueQuietly(index, currentValue);
            } else {
                parameter.setParameterValue(index, currentValue);
            }
            likelihood.makeDirty();
        }
    }
}

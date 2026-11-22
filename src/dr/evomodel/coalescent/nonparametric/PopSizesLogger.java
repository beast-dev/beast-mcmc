/*
 * PopSizesLogger.java
 *
 * Copyright © 2002-2026, the BEAST Development Team.
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
 */

package dr.evomodel.coalescent.nonparametric;

import dr.evomodel.coalescent.legacy.GMRFMultilocusSkyrideLikelihood;
import dr.evomodel.coalescent.skyride.GMRFSkygridLikelihood;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.Loggable;
import dr.inference.loggers.NumberColumn;
import dr.inference.model.AbstractModel;
import dr.util.Transform;

public class PopSizesLogger implements Loggable {
//    private final MultilocusNonparametricCoalescentLikelihood likelihood;
    private final Transform transform;
    private final String order;
    private final int popSizeDimension;
    private AbstractModel likelihood;

    public PopSizesLogger(NonParametricCoalescentLikelihood likelihood, Transform transform, String order) {
        this.transform = transform;
        this.order = order; // not yet implemented
        this.popSizeDimension = likelihood.getPopSizeDimension();
        this.likelihood = likelihood;
    }

    // TODO this is a temporary constructor
    public PopSizesLogger(GMRFSkygridLikelihood likelihood, Transform transform, String order) {
        this.transform = transform;
        this.order = order; // not yet implemented
        this.popSizeDimension = likelihood.getDimension();
        this.likelihood = likelihood;
    }

    // TODO this is a temporary constructor
    public PopSizesLogger(GMRFMultilocusSkyrideLikelihood likelihood, Transform transform, String order) {
        this.transform = transform;
        this.order = order; // not yet implemented
        this.popSizeDimension = likelihood.getDimension();
        this.likelihood = likelihood;
    }

    private double getPopSize(int i) {
        if (likelihood instanceof NonParametricCoalescentLikelihood) {
            return ((NonParametricCoalescentLikelihood) likelihood).getPopulationSize(i);
        } else if (likelihood instanceof GMRFSkygridLikelihood) {
            double x = ((GMRFSkygridLikelihood) likelihood).getParameter().getParameterValue(i);
            return Math.exp(x);
        } else if (likelihood instanceof GMRFMultilocusSkyrideLikelihood) {
            double x = ((GMRFMultilocusSkyrideLikelihood) likelihood).getParameter().getParameterValue(i);
            return Math.exp(x);
        } else {
            throw new RuntimeException("Likelihood type not yet implemented.");
        }
    }
    @Override
    public LogColumn[] getColumns() {
        LogColumn[] columns = new LogColumn[popSizeDimension];
        for (int k = 0; k < popSizeDimension; ++k) {
            int finalK = k;
            columns[k] = new NumberColumn("PopSize" + "." + (finalK + 1)) {
                @Override
                public double getDoubleValue() {
                    if (transform != null) {
                        return transform.transform(getPopSize(finalK));
//                        return transform.transform(likelihood.getPopulationSize(finalK));
                    } else {
                        return getPopSize(finalK);
                    }
                }
            };
        }
        return columns;
    }
}

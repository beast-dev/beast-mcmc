/*
 * SampleNonActiveGibbsOperator.java
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

package dr.evomodel.coalescent.operators;

import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.Parameter;
import dr.inference.operators.GibbsOperator;
import dr.inference.operators.SimpleMCMCOperator;
import dr.math.MathUtils;

/**

 */
public class SampleNonActiveGibbsOperator extends SimpleMCMCOperator implements GibbsOperator {
    private final ParametricDistributionModel distribution;
    private final Parameter data;
    private final Parameter indicators;

    public SampleNonActiveGibbsOperator(ParametricDistributionModel distribution,
                                        Parameter data, Parameter indicators, double weight) {
        this.distribution = distribution;
        this.data = data;
        this.indicators = indicators;
        setWeight(weight);
    }


    public String getPerformanceSuggestion() {
        return null;
    }

    public String getOperatorName() {
        return "SampleNonActive(" + indicators.getId() + ")";
    }

    public double doOperation() {
        final int idim = indicators.getDimension();

        final int offset = (data.getDimension() - 1) == idim ? 1 : 0;
        assert offset == 1 || data.getDimension() == idim : "" + idim + " (?+1) != " + data.getDimension();

        // available locations for direct sampling
        int[] loc = new int[idim];
        int nLoc = 0;

        for (int i = 0; i < idim; ++i) {
            final double value = indicators.getStatisticValue(i);
            if (value == 0) {
                loc[nLoc] = i + offset;
                ++nLoc;
            }
        }

        if (nLoc > 0) {
            final int index = loc[MathUtils.nextInt(nLoc)];
            try {
                final double val = distribution.quantile(MathUtils.nextDouble());
                data.setParameterValue(index, val);
            } catch (Exception e) {
                // some distributions fail on extreme values - currently gamma
               return Double.NEGATIVE_INFINITY;
            }
        } else {
//            throw new OperatorFailedException("no non-active indicators");
            return Double.NEGATIVE_INFINITY;
        }
        return 0.0;
    }

    public int getStepCount() {
        return 0;
    }
}

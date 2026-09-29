/*
 * UpDownOperator.java
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

package dr.inference.operators;

import dr.inference.model.Parameter;
import dr.math.MathUtils;

public class FreeRateUpDownOperator extends AbstractAdaptableOperator {

    private Parameter weightsParameter = null;
    private Parameter ratesParameter = null;
    private double scaleFactor;

    public FreeRateUpDownOperator(Parameter weightsParameter, Parameter ratesParameter,
                                  double scale, double weight, AdaptationMode mode) {

        super(mode);
        setWeight(weight);

        this.weightsParameter = weightsParameter;
        this.ratesParameter = ratesParameter;
        this.scaleFactor = scale;

        if(weightsParameter.getDimension() != ratesParameter.getDimension()) {
            throw new IllegalArgumentException("Rate and weight parameters should have the same dimension");
        }
    }

    public final double getScaleFactor() {
        return scaleFactor;
    }

    public final void setScaleFactor(double sf) {
        if( (sf > 0.0) && (sf < 1.0) ) {
            scaleFactor = sf;
        } else {
            throw new IllegalArgumentException("scale must be between 0 and 1");
        }
    }

    /**
     * change the parameter and return the hastings ratio.
     */
    public final double doOperation() {

        final int dim = weightsParameter.getDimension();
        final int dim1 = MathUtils.nextInt(dim);
        int dim2 = dim1;
        while (dim1 == dim2) {
            dim2 = MathUtils.nextInt(dim);
        }

        double weight1 = weightsParameter.getParameterValue(dim1);
        double weight2 = weightsParameter.getParameterValue(dim2);

        double rate1 = ratesParameter.getParameterValue(dim1);
        double rate2 = ratesParameter.getParameterValue(dim2);

        final double scale1 = (scaleFactor + (MathUtils.nextDouble() * ((1.0 / scaleFactor) - scaleFactor)));
        final double scale2 = (1-scale1) * rate1*weight1/(rate2*weight2) + 1;

        rate1 = rate1 * scale1;
        rate2 = rate2 * scale2;

        if(rate1 <=0 | rate2 <=0) {
            return Double.NEGATIVE_INFINITY;
        }

        ratesParameter.setParameterValue(dim1, rate1);
        ratesParameter.setParameterValue(dim2, rate2);

        return -Math.log(scale1);
    }

    public String getAdaptableParameterName() {
        return "scaleFactor";
    }

    public final String getOperatorName() {
        return "freeRateUpDownOperator(weights:" + weightsParameter.getParameterName() + ", rates: " + ratesParameter + ")";
    }

    @Override
    protected double getAdaptableParameterValue() {
        return Math.log(1.0 / scaleFactor - 1.0) / Math.log(10);
    }

    @Override
    public void setAdaptableParameterValue(double value) {
        scaleFactor = 1.0 / (Math.pow(10.0, value) + 1.0);
    }

    @Override
    public double getRawParameter() {
        return scaleFactor;
    }

}



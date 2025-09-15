/*
 * DeltaExchangeOperator.java
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

import dr.inference.model.Bounds;
import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.math.distributions.BetaDistribution;
import dr.math.distributions.NormalDistribution;
import dr.util.Transform;
import org.apache.commons.math.MathException;

/**
 * A generic operator for use with a sum-constrained (possibly weighted) vector parameter.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class FreeRateDeltaExchangeOperator extends AbstractAdaptableOperator {

    public FreeRateDeltaExchangeOperator(Parameter weightsParameter, Parameter ratesParameter, double delta) {

        super(AdaptationMode.ADAPTATION_ON);

        this.weightsParameter = weightsParameter;
        this.ratesParameter = ratesParameter;
        this.delta = delta;
        if(delta <= 0){
            throw new IllegalArgumentException("Delta must be greater than 2");
        }
        setWeight(1.0);

        if(weightsParameter.getDimension() != ratesParameter.getDimension()) {
            throw new IllegalArgumentException("Rate and weight parameters should have the same dimension");
        }
    }

    public FreeRateDeltaExchangeOperator(Parameter weightsParameter, Parameter ratesParameter, double delta, double weight, AdaptationMode mode) {

        super(mode);

        this.weightsParameter = weightsParameter;
        this.ratesParameter = ratesParameter;
        this.delta = delta;
        setWeight(weight);

    }

    /**
     * @return the parameter this operator acts on.
     */
    public Parameter getWeightsParameter() {
        return weightsParameter;
    }

    public Parameter getRatesParameter() {  return ratesParameter; }

    /**
     * change the parameter and return the hastings ratio.
     * performs a delta exchange operation between two scalars in the vector
     * and return the hastings ratio.
     */
    public final double doOperation() {

        // get two dimensions
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

        // logit-transform weight1

        double logNum = Math.log(weight1);

        if(logNum == 1) {
            logNum = Math.log1p(weight1 - 1);
        }

        double logDenom = Math.log(1-weight1);

        if(logDenom == 1){
            logDenom = Math.log1p(-weight1);
        }

        double logitWeight1 = logNum - logDenom;

        NormalDistribution drawDist = new NormalDistribution(logitWeight1, delta);

        double draw = (double) drawDist.nextRandom();
        double newWeight1 = 1/(1+Math.exp(-draw));

        double change = newWeight1 - weight1;
        double newWeight2 = weight2 - change;

        double newRate1 = rate1 + change*(rate2-rate1)/(weight1 + change);

        Bounds<Double> weightsBounds = weightsParameter.getBounds();
        Bounds<Double> ratesBounds = ratesParameter.getBounds();

        if (newWeight1 <= weightsBounds.getLowerLimit(dim1) ||
                newWeight1 >= weightsBounds.getUpperLimit(dim1) ||
                newWeight2 <= weightsBounds.getLowerLimit(dim2) ||
                newWeight2 >= weightsBounds.getUpperLimit(dim2)) {
            return Double.NEGATIVE_INFINITY;
        }

        if (newRate1 < ratesBounds.getLowerLimit(dim1) ||
                newRate1 > ratesBounds.getUpperLimit(dim1)) {
            return Double.NEGATIVE_INFINITY;
        }

        ratesParameter.setParameterValue(dim1, newRate1);
        weightsParameter.setParameterValue(dim1, newWeight1);
        weightsParameter.setParameterValue(dim2, newWeight2);

        // Need to adjust by the Jacobian

        double hr = -(1-(1+Math.exp(-draw))*weight1) + 1;
        hr = hr*Math.exp(-draw);
        hr = hr/Math.pow(1+Math.exp(-draw), 2);
        hr = hr*((1/weight1) + 1/(1-weight1));
        hr = Math.abs(hr);
        double loghr = Math.log(hr);
  ;

        return loghr;
    }

    // Interface MCMCOperator
    public final String getOperatorName() {
        return "freeRateDeltaExchange(weights:" + weightsParameter.getParameterName() + ", rates: " + ratesParameter + ")";
    }

    @Override
    protected double getAdaptableParameterValue() {
        return Math.log(delta);
    }

    public void setAdaptableParameterValue(double value) {
        delta = Math.exp(value);
    }

    public double getRawParameter() {
        return delta;
    }

    public String getAdaptableParameterName() {
        return "delta";
    }

    public String toString() {
        return getOperatorName() + "(windowsize=" + delta + ")";
    }

    // Private instance variables

    private Parameter weightsParameter = null;
    private Parameter ratesParameter = null;
    private double delta = 0.02;
}

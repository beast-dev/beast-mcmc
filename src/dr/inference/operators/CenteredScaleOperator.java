/*
 * CenteredScaleOperator.java
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

package dr.inference.operators;

import dr.inference.model.Parameter;
import dr.math.MathUtils;

/**
 * A generic operator for use with a sum-constrained vector parameter.
 *
 * @author Alexei Drummond
 * @version $Id: CenteredScaleOperator.java,v 1.20 2005/06/14 10:40:34 rambaut Exp $
 */

// AR - I don't see how this can work...
@Deprecated
public class CenteredScaleOperator extends AbstractCoercableOperator {

    public CenteredScaleOperator(Parameter parameter) {
        super(CoercionMode.DEFAULT);
        this.parameter = parameter;
    }

    public CenteredScaleOperator(Parameter parameter, double scale, int weight, CoercionMode mode) {
        super(mode);
        this.parameter = parameter;
        this.scaleFactor = scale;
        setWeight(weight);
    }

    /**
     * @return the parameter this operator acts on.
     */
    public Parameter getParameter() {
        return parameter;
    }

    /**
     * Change the parameter and return the hastings ratio.
     * Performs a centered scale operation on the vector
     * and returns the hastings ratio.
     * This operator changes the variance but maintains the order
     * of the scalars.
     */
    public final double doOperation() {

        double total = 0.0;

        for (int i = 0; i < parameter.getDimension(); i++) {
            total += parameter.getParameterValue(i);
        }
        double mean = total / parameter.getDimension();
        double scaleFactor = getRandomScaleFactor();
        double logq = parameter.getDimension() * Math.log(1.0 / scaleFactor);

        for (int i = 0; i < parameter.getDimension(); i++) {

            double newScalar = (parameter.getParameterValue(i) - mean) * scaleFactor + mean;
            if (newScalar < parameter.getBounds().getLowerLimit(i) || newScalar > parameter.getBounds().getUpperLimit(i)) {
//                throw new OperatorFailedException("Proposed value out of bounds");
                return Double.NEGATIVE_INFINITY;
            }
            parameter.setParameterValue(i, newScalar);
        }

        // non-symmetrical move
        return logq;
    }

    public final double getRandomScaleFactor() {
        return scaleFactor + (MathUtils.nextDouble() * ((1 / scaleFactor) - scaleFactor));
    }

    // Interface MCMCOperator
    public final String getOperatorName() {
        return parameter.getParameterName();
    }

    public double getCoercableParameter() {
        return Math.log(1.0 / scaleFactor - 1.0);
    }

    public void setCoercableParameter(double value) {
        scaleFactor = 1.0 / (Math.exp(value) + 1.0);
    }

    public double getRawParameter() {
        return scaleFactor;
    }

    public double getTargetAcceptanceProbability() {
        return 0.234;
    }

    public final String getPerformanceSuggestion() {

        double prob = MCMCOperator.Utils.getAcceptanceProbability(this);
        double targetProb = getTargetAcceptanceProbability();
        dr.util.NumberFormatter formatter = new dr.util.NumberFormatter(5);
        double sf = OperatorUtils.optimizeScaleFactor(scaleFactor, prob, targetProb);
        if (prob < getMinimumGoodAcceptanceLevel()) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else if (prob > getMaximumGoodAcceptanceLevel()) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else return "";
    }

    public String toString() {
        return getOperatorName() + "(scaleFactor=" + scaleFactor + ")";
    }

    // Private instance variables

    private Parameter parameter = null;
    public double scaleFactor = 0.5;
}

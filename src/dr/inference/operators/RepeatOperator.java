/*
 * Repeat.java
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

import dr.math.Poisson;



/**
 *   An operator adaptable operator that repeats an operation a random
 *   number of times drawn from a shifted Poisson distribution
 * @author jtmccrone
 */
public class RepeatOperator  extends AbstractAdaptableOperator {
    private final SimpleMCMCOperator operator;
    private double nonshiftedMean;
    public RepeatOperator(double weight, double targetAcceptanceProbability,SimpleMCMCOperator operator,double nonshiftedMean) {
        super(AdaptationMode.DEFAULT, targetAcceptanceProbability);
        this.operator = operator;
        this.nonshiftedMean = nonshiftedMean;
        setWeight(weight);
    }

    @Override
    public String getOperatorName() {
        return "Repeat " + operator.getOperatorName();
    }

    public final double doOperation() {
        double logP = 0;
        int count = Poisson.nextPoisson(this.nonshiftedMean) +1;
        for (int i =0;i<count;i++) {
            logP += operator.doOperation();
        }
        return logP;
    }


    /**
     * Sets the adaptable parameter value.
     *
     * @param value the value to set the adaptable parameter to
     */
    @Override
    protected void setAdaptableParameterValue(double value) {
        nonshiftedMean = Math.exp(value);
    }

    /**
     * Gets the adaptable parameter value.
     *
     * @returns the value
     */
    @Override
    protected double getAdaptableParameterValue() {
        return Math.log(nonshiftedMean);
    }


    /**
     * @return the underlying tuning parameter value
     */
    @Override
    public double getRawParameter() {
        return nonshiftedMean;
    }

    @Override
    public String getAdaptableParameterName() {
        return "Repeat" + operator.getOperatorName();
    }
}

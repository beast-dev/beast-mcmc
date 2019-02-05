/*
 * AbstractAdaptableOperator.java
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

import java.util.ArrayDeque;
import java.util.Deque;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public abstract class AbstractAdaptableOperator extends SimpleMCMCOperator implements AdaptableMCMCOperator {

    public final AdaptationMode mode;
    private final double targetAcceptanceProbability;
    private long adaptationCount = 0;

    public AbstractAdaptableOperator() {
        this(AdaptationMode.DEFAULT, DEFAULT_ADAPTATION_TARGET);
    }

    public AbstractAdaptableOperator(AdaptationMode mode) {
        this(mode, DEFAULT_ADAPTATION_TARGET);
    }

    public AbstractAdaptableOperator(AdaptationMode mode, double targetAcceptanceProbability) {
        this.mode = mode;
        if (System.getProperty("mcmc.adaptation_target") != null) {
            this.targetAcceptanceProbability = Double.parseDouble(System.getProperty("mcmc.adaptation_target"));
        } else {
            this.targetAcceptanceProbability = targetAcceptanceProbability;
        }
    }

    @Override
    public final double getTargetAcceptanceProbability() {
        return targetAcceptanceProbability;
    }

    @Override
    public final long getAdaptationCount() {
        return adaptationCount;
    }

    @Override
    public final void setAdaptableParameter(double value) {
        setAdaptableParameterValue(value);
        adaptationCount ++;
    }

    @Override
    public final double getAdaptableParameter() {
        return getAdaptableParameterValue();
    }

    /**
     * Sets the adaptable parameter value.
     *
     * @param value the value to set the adaptable parameter to
     */
    protected abstract void setAdaptableParameterValue(double value);

    /**
     * Gets the adaptable parameter value.
     *
     * @returns the value
     */
    protected abstract double getAdaptableParameterValue();


    public double getMinimumAcceptanceLevel() {
        return MINIMUM_ACCEPTANCE_LEVEL;
    }

    public double getMaximumAcceptanceLevel() {
        return MAXIMUM_ACCEPTANCE_LEVEL;
    }

    public double getMinimumGoodAcceptanceLevel() {
        return MINIMUM_GOOD_ACCEPTANCE_LEVEL;
    }

    public double getMaximumGoodAcceptanceLevel() {
        return MAXIMUM_GOOD_ACCEPTANCE_LEVEL;
    }


    public String getPerformanceSuggestion() {
        //double d = OperatorUtils.optimizeWindowSize(getRawParameter(), parameter.getParameterValue(0) * 2.0, prob, targetProb);

        return getPerformanceSuggestion(
                getAcceptanceProbability(),
                getTargetAcceptanceProbability(),
                getRawParameter(),
                getAdaptableParameterName());
    }

    public AdaptationMode getMode() {
        return mode;
    }

    public static String getPerformanceSuggestion(double acceptanceProbability,
                                                  double targetAcceptanceProbability,
                                                  double adaptationParameterValue,
                                                  String adaptationParameterName) {

        dr.util.NumberFormatter formatter = new dr.util.NumberFormatter(5);

        double sf = OperatorUtils.optimizeWindowSize(adaptationParameterValue, acceptanceProbability, targetAcceptanceProbability);

        //double d = OperatorUtils.optimizeWindowSize(getRawParameter(), parameter.getParameterValue(0) * 2.0, prob, targetProb);


        if (acceptanceProbability < MINIMUM_ACCEPTANCE_LEVEL) {
            return "Try setting " + adaptationParameterName + " to about " + formatter.format(sf);
        } else if (acceptanceProbability > MAXIMUM_ACCEPTANCE_LEVEL) {
            return "Try setting " + adaptationParameterName + " to about " + formatter.format(sf);
        } else return "";
    }


}

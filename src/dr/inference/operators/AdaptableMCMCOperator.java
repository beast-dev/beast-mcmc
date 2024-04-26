/*
 * AdaptableMCMCOperator.java
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

/**
 * An MCMC operator that can be coerced to produce a target acceptance probability.
 *
 * @author Alexei Drummond
 * @author Andrew Rambaut
 * @version $Id: AdaptableMCMCOperator.java,v 1.3 2005/05/24 20:26:00 rambaut Exp $
 */
public interface AdaptableMCMCOperator extends MCMCOperator {

    String ADAPTABLE = "adaptable";
    String AUTO_OPTIMIZE = "autoOptimize";

    double DEFAULT_ADAPTATION_TARGET = 0.234;
    double MINIMUM_ACCEPTANCE_LEVEL = 0.1;
    double MAXIMUM_ACCEPTANCE_LEVEL = 0.4;
    double MINIMUM_GOOD_ACCEPTANCE_LEVEL = 0.2;
    double MAXIMUM_GOOD_ACCEPTANCE_LEVEL = 0.3;

    /**
     * An adaptable parameter must have a range from -infinity to +infinity with a preference for
     * small numbers.
     * <p/>
     * If operator acceptance is too high, BEAST increases the parameter; if operator acceptance is
     * too low, BEAST decreases the parameter.
     * <p/>
     * From MarkovChain.adaptAcceptanceProbability:
     * <p/>
     * new parameter = old parameter + 1/(1+N) * (current-step acceptance probability - target probability),
     * <p/>
     * where N is some function of the number of operator trials.
     *
     * @return an "adaptable" parameter
     */
    double getAdaptableParameter();

    double getTargetAcceptanceProbability();

    /**
     * Sets the adaptable parameter value.
     *
     * @param value the value to set the adaptable parameter to
     */
    void setAdaptableParameter(double value);

    /**
     * returns the number of times the setAdaptableParameter method has been called
     * @return the count
     */
    long getAdaptationCount();

    /**
     * sets the number of times the setAdaptableParameter method has been called (checkpointing)
     * @param count
     */
    void setAdaptationCount(long count);

    /**
     * @return the underlying tuning parameter value
     */
    double getRawParameter();

    String getAdaptableParameterName();

    double getMinimumAcceptanceLevel();

    double getMaximumAcceptanceLevel();

    double getMinimumGoodAcceptanceLevel();

    double getMaximumGoodAcceptanceLevel();

    /**
     * @return a short descriptive message of the performance of this operator.
     */
    String getPerformanceSuggestion();

    /**
     * @return the mode of this operator.
     */
    AdaptationMode getMode();

}

/*
 * AbstractAdaptableTreeOperator.java
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

package dr.evomodel.operators;

import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.operators.SubtreeLeapOperatorParser;
import dr.evomodelxml.operators.TipLeapOperatorParser;
import dr.inference.operators.AdaptableMCMCOperator;
import dr.inference.operators.AdaptationMode;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public abstract class AbstractAdaptableTreeOperator extends AbstractTreeOperator implements AdaptableMCMCOperator {

    public final AdaptationMode mode;
    private final double targetAcceptanceProbability;
    private long adaptationCount = 0;

    public AbstractAdaptableTreeOperator(AdaptationMode mode) {
        this(mode, DEFAULT_ADAPTATION_TARGET);
    }

    public AbstractAdaptableTreeOperator(AdaptationMode mode, double targetAcceptanceProbability) {
        this.mode = mode;
        if (System.getProperty("mcmc.adaptation_target") != null) {
            this.targetAcceptanceProbability = Double.parseDouble(System.getProperty("mcmc.adaptation_target"));
        } else {
            this.targetAcceptanceProbability = targetAcceptanceProbability;
        }
    }

    @Override
    public void setAdaptableParameter(double value) {
        adaptationCount ++;
        setAdaptableParameterValue(value);
    }

    @Override
    public double getAdaptableParameter() {
        return getAdaptableParameterValue();
    }

    @Override
    public long getAdaptationCount() {
        return adaptationCount;
    }

    @Override
    public void setAdaptationCount(long count) {
        adaptationCount = count;
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

    @Override
    public AdaptationMode getMode() {
        return mode;
    }

    @Override
    public double getMinimumAcceptanceLevel() {
        return 0.1;
    }

    @Override
    public double getMaximumAcceptanceLevel() {
        return 0.4;
    }

    @Override
    public double getMinimumGoodAcceptanceLevel() {
        return 0.20;
    }

    @Override
    public double getMaximumGoodAcceptanceLevel() {
        return 0.30;
    }

    @Override
    public final String getPerformanceSuggestion() {
        return null;
    }

    @Override
    public double getTargetAcceptanceProbability() {
        return targetAcceptanceProbability;
    }
}

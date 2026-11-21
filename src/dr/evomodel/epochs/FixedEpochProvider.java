/*
 * FixedEpochProvider.java
 *
 * Copyright © 2002-2026 the BEAST Development Team
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

package dr.evomodel.epochs;

import dr.inference.model.AbstractModel;
import dr.inference.model.Model;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;

/**
 * Epochs with transition times that are given when it is made and do not change.
 *
 * @author Andrew Rambaut
 */
public class FixedEpochProvider extends AbstractModel implements EpochProvider {

    public static final String FIXED_EPOCHS = "fixedEpochs";
    public static final String TRANSITION_TIMES = "transitionTimes";

    private final double[] transitionTimes;
    private final Parameter transitionTimeParameter;

    /**
     * @param transitionTimes the times at which one epoch gives way to the next: positive heights in
     *                        increasing order, in years. There is one more epoch
     *                        than transition times.
     */
    public FixedEpochProvider(double[] transitionTimes) {
        super(FIXED_EPOCHS);

        for (int i = 0; i < transitionTimes.length; i++) {
            final double previous = i == 0 ? 0.0 : transitionTimes[i - 1];
            // written so that a NaN fails too
            if (!(transitionTimes[i] > previous) || Double.isInfinite(transitionTimes[i])) {
                throw new IllegalArgumentException("The transition times of the epochs must be positive, finite " +
                        "and in increasing order (transition time " + (i + 1) + " is " + transitionTimes[i] +
                        (i == 0 ? ")" : ", after " + previous + ")"));
            }
        }

        this.transitionTimes = transitionTimes.clone();
        this.transitionTimeParameter = new Parameter.Default(TRANSITION_TIMES, this.transitionTimes.clone());
    }

    /**
     * Epochs of the same duration, apart from the last, which has no end.
     *
     * @param epochCount the number of epochs, counting the last
     * @param duration   the duration of each epoch, in years
     * @return epochs with transition times at duration, 2 x duration ... (epochCount - 1) x duration
     */
    public static FixedEpochProvider equalDurations(int epochCount, double duration) {
        return new FixedEpochProvider(equalDurationTransitionTimes(epochCount, duration));
    }

    public static double[] equalDurationTransitionTimes(int epochCount, double duration) {
        if (epochCount < 1) {
            throw new IllegalArgumentException("There must be at least one epoch");
        }
        if (!(duration > 0.0) || Double.isInfinite(duration)) {
            throw new IllegalArgumentException("The duration of the epochs must be positive and finite");
        }
        double[] transitionTimes = new double[epochCount - 1];
        for (int i = 0; i < transitionTimes.length; i++) {
            transitionTimes[i] = (i + 1) * duration;
        }
        return transitionTimes;
    }

    @Override
    public int getEpochCount() {
        return transitionTimes.length + 1;
    }

    @Override
    public double[] getTransitionTimes() {
        return transitionTimes;
    }

    @Override
    public Parameter getTransitionTimeParameter() {
        return transitionTimeParameter;
    }

    // **************************************************************
    // Model IMPLEMENTATION: nothing changes, so there is nothing to do
    // **************************************************************

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
    }

    @Override
    protected void storeState() {
    }

    @Override
    protected void restoreState() {
    }

    @Override
    protected void acceptState() {
    }
}

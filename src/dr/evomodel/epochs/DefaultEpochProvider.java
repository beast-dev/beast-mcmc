/*
 * DefaultEpochProvider.java
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
 * Epochs with transition times held in a parameter. They are given when it is made, as an array or as the
 * values of a parameter, and change only if something changes the parameter.
 *
 * Given its own parameter, the provider does not change. Given one from outside, operators can be applied
 * to it; the provider tells its listeners when it changes but does not check that the times stay positive
 * and in order, which is for the bounds and operators on the parameter to see to.
 *
 * @author Andrew Rambaut
 */
public class DefaultEpochProvider extends AbstractModel implements EpochProvider {

    public static final String EPOCHS = "epochs";
    public static final String TRANSITION_TIMES = "transitionTimes";

    private final Parameter transitionTimeParameter;

    // the values of the parameter, copied when they are asked for after a change
    private double[] transitionTimes;
    private boolean transitionTimesKnown = false;

    /**
     * @param transitionTimes the times at which one epoch gives way to the next: positive heights in
     *                        increasing order, in years. There is one more epoch
     *                        than transition times.
     */
    public DefaultEpochProvider(double[] transitionTimes) {
        this(new Parameter.Default(TRANSITION_TIMES, transitionTimes.clone()));
    }

    /**
     * @param transitionTimeParameter a parameter whose values are the transition times
     */
    public DefaultEpochProvider(Parameter transitionTimeParameter) {
        super(EPOCHS);

        checkTransitionTimes(transitionTimeParameter.getParameterValues());

        this.transitionTimeParameter = transitionTimeParameter;
        addVariable(transitionTimeParameter);
    }

    /**
     * @param transitionTimeParameter a parameter that is given the transition times: its dimension is set to
     *                                their number and its values to them
     * @param transitionTimes         the transition times
     */
    public DefaultEpochProvider(Parameter transitionTimeParameter, double[] transitionTimes) {
        this(setValues(transitionTimeParameter, transitionTimes));
    }

    private static Parameter setValues(Parameter parameter, double[] values) {
        checkTransitionTimes(values);
        if (parameter.getDimension() != values.length) {
            parameter.setDimension(values.length);
        }
        for (int i = 0; i < values.length; i++) {
            parameter.setParameterValueQuietly(i, values[i]);
        }
        return parameter;
    }

    private static void checkTransitionTimes(double[] transitionTimes) {
        for (int i = 0; i < transitionTimes.length; i++) {
            final double previous = i == 0 ? 0.0 : transitionTimes[i - 1];
            // written so that a NaN fails too
            if (!(transitionTimes[i] > previous) || Double.isInfinite(transitionTimes[i])) {
                throw new IllegalArgumentException("The transition times of the epochs must be positive, finite " +
                        "and in increasing order (transition time " + (i + 1) + " is " + transitionTimes[i] +
                        (i == 0 ? ")" : ", after " + previous + ")"));
            }
        }
    }

    /**
     * Epochs of the same duration, apart from the last, which has no end.
     *
     * @param epochCount the number of epochs, counting the last
     * @param duration   the duration of each epoch, in years
     * @return epochs with transition times at duration, 2 x duration ... (epochCount - 1) x duration
     */
    public static DefaultEpochProvider equalDurations(int epochCount, double duration) {
        return new DefaultEpochProvider(equalDurationTransitionTimes(epochCount, duration));
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
        return transitionTimeParameter.getDimension() + 1;
    }

    @Override
    public double[] getTransitionTimes() {
        if (!transitionTimesKnown) {
            transitionTimes = transitionTimeParameter.getParameterValues();
            transitionTimesKnown = true;
        }
        return transitionTimes;
    }

    @Override
    public Parameter getTransitionTimeParameter() {
        return transitionTimeParameter;
    }

    // **************************************************************
    // Model IMPLEMENTATION
    // **************************************************************

    @Override
    protected void handleModelChangedEvent(Model model, Object object, int index) {
        // there are no models
    }

    @Override
    protected void handleVariableChangedEvent(Variable variable, int index, Variable.ChangeType type) {
        // the transition times have changed (AbstractModel tells the listeners)
        transitionTimesKnown = false;
    }

    @Override
    protected void storeState() {
        // the parameter stores its own values
    }

    @Override
    protected void restoreState() {
        transitionTimesKnown = false;
    }

    @Override
    protected void acceptState() {
    }
}

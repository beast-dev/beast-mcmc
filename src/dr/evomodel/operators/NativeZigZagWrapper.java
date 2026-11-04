/*
 * NativeZigZagWrapper.java
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

import dr.inference.hmc.PrecisionColumnProvider;
import dr.inference.operators.hmc.MinimumTravelInformation;
import dr.inference.operators.hmc.MinimumTravelInformationBinary;

/**
 * @author Marc A. Suchard
 */
public class NativeZigZagWrapper {

    private final int instanceNumber;

    public NativeZigZagWrapper(int dimension,
                               NativeZigZagOptions options,
                               double[] mask,
                               double[] observed,
                               double[] parameterSign,
                               double[] lb,
                               double[] ub) {
        this.instanceNumber = NativeZigZag.INSTANCE.createInstance(dimension, options, mask, observed, parameterSign,
                lb, ub);
    }

    public void operate(PrecisionColumnProvider columnProvider,
                        double[] position,
                        double[] velocity,
                        double[] action,
                        double[] gradient,
                        double[] moment,
                        double time) {
        NativeZigZag.INSTANCE.operate(instanceNumber, columnProvider, position, velocity, action, gradient, moment,
                time);
    }

    public MinimumTravelInformationBinary getNextReversibleEvent(double[] position,
                                                                 double[] velocity,
                                                                 double[] action,
                                                                 double[] gradient,
                                                                 double[] momentum) {
        return NativeZigZag.INSTANCE.getNextEvent(instanceNumber, position, velocity, action, gradient, momentum);
    }

    public MinimumTravelInformation getNextIrreversibleEvent(double[] position,
                                                             double[] velocity,
                                                             double[] action,
                                                             double[] gradient) {
        throw new RuntimeException("not implemented yet");
    }

    public void updateReversibleDynamics(double[] position,
                                         double[] velocity,
                                         double[] action,
                                         double[] gradient,
                                         double[] momentum,
                                         double[] column,
                                         double eventTime, int eventIndex, int eventType) {
        NativeZigZag.INSTANCE.updateDynamics(instanceNumber, position, velocity, action, gradient, momentum,
                column, eventTime, eventIndex, eventType);
    }

    public void updateIrreversibleDynamics(double[] position,
                                           double[] velocity,
                                           double[] action,
                                           double[] gradient,
                                           double[] column,
                                           double eventTime, int eventIndex, int eventType) {
        NativeZigZag.INSTANCE.updateDynamics(instanceNumber, position, velocity, action, gradient, null,
                column, eventTime, eventIndex, eventType);
    }

    @SuppressWarnings("unused")
    public int enterCriticalRegion(
            double[] position,
            double[] velocity,
            double[] action,
            double[] gradient,
            double[] momentum) {
        return NativeZigZag.INSTANCE.enterCriticalRegion(instanceNumber, position, velocity, action, gradient,
                momentum);
    }

    @SuppressWarnings("unused")
    public int exitCriticalRegion() {
        return NativeZigZag.INSTANCE.exitCriticalRegion(instanceNumber);
    }

    @SuppressWarnings("unused")
    public boolean inCriticalRegion() {
        return NativeZigZag.INSTANCE.inCriticalRegion(instanceNumber);
    }

    @SuppressWarnings("unused")
    public MinimumTravelInformation getNextEventInCriticalRegion() {
        return NativeZigZag.INSTANCE.getNextEventInCriticalRegion(instanceNumber);
    }

    @SuppressWarnings("unused")
    public void innerBounce(double[] position,
                            double[] velocity,
                            double[] action,
                            double[] gradient,
                            double[] momentum,
                            double eventTime, int eventIndex, int eventType) {
        NativeZigZag.INSTANCE.innerBounce(instanceNumber, position, velocity, action, gradient, momentum,
                eventTime, eventIndex, eventType);
    }

    @SuppressWarnings("unused")
    public void innerBounceCriticalRegion(double eventTime, int eventIndex, int eventType) {
        NativeZigZag.INSTANCE.innerBounceCriticalRegion(instanceNumber, eventTime, eventIndex, eventType);
    }
}

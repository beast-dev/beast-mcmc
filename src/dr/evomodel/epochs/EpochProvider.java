/*
 * EpochProvider.java
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

package dr.evomodel.epochs;

import dr.inference.model.Model;
import dr.inference.model.Parameter;

/**
 * A division of time into epochs, for the models that are piecewise in time (epochs of a birth-death
 * model, of a branch rate model, the pieces of a non-parametric coalescent model).
 *
 * Time is measured as height: back from zero, the time of the most recent tip, in years, the units of time of
 * the tree. The epochs are numbered from 0, the most recent, backwards. With N epochs there are N - 1
 * transition times: epoch 0 runs from height zero to the first transition time, epoch i from transition
 * time i - 1 to transition time i, and the last epoch from the last transition time back without end. A
 * height that is exactly a transition time is in the older of the two epochs.
 *
 * It is a Model so that an implementation whose transition times can change tells its listeners.
 *
 * @author Andrew Rambaut
 */
public interface EpochProvider extends Model {

    /**
     * @return the number of epochs, one more than the number of transition times
     */
    int getEpochCount();

    /**
     * @return the transition times, in increasing order. This may be the provider's own array and must
     * not be changed.
     */
    double[] getTransitionTimes();

    /**
     * @return the transition times as a parameter, for the models that take them as one
     */
    Parameter getTransitionTimeParameter();

    /**
     * @param height a time
     * @return the index of the epoch the time is in
     */
    default int getEpochIndex(double height) {
        final double[] transitionTimes = getTransitionTimes();
        // the number of transition times at or before the height, by binary search
        int low = 0;
        int high = transitionTimes.length;
        while (low < high) {
            final int mid = (low + high) >>> 1;
            if (transitionTimes[mid] <= height) {
                low = mid + 1;
            } else {
                high = mid;
            }
        }
        return low;
    }
}

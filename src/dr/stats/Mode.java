/*
 * Mode.java
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

package dr.stats;


/**
 * simple statistics related to mode.
 * the input is a frequency counter,
 * keys are unique values, Integer values of map are counts of that key
 *
 * @author Walter Xie
 */
public class Mode <T> {
    protected T mode;
    protected int freqOfMode = 0;

    public Mode(FrequencyCounter<T> frequencyCounter) {
        calculateMode(frequencyCounter);
    }

    private void calculateMode(FrequencyCounter<T> frequencyCounter) {
        for (T key : frequencyCounter.uniqueValues()) {
            if (freqOfMode < frequencyCounter.getCount(key)) {
                freqOfMode = frequencyCounter.getCount(key);
                mode = key;
            }
        }
    }

    public T getMode() {
        return mode;
    }

    public int getFrequencyOfMode() {
        return freqOfMode;
    }
}

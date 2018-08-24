/*
 * FrequencyCounter.java
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


import dr.util.Pair;

import java.util.*;

/**
 * A utility class to store the frequency of a set of discrete values
 * and provide some basic statistics on them (mode, credible sets etc).
 *
 * @author Walter Xie
 * @author Andrew Rambaut
 */
public final class FrequencyCounter<T extends Comparable> {

    private final List<T> uniqueValues;
    private final Map<T, Integer> frequencies;
    private final Set<T> credibleSet = new LinkedHashSet<T>();
    private final Set<T> incredibleSet = new LinkedHashSet<T>();

    private final int total;
    private final int min;
    private final int max;
    private final T mode;

    public FrequencyCounter(List<T> values) {
        this(values, 0.0);
    }
    
    public FrequencyCounter(List<T> values, double probabilityThreshold) {

        // use a treemap so the keys are sorted by their value
        frequencies = new TreeMap<T, Integer>();

        for (T value : values) {
            if (frequencies.containsKey(value)) {
                int i = frequencies.get(value) + 1;
                frequencies.put(value, i);
            } else {
                frequencies.put(value, 1);
            }
        }

        this.uniqueValues = new ArrayList<T>(frequencies.keySet());

        total = calculateTotalFrequency();
        int[] minMax = calculateMinMaxFrequency();
        min = minMax[0];
        max = minMax[1];
        mode = calculateMode();

        if (probabilityThreshold > 0.0) {
            calculateCredibleSet(probabilityThreshold);
        }
    }

    public Map<T, Integer> getFrequencies() {
        return frequencies;
    }

    public int getSize() {
        return frequencies.size();
    }

    /**
     * sort counter by counts to calculate correct credibility set
     */
    public List<Integer> getOrderByFrequency() {
        List<Pair<T, Integer>> values = new ArrayList<Pair<T, Integer>>();
        int i = 0;
        for (T value : uniqueValues) {
            values.add(new Pair<T, Integer>(value, i));
            i++;
        }
        Collections.sort(values, new Comparator<Pair<T, Integer>>() {
            public int compare(Pair<T, Integer> value1, Pair<T, Integer> value2) {
                return getFrequency(value2.fst) - getFrequency(value1.fst);
            }
        });
        List<Integer> order = new ArrayList<Integer>();
        for (Pair<T, Integer> value : values) {
            order.add(value.snd);
        }
        return order;
    }

    /**
     * the unique values in a frequency counter,
     * which are also the keys of the map
     *
     * @return
     */
    public List<T> getUniqueValues() {
        List<T> values = new ArrayList<T>(uniqueValues);
        Collections.sort(values);
        return values;
    }


    public int getMinFrequency() {
        return min;
    }

    public int getMaxFrequency() {
        return max;
    }

    public T getMode() {
        return mode;
    }


    /**
     * get the count from counter given a value
     *
     * @param value
     * @return
     */
    public int getFrequency(T value) {
        if (frequencies.containsKey(value)) {
            return frequencies.get(value);
        } else {
            return 0;
        }
    }

    /**
     * the total counts of this frequency counter
     *
     * @return
     */
    public int getTotalFrequency() {
        return total;
    }

    public double getProbability(T value) {
        return (double) getFrequency(value) / (double)total;
    }

    /**
     * rescale the frequency to make maximum count equal to 1.
     * <code>count / max_count</code>.
     * @param value
     * @return
     */
    public double getProportionalFrequency(T value) {
        return (double) getFrequency(value) / (double)max;
    }

    public Set<T> getCredibleSet() {
        return credibleSet;
    }

    public Set<T> getIncredibleSet() {
        return incredibleSet;
    }

    private int calculateTotalFrequency() {
        int total = 0;
        for (int value : frequencies.values()) {
            total += value;
        }
        return total;
    }

    private int[] calculateMinMaxFrequency() {
        int min = Collections.min(frequencies.values());
        int max = Collections.max(frequencies.values());
        return new int[]{min, max};
    }

    private T calculateMode() {
        int maxFreq = 0;
        T mode = null;
        for (T value : getUniqueValues()) {
            if (getFrequency(value) > maxFreq) {
                maxFreq = getFrequency(value);
                mode = value;
            }
        }
        return mode;
    }

    private void calculateCredibleSet(double probabilityThreshold) {
        double totalProbability = 0;
        for (int index : getOrderByFrequency()) {
            T value = uniqueValues.get(index);
            double probability = getProbability(value);
            // include the last one to make totPr >= probability
            if (totalProbability < probabilityThreshold) {
                credibleSet.add(value);
            } else {
                incredibleSet.add(value);
            }
            totalProbability += probability;
        }
    }

}

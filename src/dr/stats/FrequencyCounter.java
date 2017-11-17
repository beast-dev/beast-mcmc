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


import java.util.*;

/**
 * frequency counter
 *
 * @author Walter Xie
 * @author Andrew Rambaut
 */
public final class FrequencyCounter<T> {

    private Map<T, Integer> frequencies;
    private List<T> sortedValues = new ArrayList<T>();
    private Set<T> credibleSet = new TreeSet<T>();
    private Set<T> incredibleSet = new TreeSet<T>();

    private int total;
    private int min;
    private int max;
    private T mode;

    public FrequencyCounter(List<T> values) {
        this(values, 0.05);
    }
    
    public FrequencyCounter(List<T> values, double probabilityThreshold) {
        // http://stackoverflow.com/questions/12998568/hashmap-vs-linkedhashmap-performance-in-iteration-over-values
        frequencies = new LinkedHashMap<T, Integer>();

        for (T value : values) {
            if (frequencies.containsKey(value)) {
                int i = frequencies.get(value) + 1;
                frequencies.put(value, i);
            } else {
                frequencies.put(value, 1);
            }
        }

        total = calculateTotalFrequency();
        int[] minMax = calculateMinMaxFrequency();
        min = minMax[0];
        max = minMax[1];
        mode = calculateMode();
        calculateCredibleSet(probabilityThreshold);
    }

    public Map<T, Integer> getFrequencies() {
        return frequencies;
    }

    public int getCounterSize() {
        return frequencies.size();
    }

//    /**
//     * sort counter by counts to calculate correct credibility set
//     */
//    public void sortCounterByCounts() {
//        frequencies = Utils.sortByValue(frequencies);
//    }
//
    /**
     * the unique values in a frequency counter,
     * which are also the keys of the map
     *
     * @return
     */
    public Set<T> getUniqueValues() {
        return new TreeSet<T>(frequencies.keySet());
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
        return frequencies.get(value);
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
        int min = Integer.MAX_VALUE;
        int max = Integer.MIN_VALUE;
        for (int value : frequencies.values()) {
            if (min > value)
                min = value;
            if (max < value)
                max = value;
        }
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
        for (T value : getUniqueValues()) {
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

    public static class Utils {
        // http://stackoverflow.com/questions/109383/sort-a-mapkey-value-by-values-java
        public static <K, V extends Comparable<? super V>> Map<K, V> sortByValue(Map<K, V> map) {
            List<Map.Entry<K, V>> list = new LinkedList<Map.Entry<K, V>>(map.entrySet());
            Collections.sort(list, new Comparator<Map.Entry<K, V>>() {
                @Override
                public int compare(Map.Entry<K, V> o1, Map.Entry<K, V> o2) {
                    return (o1.getValue()).compareTo(o2.getValue()) * -1;
                }
            });

            Map<K, V> result = new LinkedHashMap<K, V>();
            for (Map.Entry<K, V> entry : list) {
                result.put(entry.getKey(), entry.getValue());
            }
            return result;
        }
    }

}

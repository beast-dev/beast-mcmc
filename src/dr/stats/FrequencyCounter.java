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
 */
public class FrequencyCounter<T> {

    private int MAX_COUNTER_SIZE = 100;
    protected Map<T, Integer> frequencies;
    protected boolean sortedByCounts = false;

    protected int total;
    protected int min;
    protected int max;

    public FrequencyCounter(List<T> values, boolean sortedByCounts) {
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

        // limit the counter size, avoid expensive computation
//        if (frequencies.size() > MAX_COUNTER_SIZE)
//            throw new IllegalArgumentException("Fail to create frequency counter: " +
//                    "number of unique values must <=" + MAX_COUNTER_SIZE + " !");


        if (sortedByCounts && frequencies.size() > 0)
            sortCounterByCounts();

        // store {total min max} counts
        total = calculateTotalCount();
        int[] minMax = calculateMinMaxCount();
        min = minMax[0];
        max = minMax[1];
    }

    public Map<T, Integer> getFrequencies() {
        return frequencies;
    }

    public int getCounterSize() {
        return frequencies.size();
    }

    /**
     * sort counter by counts to calculate correct credibility set
     */
    public void sortCounterByCounts() {
        frequencies = Utils.sortByValue(frequencies);
        sortedByCounts = true;
    }

    public boolean isSortedByCounts() {
        return sortedByCounts;
    }

    /**
     * the unique values in a frequency counter,
     * which are also the keys of the map
     *
     * @param sort
     * @return
     */
    public Set<T> uniqueValues(boolean sort) {
        if (sort)
            return new TreeSet<T>(frequencies.keySet());
        else
            return frequencies.keySet();
    }

    /**
     * default to sort unique values (keys).
     *
     * @return
     */
    public Set<T> uniqueValues() {
        return uniqueValues(true);
    }

    public String uniqueValuesToString() {
        return Utils.setToString(uniqueValues());
    }


    /**
     * sort the key of frequency counter, and return the index of a given key,
     * if not exist, return -1.
     *
     * @param key
     * @return
     */
    public int getKeyIndex(T key) {
        int i = -1;
        for (T v : uniqueValues()) {
            i++;
            if (v.equals(key)) return i;
        }
        return i;
    }

    /**
     * get the count from counter given a key
     *
     * @param key
     * @return
     */
    public int getCount(T key) {
        return frequencies.get(key);
    }

    /**
     * the total counts of this frequency counter
     *
     * @return
     */
    public int calculateTotalCount() {
        int tot = 0;
        for (Map.Entry<T, Integer> entry : frequencies.entrySet()) {
            Integer count = entry.getValue();
            tot += count;
        }
        return tot;
    }

    public int getTotalCount() {
        if (total <= 0)
            total= calculateTotalCount();
        return total;
    }

    public double getProbability(T key) {
        return (double) getCount(key) / (double) getTotalCount();
    }

    /**
     * rescale the frequency to make maximum count equal to 1.
     * <code>count / max_count</code>.
     * @param key
     * @return
     */
    public double getFreqScaledMaxTo1(T key) {
        return (double) getCount(key) / (double) getMaxCount();
    }

    /**
     * the min and max count in the frequency counter
     *
     * @return <code>int[]</code>, 1st is min, 2nd is max.
     */
    public int[] calculateMinMaxCount() {
        int min = Integer.MAX_VALUE;
        int max = Integer.MIN_VALUE;
        for (Map.Entry<T, Integer> entry : frequencies.entrySet()) {
            Integer count = entry.getValue();
            if (min > count)
                min = count;
            if (max < count)
                max = count;
        }
        return new int[]{min, max};
    }

    public int getMaxCount() {
        if (max <= 0) {
            int[] minMax = calculateMinMaxCount();
            min = minMax[0];
            max = minMax[1];
        }
        return max;
    }

    /**
     * mode calculated from this frequency counter
     *
     * @return
     */
    public Mode<T> getModeStats() {
        return new Mode<T>(this);
    }

    /**
     * Include credible and incredible set calculated from this frequency counter
     *
     * @param probability
     * @return
     */
    public CredibleSetAnalysis<T> getCredibleSetAnalysis(double probability) {
        return new CredibleSetAnalysis<T>(this, probability);
    }

    public static class Utils {
        public static <T> String setToString(Set<T> aSet) {
            String line = "{";
            for (T value : aSet) {
                line = line + value + ", ";
            }
            if (line.endsWith(", ")) {
                line = line.substring(0, line.lastIndexOf(", ")) + "}";
            } else {
                line = "{}";
            }
            return line;
        }


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

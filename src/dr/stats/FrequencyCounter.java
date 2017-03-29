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
    protected Map<T, Integer> frequencyCounter;

    public FrequencyCounter(List<T> values) {
        frequencyCounter = new HashMap<T, Integer>();
        for (T value : values) {
            if (frequencyCounter.containsKey(value)) {
                int i = frequencyCounter.get(value) + 1;
                frequencyCounter.put(value, i);
            } else {
                frequencyCounter.put(value, 1);
            }
        }

        // limit the counter size, avoid expensive computation
//        if (frequencyCounter.size() > MAX_COUNTER_SIZE)
//            throw new IllegalArgumentException("Fail to create frequency counter: " +
//                    "number of unique values must <=" + MAX_COUNTER_SIZE + " !");
    }

    public Map<T, Integer> getFrequencyCounter() {
        return frequencyCounter;
    }

    public int getCounterSize() {
        return frequencyCounter.size();
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
            return new TreeSet<T>(frequencyCounter.keySet());
        else
            return frequencyCounter.keySet();
    }

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
        return frequencyCounter.get(key);
    }

    /**
     * the total counts of this frequency counter
     *
     * @return
     */
    public int getTotalCount() {
        int tot = 0;
        for (Map.Entry<T, Integer> entry : frequencyCounter.entrySet()) {
            Integer count = entry.getValue();
            tot += count;
        }
        return tot;
    }

    /**
     * the min and max count in the frequency counter
     *
     * @return <code>int[]</code>, 1st is min, 2nd is max.
     */
    public int[] getMinMaxCount() {
        int min = Integer.MAX_VALUE;
        int max = Integer.MIN_VALUE;
        for (Map.Entry<T, Integer> entry : frequencyCounter.entrySet()) {
            Integer count = entry.getValue();
            if (min > count)
                min = count;
            if (max < count)
                max = count;
        }
        return new int[]{min, max};
    }

    /**
     * mode calculated from this frequency counter
     *
     * @return
     */
    public Mode<T> getMode() {
        return new Mode<T>(this);
    }

    /**
     * credible and incredible set calculated from this frequency counter
     *
     * @param probability
     * @return
     */
    public CredibleSet<T> getCredibleSet(double probability) {
        return new CredibleSet<T>(this, probability);
    }

    static class Utils {
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
    }

}

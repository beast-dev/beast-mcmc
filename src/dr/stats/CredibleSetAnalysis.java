/*
 * CredibleSet.java
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


import java.util.Set;
import java.util.TreeSet;

/**
 * credible set
 *
 * @author Walter Xie
 */
public class CredibleSetAnalysis<T> {
    protected double probability;
    protected Set<T> credibleSet = new TreeSet<T>();

    protected Set<T> incredibleSet = new TreeSet<T>();

    public CredibleSetAnalysis(FrequencyCounter<T> frequencyCounter, double probability) {
        calculateCredibleSet(frequencyCounter, probability);
    }

    private void calculateCredibleSet(FrequencyCounter<T> frequencyCounter, double probability) {
        this.probability = probability;

        if (! frequencyCounter.isSortedByCounts())
            frequencyCounter.sortCounterByCounts();

        double totPr = 0;
        for (T key : frequencyCounter.uniqueValues(false)) { // use the original ordering of unique values
            // frequency / total size
            double prob = frequencyCounter.getProbability(key);
            // include the last one to make totPr >= probability
            if (totPr < probability) {
                credibleSet.add(key);
            } else {
                incredibleSet.add(key);
            }
            totPr += prob;
        }

        // totPr = 1.0000000000000002
        totPr = (double) Math.round(totPr * 1000000.0) / 1000000.0;

        if (totPr != 1)
            throw new IllegalArgumentException("Probability NOT sums to 1 in credibility set analysis ! totPr = " + totPr);
    }

    public Set<T> getCredibleSet() {
        return credibleSet;
    }

    public Set<T> getIncredibleSet() {
        return incredibleSet;
    }

    public String toStringCredibleSet() {
        return FrequencyCounter.Utils.setToString(credibleSet);
    }

    public String toStringIncredibleSet() {
        return FrequencyCounter.Utils.setToString(incredibleSet);
    }
}

/*
 * PatternList.java
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

package dr.evolution.alignment;

import dr.evolution.datatype.DataType;
import dr.evolution.util.TaxonList;
import dr.util.Identifiable;

/**
 * interface for any list of patterns with weights.
 *
 * @author Andrew Rambaut
 * @version $Id: PatternList.java,v 1.12 2005/05/24 20:25:55 rambaut Exp $
 */
public interface PatternList extends TaxonList, Identifiable {
    /**
     * @return number of patterns
     */
    int getPatternCount();

    /**
     * @return number of states for this PatternList
     */
    int getStateCount();

    /**
     * Gets the length of the pattern strings which will usually be the
     * same as the number of taxa
     *
     * @return the length of patterns
     */
    int getPatternLength();

    /**
     * Gets the pattern as an array of state numbers (one per sequence)
     *
     * @param patternIndex the index of the pattern to return
     * @return the site pattern at patternIndex
     */
    int[] getPattern(int patternIndex);

    /**
     * Gets the pattern as an array of state frequency vectors (one per sequence)
     *
     * @param patternIndex the index of the pattern to return
     * @return the site pattern at patternIndex
     */
    double[][] getUncertainPattern(int patternIndex);

    /**
     * @param taxonIndex   the taxon
     * @param patternIndex the pattern
     * @return state at (taxonIndex, patternIndex)
     */
    int getPatternState(int taxonIndex, int patternIndex);

    /**
     * @param taxonIndex   the taxon
     * @param patternIndex the pattern
     * @return state frequency vector at (taxonIndex, patternIndex)
     */
    double[] getUncertainPatternState(int taxonIndex, int patternIndex);

    /**
     * Gets the weight of a site pattern
     *
     * @param patternIndex the pattern
     * @return the weight of the specified pattern
     */
    double getPatternWeight(int patternIndex);

    /**
     * @return the array of pattern weights
     */
    double[] getPatternWeights();

    /**
     * @return the DataType of this PatternList
     */
    DataType getDataType();

    /**
     * @return the frequency of each state
     */
    double[] getStateFrequencies();

    /**
     * Are the patterns only the unique ones (i.e., compressed)?
     * @return are unique?
     */
    boolean areUnique();

    /**
     * Do the patterns contain any uncertain states?
     * @return are uncertain?
     */
    boolean areUncertain();

    /**
     * Helper routines for pattern lists.
     */
    public static class Utils {
        /**
         * Returns a double array containing the empirically estimated frequencies
         * for the states of patternList. This currently calls the version that maps
         * to PAUP's estimates.
         *
         * @param patternList the pattern list to calculate the empirical state
         *                    frequencies from
         * @return the empirical state frequencies of the given pattern list
         */
        public static double[] empiricalStateFrequencies(PatternList patternList) {
            return empiricalStateFrequenciesPAUP(patternList);
        }

        /**
         * Returns a double array containing the empirically estimated frequencies
         * for the states of patternList. This version of the routine should match
         * the values produced by PAUP.
         *
         * @param patternList the pattern list to calculate the empirical state
         *                    frequencies from
         * @return the empirical state frequencies of the given pattern list
         */
        public static double[] empiricalStateFrequenciesPAUP(PatternList patternList) {
            int i, j, k;
            double total, sum, x, w, difference;

            DataType dataType = patternList.getDataType();

            int stateCount = patternList.getStateCount();
            int patternLength = patternList.getPatternLength();
            int patternCount = patternList.getPatternCount();

            double[] freqs = equalStateFrequencies(patternList);

            double[] tempFreq = new double[stateCount];
            int[] pattern;
            boolean[] state;

            int count = 0;
            do {
                for (i = 0; i < stateCount; i++)
                    tempFreq[i] = 0.0;

                total = 0.0;
                for (i = 0; i < patternCount; i++) {
                    pattern = patternList.getPattern(i);
                    w = patternList.getPatternWeight(i);

                    for (k = 0; k < patternLength; k++) {
                        state = dataType.getStateSet(pattern[k]);

                        sum = 0.0;
                        for (j = 0; j < stateCount; j++)
                            if (state[j])
                                sum += freqs[j];

                        for (j = 0; j < stateCount; j++) {
                            if (state[j]) {
                                x = (freqs[j] * w) / sum;
                                tempFreq[j] += x;
                                total += x;
                            }
                        }
                    }

                }

                difference = 0.0;
                for (i = 0; i < stateCount; i++) {
                    difference += Math.abs((tempFreq[i] / total) - freqs[i]);
                    freqs[i] = tempFreq[i] / total;
                }
                count ++;
            } while (difference > 1E-8 && count < 1000);

            return freqs;
        }

        /**
         * Returns a double array containing the empirically estimated frequencies
         * for the states of patternList. This version of the routine should match
         * the values produced by MrBayes.
         *
         * @param patternList the pattern list to calculate the empirical state
         *                    frequencies from
         * @return the empirical state frequencies of the given pattern list
         */
        public static double[] empiricalStateFrequenciesMrBayes(PatternList patternList) {

            DataType dataType = patternList.getDataType();

            int stateCount = patternList.getStateCount();
            int patternLength = patternList.getPatternLength();
            int patternCount = patternList.getPatternCount();

            double[] freqs = equalStateFrequencies(patternList);

            double sumTotal = 0.0;

            double[] sumFreq = new double[stateCount];

            for (int i = 0; i < patternCount; i++) {
                int[] pattern = patternList.getPattern(i);
                double w = patternList.getPatternWeight(i);

                for (int k = 0; k < patternLength; k++) {
                    boolean[] state = dataType.getStateSet(pattern[k]);

                    double sum = 0.0;
                    for (int j = 0; j < stateCount; j++) {
                        if (state[j]) {
                            sum += freqs[j];
                        }
                    }

                    for (int j = 0; j < stateCount; j++) {
                        if (state[j]) {
                            double x = (freqs[j] * w) / sum;
                            sumFreq[j] += x;
                            sumTotal += x;
                        }
                    }
                }

            }

            for (int i = 0; i < stateCount; i++) {
                freqs[i] = sumFreq[i] / sumTotal;
            }

            return freqs;
        }

        /**
         * Returns a double array containing the equal frequencies
         * for the states of patternList.
         *
         * @param patternList the pattern list
         * @return return equal state frequencies based on the data type of
         *         the patternlist
         */
        public static double[] equalStateFrequencies(PatternList patternList) {
            int i, n = patternList.getStateCount();
            double[] freqs = new double[n];
            double f = 1.0 / n;

            for (i = 0; i < n; i++)
                freqs[i] = f;

            return freqs;
        }
    }
}

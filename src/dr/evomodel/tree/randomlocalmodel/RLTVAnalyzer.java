/*
 * RLTVAnalyzer.java
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

package dr.evomodel.tree.randomlocalmodel;

import dr.stats.DiscreteStatistics;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.*;

/**
 * @author Alexei Drummond
 */
public class RLTVAnalyzer {

    static TreeMap<Integer, List<Double>> nodeRates = new TreeMap<Integer, List<Double>>();
    static Map<String, Integer> comboCounts = new HashMap<String, Integer>();
    static Map<String, Integer> pairCounts = new HashMap<String, Integer>();

    static NumberFormat format = NumberFormat.getNumberInstance(Locale.ENGLISH);

    public static void main(String[] args) throws IOException {

        int total = 0;

        double priorChange = Double.parseDouble(args[1]);

        format.setMaximumFractionDigits(5);

        BufferedReader reader = new BufferedReader(new FileReader(args[0]));

        // read header line
        String line = reader.readLine();
        while (line.toLowerCase().startsWith("state")) {
            line = reader.readLine();
        }

        while (line != null) {
            StringTokenizer tokens = new StringTokenizer(line);

            // state
            tokens.nextToken();

            List<Integer> changes = new ArrayList<Integer>();
            while (tokens.hasMoreTokens()) {
                int nodeNumber = Integer.parseInt(tokens.nextToken());
                double rate = Double.parseDouble(tokens.nextToken());

                changes.add(nodeNumber);

                putRate(nodeNumber, rate);
            }
            putCombo(changes);
            total += 1;
            line = reader.readLine();
        }

        List<RateCombo> rateCombos = new ArrayList<RateCombo>();
        for (Map.Entry<String, Integer> stringIntegerEntry1 : comboCounts.entrySet()) {

            int count = stringIntegerEntry1.getValue();
            if (count > 1) rateCombos.add(new RateCombo(stringIntegerEntry1.getKey(), count));
        }

        Collections.sort(rateCombos, rateComboComparator);

        System.out.println("combo\tPr");
        for (RateCombo rateCombo : rateCombos) {
            double posterior = (double) rateCombo.count / (double) total;

            String[] nodes = rateCombo.combo.split("-");

            System.out.println(rateCombo.combo + "\t" + format.format(posterior));
        }

        List<RateCombo> pairRateCombos = new ArrayList<RateCombo>();
        for (Map.Entry<String, Integer> stringIntegerEntry : pairCounts.entrySet()) {
            int count = stringIntegerEntry.getValue();
            if (count > 1) pairRateCombos.add(new RateCombo(stringIntegerEntry.getKey(), count));
        }

        Collections.sort(pairRateCombos, rateComboComparator);

        System.out.println();
        System.out.println("pair\tPr");
        for (RateCombo rateCombo : pairRateCombos) {

            double posterior = (double) rateCombo.count / (double) total;

            System.out.println(rateCombo.combo + "\t" + format.format(posterior));
        }

        List<RateCombo> singleRateCombos = new ArrayList<RateCombo>();
        for (int nodeNumber : nodeRates.keySet()) {
            int count = nodeRates.get(nodeNumber).size();

            if (count > 1) {
                double meanRate = getMeanRate(nodeNumber);

                double[] cpd = getRateCPD(nodeNumber);

                String details = format.format(meanRate) +
                        "\t" + format.format(cpd[0]) + "\t" + format.format(cpd[1]);

                singleRateCombos.add(
                        new RateCombo(nodeNumber + "\t" + details, count));
            }
        }

        Collections.sort(singleRateCombos, rateComboComparator);

        System.out.println();
        System.out.println("node change\trate\t95% lower\t95% upper\tPr\tBF");
        for (RateCombo rateCombo : singleRateCombos) {

            double posterior = (double) rateCombo.count / (double) total;
            double bf = bayesFactor(posterior, priorChange);

            System.out.println(rateCombo.combo + "\t" + format.format(posterior) + "\t" + format.format(bf));
        }
    }

    /**
     * @param posterior the posterior probability
     * @param prior     the prior probability of a *single* change
     * @return
     */
    private static double bayesFactor(double posterior, double prior) {

        double like1 = posterior / prior;
        double like2 = (1 - posterior) / (1 - prior);

        return like1 / like2;
    }

    private static void putCombo(List<Integer> changes) {

        if (changes.size() == 0) return;

        // put pairs first
        for (int i = 0; i < changes.size(); i++) {
            for (int j = i + 1; j < changes.size(); j++) {
                putPair(changes.get(i), changes.get(j));
            }
        }

        StringBuilder builder = new StringBuilder();
        builder.append(changes.get(0));
        for (int i = 1; i < changes.size(); i++) {
            builder.append("-");
            builder.append(changes.get(i));
        }

        String combo = builder.toString();

        Integer i = comboCounts.get(combo);
        if (i == null) {
            comboCounts.put(combo, 1);
        } else {
            comboCounts.put(combo, i + 1);
        }
    }

    private static void putPair(int i, int j) {

        String pair = "{" + i + "," + j + "}";
        Integer count = pairCounts.get(pair);
        if (count == null) {
            pairCounts.put(pair, 1);
        } else {
            pairCounts.put(pair, count + 1);
        }
    }

    private static void putRate(int nodeNumber, double rate) {

        List<Double> rates = nodeRates.get(nodeNumber);

        if (rates == null) {
            rates = new ArrayList<Double>();
            nodeRates.put(nodeNumber, rates);
        }

        rates.add(rate);
    }

    private static double getStdevRate(int nodeNumber) {
        List<Double> rates = nodeRates.get(nodeNumber);

        double[] r = new double[rates.size()];
        for (int i = 0; i < r.length; i++) {
            r[i] = rates.get(i);
        }

        return DiscreteStatistics.stdev(r);
    }

    private static double[] getRateCPD(int nodeNumber) {
        List<Double> rates = nodeRates.get(nodeNumber);

        double[] r = new double[rates.size()];
        for (int i = 0; i < r.length; i++) {
            r[i] = rates.get(i);
        }

        return new double[]{DiscreteStatistics.quantile(0.025, r), DiscreteStatistics.quantile(0.975, r)};
    }

    private static double getMeanRate(int nodeNumber) {
        List<Double> rates = nodeRates.get(nodeNumber);
        double meanRate = 0.0;
        for (double rate : rates) {
            meanRate += rate;
        }
        meanRate /= rates.size();
        return meanRate;
    }

    static class RateCombo {

        public RateCombo(String combo, int count) {
            this.combo = combo;
            this.count = count;
        }

        final String combo;
        final int count;

        public boolean equals(Object o) {
            assert o instanceof RateCombo;
            RateCombo c2 = (RateCombo) o;
            return (combo.equals(c2.combo) && count == c2.count);
        }
    }

    static Comparator<RateCombo> rateComboComparator = new Comparator<RateCombo>() {

        public int compare(RateCombo r1, RateCombo r2) {
            return r2.count - r1.count;
        }
    };
}

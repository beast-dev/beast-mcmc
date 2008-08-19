package dr.evomodel.randomYule;

import dr.stats.DiscreteStatistics;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.*;

/**
 * @author Alexei Drummond
 */
public class RandomYuleAnalyzer {

    static TreeMap<Integer, List<Double>> nodeRates = new TreeMap<Integer, List<Double>>();
    static Map<String, Integer> comboCounts = new HashMap<String, Integer>();
    static Map<String, Integer> pairCounts = new HashMap<String, Integer>();

    static NumberFormat format = NumberFormat.getNumberInstance(Locale.ENGLISH);

    public static void main(String[] args) throws IOException {

        format.setMaximumFractionDigits(3);

        BufferedReader reader = new BufferedReader(new FileReader(args[0]));

        // read header line
        reader.readLine();

        String line = reader.readLine();

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
            line = reader.readLine();
        }

        List<RateCombo> rateCombos = new ArrayList<RateCombo>();
        for (String combo : comboCounts.keySet()) {

            int count = comboCounts.get(combo);
            if (count > 1) rateCombos.add(new RateCombo(combo, count));
        }

        Collections.sort(rateCombos, rateComboComparator);

        System.out.println("combo\tfrequency");
        for (RateCombo rateCombo : rateCombos) {
            System.out.println(rateCombo.combo + "\t" + rateCombo.count);
        }

        List<RateCombo> pairRateCombos = new ArrayList<RateCombo>();
        for (String pair : pairCounts.keySet()) {
            int count = pairCounts.get(pair);
            if (count > 1) pairRateCombos.add(new RateCombo(pair, count));
        }

        Collections.sort(pairRateCombos, rateComboComparator);

        System.out.println();
        System.out.println("pair\tfrequency");
        for (RateCombo rateCombo : pairRateCombos) {
            System.out.println(rateCombo.combo + "\t" + rateCombo.count);
        }

        List<RateCombo> singleRateCombos = new ArrayList<RateCombo>();
        for (int nodeNumber : nodeRates.keySet()) {
            int count = nodeRates.get(nodeNumber).size();

            if (count > 1) {
                double meanRate = getMeanRate(nodeNumber);
                double stdev = getStdevRate(nodeNumber);

                String details = format.format(meanRate) + "\t" + format.format(stdev);

                singleRateCombos.add(
                        new RateCombo(nodeNumber + "\t" + details, count));
            }
        }

        Collections.sort(singleRateCombos, rateComboComparator);

        System.out.println();
        System.out.println("node change\trate\tstdev\tfrequency");
        for (RateCombo rateCombo : singleRateCombos) {
            System.out.println(rateCombo.combo + "\t" + rateCombo.count);
        }
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
            RateCombo c2 = (RateCombo) o;
            return (combo.equals(c2.combo) && count == c2.count);
        }
    }

    static Comparator rateComboComparator = new Comparator() {

        public int compare(Object o, Object o1) {
            RateCombo r1 = (RateCombo) o;
            RateCombo r2 = (RateCombo) o1;

            return r2.count - r1.count;
        }
    };
}

/**
 *
 */
package dr.evomodel.tree;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

/**
 * @author shhn001
 */
public class Convergence {


    double distance = Double.MAX_VALUE;

    public final static String LOG_CONVERGENCE = "LogConvergence";

    public final static String TREE_FILE_NAME = "treeFilename";

    public final static String REFERENCE_FILE_NAME = "referenceFilename";

    public final static String CHECK_EVERY = "checkEvery";

    private Tree tree = null;

    private long trees = 0;

    private HashMap<BitSet, Double> referenceCladeFrequencies = null;

    private HashMap<BitSet, Double> cladeOccurences = null;

    private HashMap<String, Integer> taxonMap;

    private String[] taxonArray;

    private String referenceTreeFileName;

    private double logEvery;

    /**
     *
     */
    public Convergence(Tree tree, int logEvery, String referenceTreeFileName) {

        this.logEvery = logEvery;
        this.referenceTreeFileName = referenceTreeFileName;
        this.tree = tree;

        taxonMap = new HashMap<String, Integer>();
        cladeOccurences = new HashMap<BitSet, Double>();

        startLogging();
    }

    private void startLogging() {

        parseReferenceFile(referenceTreeFileName);

        taxonArray = new String[taxonMap.size()];
        Set<String> taxon = taxonMap.keySet();
        for (String s : taxon) {
            taxonArray[taxonMap.get(s)] = s;
        }
    }

    private void parseReferenceFile(String log) {
        FileReader fr;
        BufferedReader br;
        String line;
        referenceCladeFrequencies = new HashMap<BitSet, Double>();

        try {
            fr = new FileReader(log);
            br = new BufferedReader(fr);

            line = br.readLine();
            String[] tokens = line.split("\\s++");
            long total = Long.parseLong(tokens[0]);

            while ((line = br.readLine()) != null) {
                String[] taxa = getTaxa(line);
                BitSet split = new BitSet();
                for (String t : taxa) {
                    if (!taxonMap.containsKey(t)) {
                        taxonMap.put(t, taxonMap.size());
                    }
                    split.set(taxonMap.get(t));
                }
                tokens = line.split("\\s++");
                double frequency = (100.0 * Long.parseLong(tokens[0])) / (double) total;
                referenceCladeFrequencies.put(split,
                        frequency);
            }

            br.close();
            fr.close();
        } catch (FileNotFoundException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    private String[] getTaxa(String line) {
        line = line.substring(line.indexOf('(') + 1, line.indexOf(')'));
        String[] taxa = line.split(",");
        return taxa;
    }

    public double log(long state) {

        if (logEvery <= 0 || ((state % logEvery) == 0)) {

            addTree(tree);

            distance = getMaxCladeDistance();
        }

        return distance;
    }

    private void addTree(Tree tree) {
        trees++;
        List<BitSet> splits = getSplits(tree);

        for (BitSet s : splits) {
            if (!cladeOccurences.containsKey(s)) {
                cladeOccurences.put(s, 1.0);
            } else {
                cladeOccurences.put(s, cladeOccurences.get(s) + 1);
            }
        }
    }

    private List<BitSet> getSplits(Tree tree) {
        List<BitSet> splits = new ArrayList<BitSet>();

        NodeRef root = tree.getRoot();
        fillSplits(splits, root, tree);

        return splits;
    }

    private BitSet fillSplits(List<BitSet> splits, NodeRef root, Tree tree) {
        BitSet split = new BitSet();

        if (!tree.isExternal(root)) {
            split.or(fillSplits(splits, tree.getChild(root, 0), tree));
            split.or(fillSplits(splits, tree.getChild(root, 1), tree));
            splits.add(split);
        } else {
            Taxon taxon = tree.getNodeTaxon(root);
            String name = taxon.getId();
            split.set(taxonMap.get(name));
        }

        return split;

    }

    private double getMaxCladeDistance() {

        return getMaxDeviation();
    }

    private double getMaxDeviation() {
        double[] deviation = getDeviations();

        double max = 0;
        for (double m : deviation) {
            if (m > max) {
                max = m;
            }
        }
        return max;
    }

    private double[] getDeviations() {
        double[] deviations = new double[Math.max(referenceCladeFrequencies
                .size(), cladeOccurences.size())];

        Set<BitSet> keys = null;
        if (referenceCladeFrequencies.size() < cladeOccurences.size()) {
            keys = cladeOccurences.keySet();
        } else {
            keys = referenceCladeFrequencies.keySet();
        }

        Iterator<BitSet> it = keys.iterator();
        for (int i = 0; i < deviations.length; i++) {
            BitSet k = it.next();
            double r = 0.0;
            if (referenceCladeFrequencies.containsKey(k)) {
                r = referenceCladeFrequencies.get(k);
            }
            double o = 0.0;
            if (cladeOccurences.containsKey(k)) {
                o = cladeOccurences.get(k);
            }
            deviations[i] = Math.abs(r - ((o * 100.0) / (double) trees));
        }

        return deviations;
    }

}

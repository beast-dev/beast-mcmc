/**
 *
 */
package dr.evomodel.tree;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.inference.loggers.LogFormatter;
import dr.inference.loggers.MCLogger;
import dr.inference.loggers.MLLogger;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.inferencexml.loggers.LoggerParser;
import dr.xml.*;

import java.io.*;
import java.util.*;

/**
 * @author shhn001
 */
public class ConvergenceLogger extends MCLogger {

    public final static String LOG_CONVERGENCE = "LogConvergence";

    public final static String TREE_FILE_NAME = "treeFilename";

    public final static String REFERENCE_FILE_NAME = "referenceFilename";

    public final static String CHECK_EVERY = "checkEvery";

    private Tree tree = null;

    private BufferedWriter bwDistances = null;

    private long trees = 0;

    private HashMap<BitSet, Double> referenceCladeFrequencies = null;

    private HashMap<BitSet, Double> cladeOccurences = null;

    private String outputFilename;

    private HashMap<String, Integer> taxonMap;

    private String[] taxonArray;

    private String referenceTreeFileName;

    public ConvergenceLogger(Tree tree, LogFormatter formatter, int logEvery,
                             String outputFilename, String referenceTreeFileName) {

        super(formatter, logEvery, false);

        this.referenceTreeFileName = referenceTreeFileName;
        this.tree = tree;
        this.outputFilename = outputFilename.substring(0, outputFilename
                .indexOf('.'))
                + ".dist";

        taxonMap = new HashMap<String, Integer>();
        cladeOccurences = new HashMap<BitSet, Double>();
    }

    public void startLogging() {

        File f = new File(outputFilename.substring(0, outputFilename
                .indexOf('.'))
                + ".dist");
        if (f.exists()) {
            f.delete();
        }

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

//			do {
//				line = br.readLine();
//				if (line != null) {
//					line = line.trim();
//					if (cladeSection && !line.equals("")) {
//						if (line.startsWith("Clade")) {
//							cladeSection = false;
//						} else {
//							String[] taxa = getTaxa(line);
//							BitSet split = new BitSet();
//							for (String t : taxa) {
//								if (t.indexOf(',') != -1) {
//									t = t.substring(0, t.indexOf(','));
//								}
//								if (!taxonMap.containsKey(t)) {
//									taxonMap.put(t, taxonMap.size());
//								}
//								split.set(taxonMap.get(t));
//							}
//							referenceSplitFrequencies.put(split,
//									getPercentage(line));
//						}
//					} else if (line.contains("%-rule clades")) {
//						cladeSection = true;
//					}
//				}
//			} while (line != null);

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

    public void log(long state) {

        if (logEvery <= 0 || ((state % logEvery) == 0)) {

            addTree(tree);

            if (state % (logEvery * 10) == 0) {
                double distance = getMaxCladeDistance();
//				String summary = getSummary();
                try {
                    if (bwDistances == null) {
                        bwDistances = new BufferedWriter(new FileWriter(
                                outputFilename, true));
                    }
//					bwDistances.append((state) + "\t" + distance + "\t" + summary);
                    bwDistances.append((state) + "\t" + distance);
                    bwDistances.newLine();
                    bwDistances.flush();
                } catch (IOException e) {
                    // TODO Auto-generated catch block
                    e.printStackTrace();
                }
            }
        }
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

//	private String getSummary() {
//		double maxDev = 0;
//		double expected = 0;
//		double achieved = 0;
//		BitSet maxDevClade = null;
//
//		Set<BitSet> keys = null;
//		if (referenceCladeFrequencies.size() < cladeOccurences.size()) {
//			keys = cladeOccurences.keySet();
//		} else {
//			keys = referenceCladeFrequencies.keySet();
//		}
//
//		double count = Math.max(referenceCladeFrequencies
//				.size(), cladeOccurences.size());
//		Iterator<BitSet> it = keys.iterator();
//		for (int i = 0; i < count; i++) {
//			BitSet k = it.next();
//			double r = 0.0;
//			if (referenceCladeFrequencies.containsKey(k)) {
//				r = referenceCladeFrequencies.get(k);
//			}
//			double o = 0.0;
//			if (cladeOccurences.containsKey(k)) {
//				o = cladeOccurences.get(k);
//			}
//			double tmp = Math.abs(r - ((o * 100.0) / (double) trees));
//			if (tmp > maxDev){
//				expected = r;
//				achieved = ((o * 100.0) / (double) trees);
//				maxDev = tmp;
//				maxDevClade = k;
//			}
//		}
//		
//		String summary = "(";
//		int oldIndex = 0;
//		for (int i=0; i<maxDevClade.cardinality(); i++){
//			int index = maxDevClade.nextSetBit(oldIndex);
//			oldIndex = index+1;
//			summary += taxonArray[index];
//			
//			if (i+1 < maxDevClade.cardinality()){
//				summary += ",";
//			}
//		}
//		expected = ((int)(expected * 100))/100.0;
//		achieved = ((int)(achieved * 100))/100.0;
//		summary += ") - expected: " + expected + " achieved: " + achieved;
//
//		return summary;
//	}

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

    public void stopLogging() {

        logLine("End;");
        super.stopLogging();
    }

    public static XMLObjectParser PARSER = new LoggerParser() {

        public String getParserName() {
            return LOG_CONVERGENCE;
        }

        /**
         * @return an object based on the XML element it was passed.
         */
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String treeFileName = xo.getStringAttribute(TREE_FILE_NAME);

            String referenceFileName = xo
                    .getStringAttribute(REFERENCE_FILE_NAME);

            int checkEvery = 1;
            if (xo.hasAttribute(CHECK_EVERY)) {
                checkEvery = xo.getIntegerAttribute(CHECK_EVERY);
            }

            Tree tree = (Tree) xo.getChild(Tree.class);
            PrintWriter pw = getLogFile(xo, getParserName());
            LogFormatter formatter = new TabDelimitedFormatter(pw);

            ConvergenceLogger logger = new ConvergenceLogger(tree, formatter,
                    checkEvery, treeFileName, referenceFileName);

            return logger;
        }

        // ************************************************************************
        // AbstractXMLObjectParser implementation
        // ************************************************************************
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new StringAttributeRule(TREE_FILE_NAME,
                        "name of a tree log file", "trees.log"),
                new StringAttributeRule(REFERENCE_FILE_NAME,
                        "name of a reference tree file", "trees.log"),
                AttributeRule.newIntegerRule(CHECK_EVERY, true),
                new ElementRule(TreeModel.class)};

        public String getParserDescription() {
            return "Checks the convergence in terms of distance to the reference run.";
        }

        public String getExample() {
            return "<!-- The " + getParserName()
                    + " element takes a treeModel to be logged -->\n" + "<"
                    + getParserName() + " " + LOG_EVERY + "=\"100\" "
                    + TREE_FILE_NAME + "=\"log.trees\" " + REFERENCE_FILE_NAME
                    + "=\"log.trees\" "
                    + "	<treeModel idref=\"treeModel1\"/>\n" + "</"
                    + getParserName() + ">\n";
        }

        public Class getReturnType() {
            return MLLogger.class;
        }
    };

}

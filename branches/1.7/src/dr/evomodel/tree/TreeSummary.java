/**
 *
 */
package dr.evomodel.tree;

import dr.evolution.tree.Tree;
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
public class TreeSummary extends MCLogger {

    public final static String LOG_TREE = "TreeSummary";

    public final static String OUTPUT_FILE_NAME = "file";

    public final static String CHECK_EVERY = "checkEvery";

    private Tree tree = null;

    private HashMap<String, Integer> taxa = null;

    private HashMap<String, Integer> treeOccurences = null;

    private String outputFilename;

    public TreeSummary(Tree tree, LogFormatter formatter, int logEvery,
                       String outputFilename) {

        super(formatter, logEvery, false);
        this.tree = tree;
        this.outputFilename = outputFilename;

        treeOccurences = new HashMap<String, Integer>();
    }

    public void startLogging() {

        File f = new File(outputFilename);
        if (f.exists()) {
            f.delete();
        }
    }

    private HashMap<String, Integer> getTaxa() {
        int n = tree.getTaxonCount();

        List<String> l = new ArrayList<String>();
        for (int i = 0; i < n; i++) {
            l.add(tree.getTaxonId(i));
        }
        Collections.sort(l);

        HashMap<String, Integer> map = new HashMap<String, Integer>();
        for (int i = 1; i <= n; i++) {
            map.put(l.get(i - 1), i);
        }

        return map;
    }

    public void log(long state) {

        if (logEvery <= 0 || ((state % logEvery) == 0)) {
            if (state == 0) {
                taxa = getTaxa();
            }

            addTree(tree);

        }
    }

    private void addTree(Tree tree) {
        String newick = Tree.Utils.uniqueNewick(tree, tree.getRoot());
        newick = replaceTaxa(newick);

        if (treeOccurences.containsKey(newick)) {
            treeOccurences.put(newick, treeOccurences.get(newick) + 1);
        } else {
            treeOccurences.put(newick, 1);
        }
    }

    private String replaceTaxa(String newick) {
        String s = newick;
        int start = 0;
        int end;

        for (int i = 0; i < s.length(); i++) {
            char c = s.charAt(i);

            if (c == '(') {
                start = i + 1;
            } else if (c == ')') {
                end = i;
                if (end - start > 0) {
                    String taxon = s.substring(start, end);
                    int index = taxa.get(taxon);
                    s = s.substring(0, start) + index + s.substring(end);
                    i = start + (int) Math.log10(index);
                }
                start = i + 1;
            } else if (c == ',') {
                end = i;
                if (end - start > 0) {
                    String taxon = s.substring(start, end);
                    int index = taxa.get(taxon);
                    s = s.substring(0, start) + index + s.substring(end);
                    i = start + (int) Math.log10(index);
                }

                start = i + 1;
            }
        }

        return s;
    }

    public void stopLogging() {
        printTrees();

        logLine("End;");
        super.stopLogging();
    }

    private void printTrees() {
        Set<String> keys = treeOccurences.keySet();

        List<Sample> samples = new ArrayList<Sample>();
        for (String t : keys) {
            samples.add(new Sample(t, treeOccurences.get(t)));
        }
        Collections.sort(samples);

        try {
            FileWriter fw = new FileWriter(outputFilename);
            BufferedWriter writer = new BufferedWriter(fw);

            writer.write("Taxa");
            writer.newLine();
            writer.newLine();

            Set<String> taxon = taxa.keySet();
            for (String t : taxon) {
                int i = taxa.get(t);
                writer.write(i + "\t=\t" + t);
                writer.newLine();
                writer.flush();
            }

            writer.newLine();
            writer.newLine();
            writer.newLine();

            for (Sample s : samples) {
                writer.write(s.samples + "\t" + s.tree);
                writer.newLine();
                writer.flush();
            }

            writer.close();
            fw.close();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    public static XMLObjectParser PARSER = new LoggerParser() {

        public String getParserName() {
            return LOG_TREE;
        }

        /**
         * @return an object based on the XML element it was passed.
         */
        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String outputFile = xo.getStringAttribute(OUTPUT_FILE_NAME);

            int checkEvery = 1;
            if (xo.hasAttribute(CHECK_EVERY)) {
                checkEvery = xo.getIntegerAttribute(CHECK_EVERY);
            }

            Tree tree = (Tree) xo.getChild(Tree.class);

            final PrintWriter pw = getLogFile(xo, getParserName());

            LogFormatter formatter = new TabDelimitedFormatter(pw);

            return new TreeSummary(tree, formatter, checkEvery, outputFile);
        }

        // ************************************************************************
        // AbstractXMLObjectParser implementation
        // ************************************************************************
        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new StringAttributeRule(OUTPUT_FILE_NAME,
                        "name of a tree log file", "ds.trees"),
                AttributeRule.newIntegerRule(CHECK_EVERY, true),};

        public String getParserDescription() {
            return "Calculates the tree probabilities on the flow.";
        }

        public String getExample() {
            return "<!-- The " + getParserName()
                    + " element takes a treeModel to be logged -->\n" + "<"
                    + getParserName() + " " + LOG_EVERY + "=\"100\" "
                    + OUTPUT_FILE_NAME + "=\"log.trees\" "
                    + "	<treeModel idref=\"treeModel1\"/>\n" + "</"
                    + getParserName() + ">\n";
        }

        public Class getReturnType() {
            return MLLogger.class;
        }
    };

    class Sample implements Comparable<Sample> {

        String tree;
        int samples;

        public Sample(String tree, int samples) {
            super();
            this.tree = tree;
            this.samples = samples;
        }

        /* (non-Javadoc)
           * @see java.lang.Comparable#compareTo(java.lang.Object)
           */
        public int compareTo(Sample o) {
            return o.samples - samples;
        }

    }

}

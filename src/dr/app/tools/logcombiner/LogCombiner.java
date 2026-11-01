/*
 * LogCombiner.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.app.tools.logcombiner;

import dr.app.beast.BeastVersion;
import dr.evolution.tree.TreeUtils;
import dr.app.util.Arguments;
import dr.app.util.MacDesktopIntegration;
import dr.app.util.OSType;
import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.TreeImporter;
import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.util.Version;
import jam.console.ConsoleApplication;

import javax.swing.*;
import java.io.*;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Locale;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * @author Andrew Rambaut
 * @author Alexei Drummond
 */
public class LogCombiner {

    private final static Version version = new BeastVersion();

    /**
     * How the burn-in value given for each file is to be interpreted.
     */
    public enum BurninType {
        /** the burn-in is a number of MCMC states: samples with a state number below it are discarded */
        STATES("states"),
        /** the burn-in is a number of samples (trees or log lines) to discard from the start of each file */
        SAMPLES("samples"),
        /** the burn-in is a percentage of the samples in each file to discard from its start */
        PERCENT("percent");

        BurninType(String label) {
            this.label = label;
        }

        public String getLabel() {
            return label;
        }

        private final String label;
    }

    public LogCombiner(long[] burnins, long resample, String[] inputFileNames, String outputFileName, boolean treeFiles,
                       boolean convertToDecimal, boolean stripAnnotations,
                       boolean renumberOutput, boolean useScale, double scale) throws IOException {
        this(burnins, BurninType.STATES, resample, inputFileNames, outputFileName, treeFiles, convertToDecimal,
                stripAnnotations, renumberOutput, useScale, scale);
    }

    public LogCombiner(long[] burnins, BurninType burninType, long resample, String[] inputFileNames, String outputFileName,
                       boolean treeFiles, boolean convertToDecimal, boolean stripAnnotations,
                       boolean renumberOutput, boolean useScale, double scale) throws IOException {

        System.out.println("Creating combined " + (treeFiles ? "tree" : "log") + " file: '" + outputFileName + "'");

        if (convertToDecimal) {
            System.out.println("Converting to decimal.");
        }

        if (treeFiles && stripAnnotations) {
            System.out.println("Stripping annotations from trees.");
        }

        if (renumberOutput) {
            System.out.println("Renumbering output.");
        }

        if (useScale) {
            System.out.println("Rescaling using scale factor: " + scale);
        }

        System.out.println();

        PrintWriter writer = new PrintWriter(new FileOutputStream(outputFileName));

        boolean firstFile = true;
        boolean firstTree = true;
        long stateCount = (renumberOutput ? -1 : 0);
        long stateStep = -1;
        int columnCount = 0;

        String[] titles = null;

        System.out.println();
        for (int i = 0; i < inputFileNames.length; i++) {
            File inputFile = new File(inputFileNames[i]);

            if (!inputFile.exists()) {
                System.err.println(inputFileNames[i] + " does not exist!");
                return;
            } else if (inputFile.isDirectory()) {
                System.err.println(inputFileNames[i] + " is a directory.");
                return;
            }

            long burnin = burnins[0];
            if (burnins.length > i) {
                burnin = burnins[i];
            }

            Burnin fileBurnin = resolveBurnin(inputFile, burnin, burninType, treeFiles);
            burnin = fileBurnin.stateThreshold;
            System.out.print("Combining file: '" + inputFileNames[i] + "'" + fileBurnin.description);

            long samplesRead = 0;
            long samplesKept = 0;

            if (resample > 0) {
                System.out.print(", resampling with frequency: " + resample);
            }

            if (useScale) {
                System.out.println(", rescaling by: " + scale);
            } else {
                System.out.println();
            }

            boolean processTrees = stripAnnotations || convertToDecimal;

            if (treeFiles) {

                if (processTrees) { // then we need to read each tree and write it.
                    TreeImporter importer = new NexusImporter(new FileReader(inputFile), stripAnnotations);
                    try {
                        while (importer.hasTree()) {
                            Tree tree = importer.importNextTree();
                            if (firstTree) {
                                startLog(tree, writer);
                                firstTree = false;
                            }

                            String name = tree.getId();
                            if (name == null) {
                                System.err.println("ERROR: Trees do not give state numbers as tree attributes.");
                                return;
                            }

                            // split on underscore in STATE_xxxx
                            String[] bits = name.split("_");
                            long state = Long.parseLong(bits[1]);

                            if (stateStep < 0 && state > 0) {
                                stateStep = state;
                            }

                            samplesRead++;
                            if (state >= burnin) {
                                if (stateStep > 0) {
                                    if (!renumberOutput) {
                                        stateCount += stateStep;
                                    } else {
                                        stateCount += 1;
                                    }
                                }

                                if (resample >= 0) {
                                    if (resample % stateStep != 0) {
                                        System.err.println("ERROR: Resampling frequency is not a multiple of existing sampling frequency");
                                        return;
                                    }
                                }

                                boolean logThis;
                                if (resample < 0) {
                                    // not resampling, log every state
                                    logThis = true;
                                } else if (!renumberOutput) {
                                    // resampling but not renumbering
                                    logThis = (stateCount % resample == 0);
                                } else {
                                    logThis = (stateCount * stateStep % resample == 0);
                                }

                                long stateLineEntry;
                                if (!renumberOutput) {
                                    stateLineEntry = stateCount;
                                } else {
                                    stateLineEntry = stateCount / (resample / stateStep);
                                }

                                if (logThis) {
                                    samplesKept++;
                                    writeTree(stateLineEntry, tree, convertToDecimal, writer);
                                }

                            }
                        }
                    } catch (Importer.ImportException e) {
                        System.err.println("Error Parsing Input Tree: " + e.getMessage());
                        return;
                    }
                } else {
                    BufferedReader reader = new BufferedReader(new FileReader(inputFile));

                    String line = reader.readLine();

                    // skip (or write) the headers
                    while (line != null && !line.trim().startsWith("tree ")) {
                        if (firstFile) {
                            writer.println(line);
                        }
                        line = reader.readLine();
                    }

                    Pattern pattern = Pattern.compile("tree STATE_(\\d+)(\\s.*)");

                    while (line != null) {
                        Matcher m = pattern.matcher(line);
                        if (m.matches()) {

                            long state = Long.parseLong(m.group(1));

                            String treeString = m.group(2);

                            boolean skip = false;
                            if (!skip) {
                                if (stateStep < 0 && state > 0) {
                                    stateStep = state;
                                }

                                samplesRead++;
                                if (state >= burnin) {
                                    if (stateStep > 0) {
                                        if (!renumberOutput) {
                                            stateCount += stateStep;
                                        } else {
                                            stateCount += 1;
                                        }
                                    }

                                    if (resample >= 0) {
                                        if (resample % stateStep != 0) {
                                            System.err.println("ERROR: Resampling frequency is not a multiple of existing sampling frequency");
                                            return;
                                        }
                                    }

                                    boolean logThis;
                                    if (resample < 0) {
                                        logThis = true;
                                    } else if (!renumberOutput) {
                                        logThis = (stateCount % resample == 0);
                                    } else {
                                        logThis = ((stateCount * stateStep) % resample == 0);
                                    }

                                    long stateLineEntry;
                                    if (!renumberOutput) {
                                        stateLineEntry = stateCount;
                                    } else {

//                                	System.out.println("stateCount: " + stateCount);
//                                	System.out.println("resample: " + resample);
//                                	System.out.println("stateStep: " + stateStep);
//                                	System.out.println("resample / stateStep: " + (resample / stateStep));

                                        stateLineEntry = stateCount / (resample / stateStep);
                                    }

                                    if (logThis) {
                                        samplesKept++;
                                        writer.print("tree STATE_");
                                        writer.print(stateLineEntry);
                                        writer.println(treeString);
                                    }
                                }
                            }
                        }
                        line = reader.readLine();
                        //lineCount++;
                    }
                }

                firstFile = false;

            } else {
                BufferedReader reader = new BufferedReader(new FileReader(inputFile));
                //int lineCount = 1;
                String line = reader.readLine();

                // lines starting with [ are ignored, assuming comments in MrBayes file
                // lines starting with # are ignored, assuming comments in Migrate or BEAST file
                while (line.startsWith("[") || line.startsWith("#")) {

                    line = reader.readLine();
                }


                if (firstFile) {
                    titles = line.split("\t");
                    writer.println(line);
                } else {
                    String[] newTitles = line.split("\t");
                    if (newTitles.length != titles.length) {
                        System.err.println("ERROR: The number of columns in file, " + inputFileNames[i] + ", does not match that of the first file");
                        return;
                    }
                    for (int k = 0; k < newTitles.length; k++) {
                        if (!newTitles[k].equals(titles[k])) {
                            System.err.println("WARNING: The column heading, " + newTitles[k] + " in file, " + inputFileNames[i] + ", does not match the first file's heading, " + titles[k]);
                        }
                    }
                }
                line = reader.readLine();
                //lineCount++;

                while (line != null) {
                    String[] parts = line.split("\t");

                    long state = -1;

                    boolean skip = false;
                    try {
                        state = Long.parseLong(parts[0]);
                    } catch (NumberFormatException nfe) {
                        skip = true;
                    }

                    if (!skip) {
                        samplesRead++;
                        if (columnCount == 0) {
                            // the first sample line sets the expected number of columns (it used to be set with
                            // the state step, from the second line, so the first sample of the first file was lost)
                            columnCount = parts.length;
                        }
                        if (stateStep < 0 && state > 0) {
                            stateStep = state;
                        }

                        // if the columnCount is not the same then perhaps the line is corrupt so skip it.
                        if (state >= burnin && parts.length == columnCount) {
                            for (int j = 1; j < parts.length; j++) {
                                try {
                                    if (!parts[j].startsWith("{")) {
                                        // complex log values start with a curly bracket - otherwise attempt to parse
                                        // it as a number. If it fails, skip the line as a possible corruption.

                                        // attempt to convert the column value...
                                        double value = Double.valueOf(parts[j]);
                                    }
                                } catch (NumberFormatException nfe) {
                                    skip = true;
                                    break;
                                }
                            }

                            if (!skip) {

                                if (stateStep > 0) {
                                    if (!renumberOutput) {
                                        stateCount += stateStep;
                                    } else {
                                        stateCount += 1;
                                    }
                                }

                                if (resample >= 0) {
                                    if (resample % stateStep != 0) {
                                        System.err.println("ERROR: Resampling frequency is not a multiple of existing sampling frequency");
                                        return;
                                    }
                                }

                                boolean logThis;
                                if (resample < 0) {
                                    logThis = true;
                                } else if (!renumberOutput){
                                    logThis = (stateCount % resample == 0);
                                } else {
                                    logThis = ((stateCount * stateStep) % resample == 0);
                                }

                                long stateLineEntry;
                                if (!renumberOutput){
                                    stateLineEntry = stateCount;
                                } else {

//                                	System.out.println("stateCount: " + stateCount);
//                                	System.out.println("resample: " + resample);
//                                	System.out.println("stateStep: " + stateStep);
//                                	System.out.println("resample / stateStep: " + (resample / stateStep));

                                    stateLineEntry = stateCount / (resample / stateStep);
                                }

                                if (logThis) {
                                    samplesKept++;
                                    writer.print(stateLineEntry);
                                    for (int j = 1; j < parts.length; j++) {
                                        String value = parts[j];

                                        if (useScale) {
                                            if (titles[j].equals("clock.rate") || titles[j].startsWith("skyline.popSize")) {
                                                value = reformatNumbers(value, convertToDecimal, true, 1.0 / scale);
                                            } else if (titles[j].equals("treeModel.rootHeight")) {
                                                value = reformatNumbers(value, convertToDecimal, true, scale);
                                            }
                                        } else  if (convertToDecimal) {
                                            value = reformatNumbers(value, convertToDecimal, false, 1.0);
                                        }
                                        writer.print("\t" + value);
                                    }
                                    writer.println();
                                }
                            }

                        }
                    }
                    line = reader.readLine();
                    //lineCount++;
                }
            }

            firstFile = false;

            reportFile(inputFileNames[i], fileBurnin, samplesRead, samplesKept, resample, treeFiles);
        }

        if (treeFiles) {
            stopLog(writer);
        }
        writer.close();
    }

    /**
     * The burn-in for one file, resolved to a threshold on the state number.
     */
    private static class Burnin {
        /** samples with a state number below this are discarded */
        long stateThreshold;
        /** the number of samples in the file (counted only for SAMPLES and PERCENT) */
        long sampleCount = -1;
        /** the value as given by the user */
        long value;
        /** how the value was given */
        BurninType type;
        /** the text appended to the "Combining file" message */
        String description;
    }

    /**
     * Converts the burn-in given by the user into a threshold on the state number for one file. For a burn-in
     * given as a number of samples or a percentage the file is scanned first to find the state number of the
     * first sample to be kept, so that the main loop only ever compares state numbers.
     */
    private static Burnin resolveBurnin(File inputFile, long value, BurninType type, boolean treeFiles) throws IOException {
        Burnin burnin = new Burnin();
        burnin.value = value;
        burnin.type = type;

        if (type == BurninType.STATES) {
            burnin.stateThreshold = value;
            burnin.description = (value > 0 ? " removing burnin: " + value + " states" : " without removing burnin");
            return burnin;
        }

        java.util.List<Long> states = readSampleStates(inputFile, treeFiles);
        burnin.sampleCount = states.size();

        long samplesToRemove;
        if (type == BurninType.PERCENT) {
            samplesToRemove = Math.round(states.size() * (value / 100.0));
        } else {
            samplesToRemove = value;
        }

        String unit = treeFiles ? "trees" : "samples";
        if (samplesToRemove <= 0) {
            burnin.stateThreshold = -1;
            burnin.description = " without removing burnin";
        } else if (samplesToRemove >= states.size()) {
            burnin.stateThreshold = Long.MAX_VALUE;
            burnin.description = " removing burnin: all " + states.size() + " " + unit;
        } else {
            burnin.stateThreshold = states.get((int) samplesToRemove);
            burnin.description = " removing burnin: " + samplesToRemove + " " + unit +
                    (type == BurninType.PERCENT ? " (" + value + "%)" : "") +
                    ", i.e. states below " + burnin.stateThreshold;
        }
        return burnin;
    }

    /**
     * Reads the state numbers of the samples in a log or tree file, in file order.
     */
    private static java.util.List<Long> readSampleStates(File inputFile, boolean treeFiles) throws IOException {
        java.util.List<Long> states = new java.util.ArrayList<Long>();
        BufferedReader reader = new BufferedReader(new FileReader(inputFile));
        try {
            String line = reader.readLine();
            if (treeFiles) {
                Pattern pattern = Pattern.compile("tree STATE_(\\d+)\\s.*");
                while (line != null) {
                    Matcher m = pattern.matcher(line);
                    if (m.matches()) {
                        states.add(Long.parseLong(m.group(1)));
                    }
                    line = reader.readLine();
                }
            } else {
                // skip comments and the header line, as the main loop does
                while (line != null && (line.startsWith("[") || line.startsWith("#"))) {
                    line = reader.readLine();
                }
                if (line != null) {
                    line = reader.readLine();
                }
                while (line != null) {
                    int tab = line.indexOf('\t');
                    String first = (tab < 0 ? line : line.substring(0, tab));
                    try {
                        states.add(Long.parseLong(first));
                    } catch (NumberFormatException nfe) {
                        // not a sample line
                    }
                    line = reader.readLine();
                }
            }
        } finally {
            reader.close();
        }
        return states;
    }

    /**
     * Prints how many samples were read from a file and how many were written, warning when the burn-in
     * removed nothing or everything, which is usually a sign that it was given in the wrong units.
     */
    private static void reportFile(String fileName, Burnin burnin, long samplesRead, long samplesKept, long resample, boolean treeFiles) {
        String unit = treeFiles ? "trees" : "samples";
        System.out.println("  " + samplesRead + " " + unit + " read, " + samplesKept + " written");

        if (samplesRead == 0) {
            System.out.println("  WARNING: no " + unit + " were found in '" + fileName + "'" +
                    (treeFiles ? " - the trees must be in BEAST's 'tree STATE_<n> ... = ...' format" : ""));
        } else if (samplesKept == 0) {
            System.out.println("  WARNING: the burnin removed every " + (treeFiles ? "tree" : "sample") + " in '" + fileName + "'");
        } else if (burnin.value > 0 && burnin.type == BurninType.STATES && resample < 0) {
            long removed = samplesRead - samplesKept;
            if (removed <= 1 && samplesRead > 2) {
                System.out.println("  WARNING: the burnin of " + burnin.value + " states removed " +
                        (removed == 0 ? "nothing" : "only the first " + (treeFiles ? "tree" : "sample")) +
                        " from '" + fileName + "' - the burnin is a number of MCMC states, not of " + unit +
                        " (use -burninTrees or -burninPercent to give it as " + unit + " or as a percentage)");
            }
        }
    }

    private void rescaleTree(Tree tree, double scale) {
        if (tree instanceof MutableTree) {
            MutableTree mutableTree = (MutableTree) tree;

            for (int i = 0; i < tree.getNodeCount(); i++) {
                NodeRef node = tree.getNode(i);
                if (node != tree.getRoot()) {
                    double length = tree.getBranchLength(node);
                    mutableTree.setBranchLength(node, length * scale);
                }
            }
        } else {
            throw new IllegalArgumentException("Tree not mutable");
        }
    }

    private final Map<String, Integer> taxonMap = new HashMap<String, Integer>();

    private void startLog(Tree tree, PrintWriter writer) {

        int taxonCount = tree.getTaxonCount();
        writer.println("#NEXUS");
        writer.println("");
        writer.println("Begin taxa;");
        writer.println("\tDimensions ntax=" + taxonCount + ";");
        writer.println("\tTaxlabels");
        for (int i = 0; i < taxonCount; i++) {
            writer.println("\t\t" + TreeUtils.quoteTaxonName(tree.getTaxon(i).getId()));
        }
        writer.println("\t\t;");
        writer.println("End;");
        writer.println("");
        writer.println("Begin trees;");

        // This is needed if the trees use numerical taxon labels
        writer.println("\tTranslate");
        for (int i = 0; i < taxonCount; i++) {
            int k = i + 1;
            Taxon taxon = tree.getTaxon(i);
            taxonMap.put(taxon.getId(), k);
            writer.println("\t\t" + k + " " + TreeUtils.quoteTaxonName(taxon.getId()) + (k < taxonCount ? "," : ""));
        }
        writer.println("\t\t;");
    }

    private void writeTree(long state, Tree tree, boolean convertToDecimal, PrintWriter writer) {

        StringBuilder buffer = new StringBuilder("tree STATE_");
        buffer.append(state);
//        Double lnP = (Double) tree.getAttribute("lnP");
//        if (lnP != null) {
//            buffer.append(" [&lnP=").append(lnP).append("]");
//        }

        boolean hasAttribute = false;
        Iterator iter = tree.getAttributeNames();
        while (iter != null && iter.hasNext()) {
            String name = (String) iter.next();
            Object value = tree.getAttribute(name);

            if (!hasAttribute) {
                buffer.append(" [&");
                hasAttribute = true;
            } else {
                buffer.append(",");
            }

            buffer.append(name).append("=").append(formatValue(value));
        }

        if (hasAttribute) {
            buffer.append("]");
        }

        buffer.append(" = [&R] ");

        writeTree(tree, tree.getRoot(), taxonMap, convertToDecimal, buffer);

        buffer.append(";");
        writer.println(buffer.toString());
    }

    private String formatValue(Object value) {
        if( value instanceof String ) {
            return (String) value;
        } else if (value instanceof Object[] ) {
            String val = "{";
            for( Object v : (Object[]) value ) {
                val += formatValue(v);
                val += ',';
            }
            return val.substring(0, val.length() - 1 ) + '}';
        }
        return value.toString();
    }

    private void writeTree(Tree tree, NodeRef node, Map taxonMap, boolean convertToDecimal, StringBuilder buffer) {

        NodeRef parent = tree.getParent(node);

        if (tree.isExternal(node)) {
            String taxon = tree.getNodeTaxon(node).getId();
            Integer taxonNo = (Integer) taxonMap.get(taxon);
            if (taxonNo == null) {
                throw new IllegalArgumentException("Taxon, " + taxon + ", not recognized from first tree file");
            }
            buffer.append(taxonNo);
        } else {
            buffer.append("(");
            writeTree(tree, tree.getChild(node, 0), taxonMap, convertToDecimal, buffer);
            for (int i = 1; i < tree.getChildCount(node); i++) {
                buffer.append(",");
                writeTree(tree, tree.getChild(node, i), taxonMap, convertToDecimal, buffer);
            }
            buffer.append(")");
        }

        boolean hasAttribute = false;
        Iterator iter = tree.getNodeAttributeNames(node);
        while (iter != null && iter.hasNext()) {
            String name = (String) iter.next();
            Object value = tree.getNodeAttribute(node, name);

            if (!hasAttribute) {
                buffer.append("[&");
                hasAttribute = true;
            } else {
                buffer.append(",");
            }

            buffer.append(name).append("=").append(formatValue(value));
        }

        if (hasAttribute) {
            buffer.append("]");
        }

        if (parent != null) {
            buffer.append(":");
            double length = tree.getBranchLength(node);
            buffer.append(convertToDecimal ? decimalFormatter.format(length) : scientificFormatter.format(length));
        }
    }

    private void stopLog(PrintWriter writer) {
        writer.println("End;");
    }

    private static final DecimalFormat decimalFormatter = new DecimalFormat("#.############", new DecimalFormatSymbols(Locale.US));
    private static final DecimalFormat scientificFormatter = new DecimalFormat("#.############E0", new DecimalFormatSymbols(Locale.US));

    private String reformatNumbers(String line, boolean convertDecimal, boolean useScale, double scale) {
        StringBuilder outLine = new StringBuilder();

        Pattern pattern = Pattern.compile("\\d+\\.\\d+(E[\\-\\d\\.]+)?");
        Matcher matcher = pattern.matcher(line);

        int lastEnd = 0;
        while (matcher.find()) {
            int start = matcher.start();
            String token = matcher.group();
            double value = Double.parseDouble(token);
            if (useScale) {
                value *= scale;
            }
            String outToken = (convertDecimal ? decimalFormatter.format(value) : scientificFormatter.format(value));

            outLine.append(line.substring(lastEnd, start));
            outLine.append(outToken);

            lastEnd = matcher.end();
        }
        outLine.append(line.substring(lastEnd));
        return outLine.toString();
    }

    public static void printTitle() {
        System.out.println();
        centreLine("LogCombiner " + version.getVersionString() + ", " + version.getDateString(), 60);
        centreLine("MCMC Output Combiner", 60);
        centreLine("by", 60);
        centreLine("Andrew Rambaut and Alexei J. Drummond", 60);
        System.out.println();
        centreLine("Institute of Ecology and Evolution", 60);
        centreLine("University of Edinburgh", 60);
        centreLine("a.rambaut@ed.ac.uk", 60);
        System.out.println();
        centreLine("Department of Computer Science", 60);
        centreLine("University of Auckland", 60);
        centreLine("alexei@cs.auckland.ac.nz", 60);
        System.out.println();
        System.out.println();
    }

    public static void centreLine(String line, int pageWidth) {
        int n = pageWidth - line.length();
        int n1 = n / 2;
        for (int i = 0; i < n1; i++) {
            System.out.print(" ");
        }
        System.out.println(line);
    }


    public static void printUsage(Arguments arguments) {

        arguments.printUsage("logcombiner", "<input-file-name1> [<input-file-name2> ...] <output-file-name>");
        System.out.println();
        System.out.println("  Example: logcombiner test1.log test2.log combined.log");
        System.out.println("  Example: logcombiner -burnin 10000 test1.log test2.log combined.log");
        System.out.println("  Example: logcombiner -trees -burninPercent 10 test1.trees test2.trees combined.trees");
        System.out.println();

    }

    //Main method
    public static void main(String[] args) throws IOException {

        // There is a major issue with languages that use the comma as a decimal separator.
        // To ensure compatibility between programs in the package, enforce the US locale.
        Locale.setDefault(Locale.US);

        boolean treeFiles;
        boolean convertToDecimal;
        boolean stripAnnotations = false;
        boolean renumberOutput;

        long burnin;
        long resample = -1;
        double scale = 1.0;
        boolean useScale = false;

        if (args.length == 0) {
            // the menu bar at the top of the screen, and the application's name in it when run from a jar
            System.setProperty("apple.laf.useScreenMenuBar", "true");
            System.setProperty("apple.awt.application.name", "LogCombiner");

            java.net.URL url = LogCombiner.class.getResource("/images/utility.png");
            javax.swing.Icon icon = null;

            if (url != null) {
                icon = new javax.swing.ImageIcon(url);
            }

            final String versionString = version.getVersionString();
            String nameString = "LogCombiner " + versionString;
            String aboutString = "<html><center><p>" + versionString + ", " + version.getDateString() + "</p>" +
                    "<p>by<br>" +
                    "Andrew Rambaut and Alexei J. Drummond</p>" +
                    "<p>Institute of Ecology and Evolution, University of Edinburgh<br>" +
                    "<a href=\"mailto:a.rambaut@ed.ac.uk\">a.rambaut@ed.ac.uk</a></p>" +
                    "<p>Department of Computer Science, University of Auckland<br>" +
                    "<a href=\"mailto:alexei@cs.auckland.ac.nz\">alexei@cs.auckland.ac.nz</a></p>" +
                    "<p>Part of the BEAST package:<br>" +
                    "<a href=\"http://beast.community\">http://beast.community</a></p>" +
                    "</center></html>";

            ConsoleApplication consoleApp = new ConsoleApplication(nameString, aboutString, icon, true);
            if (OSType.isMac()) {
                // the application menu and Dock icon (jam's own adapter needs com.apple.eawt)
                MacDesktopIntegration.register(consoleApp, icon);
            }

            printTitle();

            LogCombinerDialog dialog = new LogCombinerDialog(new JFrame());

            if (!dialog.showDialog("LogCombiner " + versionString)) {
                return;
            }

            treeFiles = dialog.isTreeFiles();
            convertToDecimal = dialog.convertToDecimal();
            renumberOutput = dialog.renumberOutputStates();
            if (dialog.isResampling()) {
                resample = dialog.getResampleFrequency();
            }

            String[] inputFiles = dialog.getFileNames();
            long[] burnins = dialog.getBurnins();

            String outputFileName = dialog.getOutputFileName();

            if (outputFileName == null) {
                System.err.println("No output file specified");
            }

            try {
                new LogCombiner(burnins, dialog.getBurninType(), resample, inputFiles, outputFileName, treeFiles, convertToDecimal,
                        stripAnnotations, renumberOutput, useScale, scale);

            } catch (Exception ex) {
                System.err.println("Exception: " + ex.getMessage());
                ex.printStackTrace();
            }
            System.out.println("Finished - Quit program to exit.");
            while (true) {
                try {
                    Thread.sleep(1000);
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }
        } else {
            printTitle();

            Arguments arguments = new Arguments(
                    new Arguments.Option[]{
                            new Arguments.Option("trees", "t", "use this option to combine tree log files"),
                            new Arguments.Option("decimal", "d","this option converts numbers from scientific to decimal notation"),
                            new Arguments.LongOption("burnin", "b",  "the number of states to be considered as 'burn-in'"),
                            new Arguments.LongOption("burninTrees", "bt", "the number of trees (or log samples) at the start of each file to be considered as 'burn-in'"),
                            new Arguments.RealOption("burninPercent", "bp", 0.0, 100.0, "the percentage of each file to be considered as 'burn-in'"),
                            new Arguments.LongOption("resample", "r", "resample the log files to this frequency " +
                                    "(the original sampling frequency must be a factor of this value)"),
                            new Arguments.RealOption("scale", "s", "a scaling factor that will multiply any time units by this value"),
                            new Arguments.Option("strip", null, "strip out all annotations (trees only)"),
                            new Arguments.Option("renumber", "n", "this option renumbers output states consecutively"),
                            new Arguments.Option("help", "h", "option to print this message")
                    });

            try {
                arguments.parseArguments(args);
            } catch (Arguments.ArgumentException ae) {
                System.out.println(ae);
                printUsage(arguments);
                System.exit(1);
            }

            if (arguments.hasOption("help")) {
                printUsage(arguments);
                System.exit(0);
            }

            treeFiles = arguments.hasOption("trees");

            convertToDecimal = arguments.hasOption("decimal");

            stripAnnotations = arguments.hasOption("strip");

            renumberOutput = arguments.hasOption("renumber");

            burnin = -1;
            BurninType burninType = BurninType.STATES;
            int burninOptions = 0;
            if (arguments.hasOption("burnin")) {
                burnin = arguments.getLongOption("burnin");
                burninOptions++;
            }
            if (arguments.hasOption("burninTrees")) {
                burnin = arguments.getLongOption("burninTrees");
                burninType = BurninType.SAMPLES;
                burninOptions++;
            }
            if (arguments.hasOption("burninPercent")) {
                burnin = Math.round(arguments.getRealOption("burninPercent"));
                burninType = BurninType.PERCENT;
                burninOptions++;
            }
            if (burninOptions > 1) {
                System.err.println("Only one of -burnin, -burninTrees and -burninPercent can be given");
                System.err.println();
                printUsage(arguments);
                System.exit(1);
            }

            resample = -1;
            if (arguments.hasOption("resample")) {
                resample = arguments.getLongOption("resample");
            }

            scale = 1.0;
            useScale = false;
            if (arguments.hasOption("scale")) {
                scale = arguments.getRealOption("scale");
                useScale = true;
            }

            String[] args2 = arguments.getLeftoverArguments();

            if (args2.length < 2) {
                System.err.println("Requires a minimum of 1 input filename and 1 output filename");
                System.err.println();
                printUsage(arguments);
                System.exit(1);
            }

            String[] inputFileNames = new String[args2.length - 1];
            System.arraycopy(args2, 0, inputFileNames, 0, inputFileNames.length);
            String outputFileName = args2[args2.length - 1];

            new LogCombiner(new long[]{burnin}, burninType, resample, inputFileNames, outputFileName, treeFiles, convertToDecimal,
                    stripAnnotations, renumberOutput, useScale, scale);

            System.out.println("Finished.");
        }

        System.exit(0);
    }
}


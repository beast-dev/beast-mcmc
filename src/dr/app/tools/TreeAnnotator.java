/*
 * TreeAnnotator.java
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

package dr.app.tools;

import dr.app.beast.BeastVersion;
import dr.app.util.Arguments;
import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.TreeImporter;
import dr.evolution.tree.*;
import dr.evolution.util.TaxonList;
import dr.geo.contouring.ContourMaker;
import dr.geo.contouring.ContourPath;
import dr.geo.contouring.ContourWithSynder;
import dr.stats.DiscreteStatistics;
import dr.util.HeapSort;
import dr.util.Version;
import jam.console.ConsoleApplication;
import org.rosuda.JRI.REXP;
import org.rosuda.JRI.RVector;
import org.rosuda.JRI.Rengine;

import javax.swing.*;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 */
public class TreeAnnotator {

    private final static Version version = new BeastVersion();

    private final static boolean USE_R = false;

    private static boolean forceIntegerToDiscrete = false;

    enum Target {
        MAX_CLADE_CREDIBILITY("Maximum clade credibility tree"),
        //MAX_SUM_CLADE_CREDIBILITY("Maximum sum of clade credibilities"),
        USER_TARGET_TREE("User target tree");

        String desc;

        Target(String s) {
            desc = s;
        }

        public String toString() {
            return desc;
        }
    }

    enum HeightsSummary {
        MEDIAN_HEIGHTS("Median heights"),
        MEAN_HEIGHTS("Mean heights"),
        KEEP_HEIGHTS("Keep target heights"),
        CA_HEIGHTS("Common Ancestor heights");

        String desc;

        HeightsSummary(String s) {
            desc = s;
        }

        public String toString() {
            return desc;
        }
    }


    // Messages to stderr, output to stdout
    private static PrintStream progressStream = System.err;

    private final String location1Attribute = "longLat1";
    private final String location2Attribute = "longLat2";
    private final String locationOutputAttribute = "location";

    /**
     * Burnin can be specified as the number of trees or the number of states
     * (one or other should be zero).
     * @param burninTrees
     * @param burninStates
     * @param heightsOption
     * @param posteriorLimit
     * @param hpd2D
     * @param targetOption
     * @param targetTreeFileName
     * @param inputFileName
     * @param outputFileName
     * @throws IOException
     */
    public TreeAnnotator(final int burninTrees,
                         final int burninStates,
                         HeightsSummary heightsOption,
                         double posteriorLimit,
                         double[] hpd2D,
                         Target targetOption,
                         String targetTreeFileName,
                         String inputFileName,
                         String outputFileName
    ) throws IOException {

        this.posteriorLimit = posteriorLimit;
        this.hpd2D = hpd2D;

        attributeNames.add("height");
        attributeNames.add("length");

        CladeSystem cladeSystem = new CladeSystem();

        int burnin = -1;

        totalTrees = 10000;
        totalTreesUsed = 0;

        progressStream.println("Reading trees (bar assumes 10,000 trees)...");
        progressStream.println("0              25             50             75            100");
        progressStream.println("|--------------|--------------|--------------|--------------|");

        int stepSize = totalTrees / 60;
        if (stepSize < 1) stepSize = 1;

        if (targetOption != Target.USER_TARGET_TREE) {
            cladeSystem = new CladeSystem();
            FileReader fileReader = new FileReader(inputFileName);
            TreeImporter importer = new NexusImporter(fileReader, true);
            try {
                totalTrees = 0;
                while (importer.hasTree()) {
                    Tree tree = importer.importNextTree();

                    int state = Integer.MAX_VALUE;

                    if (burninStates > 0) {
                        // if burnin has been specified in states, try to parse it out...
                        String name = tree.getId().trim();

                        if (name != null && name.length() > 0 && name.startsWith("STATE_")) {
                            state = Integer.parseInt(name.split("_")[1]);
                        }
                    }

                    if (totalTrees >= burninTrees && state >= burninStates) {
                        // if either of the two burnin thresholds have been reached...

                        if (burnin < 0) {
                            // if this is the first time this point has been reached,
                            // record the number of trees this represents for future use...
                            burnin = totalTrees;
                        }

                        cladeSystem.add(tree, false);

                        totalTreesUsed += 1;
                    }

                    if (totalTrees > 0 && totalTrees % stepSize == 0) {
                        progressStream.print("*");
                        progressStream.flush();
                    }
                    totalTrees++;
                }

            } catch (Importer.ImportException e) {
                System.err.println("Error Parsing Input Tree: " + e.getMessage());
                return;
            }
            fileReader.close();
            progressStream.println();
            progressStream.println();

            if (totalTrees < 1) {
                System.err.println("No trees");
                return;
            }
            if (totalTreesUsed <= 1) {
                if (burnin > 0) {
                    System.err.println("No trees to use: burnin too high");
                    return;
                }
            }
            cladeSystem.calculateCladeCredibilities(totalTreesUsed);

            progressStream.println("Total trees read: " + totalTrees);
            if (burninTrees > 0) {
                progressStream.println("Ignoring first " + burninTrees + " trees" +
                        (burninStates > 0 ? " (" + burninStates + " states)." : "." ));
            } else if (burninStates > 0) {
                progressStream.println("Ignoring first " + burninStates + " states (" + burnin + " trees).");
            }

            progressStream.println("Total unique clades: " + cladeSystem.getCladeMap().keySet().size());
            progressStream.println();
        }

        MutableTree targetTree = null;

        switch (targetOption) {
            case USER_TARGET_TREE: {
                if (targetTreeFileName != null) {
                    progressStream.println("Reading user specified target tree, " + targetTreeFileName);

                    NexusImporter importer = new NexusImporter(new FileReader(targetTreeFileName));
                    try {
                        Tree tree = importer.importNextTree();
                        if (tree == null) {
                            NewickImporter x = new NewickImporter(new FileReader(targetTreeFileName));
                            tree = x.importNextTree();
                        }
                        if (tree == null) {
                            System.err.println("No tree in target nexus or newick file " + targetTreeFileName);
                            return;
                        }
                        targetTree = new FlexibleTree(tree);
                    } catch (Importer.ImportException e) {
                        System.err.println("Error Parsing Target Tree: " + e.getMessage());
                        return;
                    }
                } else {
                    System.err.println("No user target tree specified.");
                    return;
                }
                break;
            }
            case MAX_CLADE_CREDIBILITY: {
                progressStream.println("Finding maximum credibility tree...");
                targetTree = new FlexibleTree(summarizeTrees(burnin, cladeSystem, inputFileName /*, false*/));
                break;
            }
//            case MAX_SUM_CLADE_CREDIBILITY: {
//                progressStream.println("Finding maximum sum clade credibility tree...");
//                targetTree = new FlexibleTree(summarizeTrees(burnin, cladeSystem, inputFileName, true));
//                break;
//            }
            default: throw new IllegalArgumentException("Unknown targetOption");
        }

        progressStream.println("Collecting node information...");
        progressStream.println("0              25             50             75            100");
        progressStream.println("|--------------|--------------|--------------|--------------|");

        stepSize = totalTrees / 60;
        if (stepSize < 1) stepSize = 1;

        FileReader fileReader = new FileReader(inputFileName);
        NexusImporter importer = new NexusImporter(fileReader);

        // this call increments the clade counts and it shouldn't
        // this is remedied with removeClades call after while loop below
        cladeSystem = new CladeSystem(targetTree);
        totalTreesUsed = 0;
        try {
            boolean firstTree = true;
            int counter = 0;
            while (importer.hasTree()) {
                Tree tree = importer.importNextTree();

                if (counter >= burnin) {
                    if (firstTree) {
                        setupAttributes(tree);
                        firstTree = false;
                    }

                    cladeSystem.collectAttributes(tree);
                    totalTreesUsed += 1;
                }
                if (counter > 0 && counter % stepSize == 0) {
                    progressStream.print("*");
                    progressStream.flush();
                }
                counter++;

            }
            cladeSystem.removeClades(targetTree, targetTree.getRoot(), true);
            //progressStream.println("totalTreesUsed=" + totalTreesUsed);
            cladeSystem.calculateCladeCredibilities(totalTreesUsed);
        } catch (Importer.ImportException e) {
            System.err.println("Error Parsing Input Tree: " + e.getMessage());
            return;
        }
        progressStream.println();
        progressStream.println();
        fileReader.close();

        progressStream.println("Annotating target tree...");

        try {
            cladeSystem.annotateTree(targetTree, targetTree.getRoot(), null, heightsOption);

            if( heightsOption == HeightsSummary.CA_HEIGHTS ) {
                setTreeHeightsByCA(targetTree, inputFileName, burnin);
            }
        } catch (Exception e) {
            System.err.println("Error annotating tree: " + e.getMessage() + "\nPlease check the tree log file format.");
            return;
        }

        progressStream.println("Writing annotated tree....");

        try {
            final PrintStream stream = outputFileName != null ?
                    new PrintStream(new FileOutputStream(outputFileName)) :
                    System.out;

            new NexusExporter(stream).exportTree(targetTree);
        } catch (Exception e) {
            System.err.println("Error to write annotated tree file: " + e.getMessage());
            return;
        }

    }

    private void setupAttributes(Tree tree) {
        for (int i = 0; i < tree.getNodeCount(); i++) {
            NodeRef node = tree.getNode(i);
            Iterator iter = tree.getNodeAttributeNames(node);
            if (iter != null) {
                while (iter.hasNext()) {
                    String name = (String) iter.next();
                    attributeNames.add(name);
                }
            }
        }

        for (TreeAnnotationPlugin plugin : plugins) {
            Set<String> claimed = plugin.setAttributeNames(attributeNames);
            attributeNames.removeAll(claimed);
        }
    }

    private Tree summarizeTrees(int burnin, CladeSystem cladeSystem, String inputFileName /*, boolean useSumCladeCredibility */)
            throws IOException {

        Tree bestTree = null;
        double bestScore = Double.NEGATIVE_INFINITY;

        progressStream.println("Analyzing " + totalTreesUsed + " trees...");
        progressStream.println("0              25             50             75            100");
        progressStream.println("|--------------|--------------|--------------|--------------|");

        int stepSize = totalTrees / 60;
        if (stepSize < 1) stepSize = 1;

        int counter = 0;
        int bestTreeNumber = 0;
        TreeImporter importer = new NexusImporter(new FileReader(inputFileName), true);
        try {
            while (importer.hasTree()) {
                Tree tree = importer.importNextTree();

                if (counter >= burnin) {
                    double score = scoreTree(tree, cladeSystem /*, useSumCladeCredibility*/);
//                    progressStream.println(score);
                    if (score > bestScore) {
                        bestTree = tree;
                        bestScore = score;
                        bestTreeNumber = counter + 1;
                    }
                }
                if (counter > 0 && counter % stepSize == 0) {
                    progressStream.print("*");
                    progressStream.flush();
                }
                counter++;
            }
        } catch (Importer.ImportException e) {
            System.err.println("Error Parsing Input Tree: " + e.getMessage());
            return null;
        }
        progressStream.println();
        progressStream.println();
        progressStream.println("Best tree: " + bestTree.getId() + " (tree number " + bestTreeNumber + ")");
//        if (useSumCladeCredibility) {
//            progressStream.println("Highest Sum Clade Credibility: " + bestScore);
//        } else {
            progressStream.println("Highest Log Clade Credibility: " + bestScore);
//        }

        return bestTree;
    }

    private double scoreTree(Tree tree, CladeSystem cladeSystem /*, boolean useSumCladeCredibility*/) {
//        if (useSumCladeCredibility) {
//            return cladeSystem.getSumCladeCredibility(tree, tree.getRoot(), null);
//        } else {
            return cladeSystem.getLogCladeCredibility(tree, tree.getRoot(), null);
//        }
    }

    private class CladeSystem {
        //
        // Public stuff
        //

        /**
         */
        public CladeSystem() {
        }

        /**
         */
        public CladeSystem(Tree targetTree) {
            this.targetTree = targetTree;
            add(targetTree, true);
        }

        /**
         * adds all the clades in the tree
         */
        public void add(Tree tree, boolean includeTips) {
            if (taxonList == null) {
                taxonList = tree;
            }

            // Recurse over the tree and add all the clades (or increment their
            // frequency if already present). The root clade is added too (for
            // annotation purposes).
            addClades(tree, tree.getRoot(), includeTips);
        }
//
//        public Clade getClade(NodeRef node) {
//            return null;
//        }

        private BitSet addClades(Tree tree, NodeRef node, boolean includeTips) {

            BitSet bits = new BitSet();

            if (tree.isExternal(node)) {

                int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
                bits.set(index);

                if (includeTips) {
                    addClade(bits);
                }

            } else {

                for (int i = 0; i < tree.getChildCount(node); i++) {

                    NodeRef node1 = tree.getChild(node, i);

                    bits.or(addClades(tree, node1, includeTips));
                }

                addClade(bits);
            }

            return bits;
        }

        private void addClade(BitSet bits) {
            Clade clade = cladeMap.get(bits);
            if (clade == null) {
                clade = new Clade(bits);
                cladeMap.put(bits, clade);
            }
            clade.setCount(clade.getCount() + 1);
        }

        public void collectAttributes(Tree tree) {
            collectAttributes(tree, tree.getRoot());
        }

        private BitSet collectAttributes(Tree tree, NodeRef node) {

            BitSet bits = new BitSet();

            if (tree.isExternal(node)) {

                int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
                if (index < 0) {
                    throw new IllegalArgumentException("Taxon, " + tree.getNodeTaxon(node).getId() + ", not found in target tree");
                }
                bits.set(index);

            } else {

                for (int i = 0; i < tree.getChildCount(node); i++) {

                    NodeRef node1 = tree.getChild(node, i);

                    bits.or(collectAttributes(tree, node1));
                }
            }

            collectAttributesForClade(bits, tree, node);

            return bits;
        }

        private void collectAttributesForClade(BitSet bits, Tree tree, NodeRef node) {
            Clade clade = cladeMap.get(bits);
            if (clade != null) {

                if (clade.attributeValues == null) {
                    clade.attributeValues = new ArrayList<Object[]>();
                }

                int i = 0;
                Object[] values = new Object[attributeNames.size()];
                for (String attributeName : attributeNames) {
                    boolean processed = false;

                    if (!processed) {
                        Object value;
                        if (attributeName.equals("height")) {
                            value = tree.getNodeHeight(node);
                        } else if (attributeName.equals("length")) {
                            value = tree.getBranchLength(node);
// AR - we deal with this once everything                             
//                        } else if (attributeName.equals(location1Attribute)) {
//                            // If this is one of the two specified bivariate location names then
//                            // merge this and the other one into a single array.
//                            Object value1 = tree.getNodeAttribute(node, attributeName);
//                            Object value2 = tree.getNodeAttribute(node, location2Attribute);
//
//                            value = new Object[]{value1, value2};
//                        } else if (attributeName.equals(location2Attribute)) {
//                            // do nothing - already dealt with this...
//                            value = null;
                        } else {
                            value = tree.getNodeAttribute(node, attributeName);
                            if (value instanceof String && ((String) value).startsWith("\"")) {
                                value = ((String) value).replaceAll("\"", "");
                            }
                        }

                        //if (value == null) {
                        //    progressStream.println("attribute " + attributeNames[i] + " is null.");
                        //}

                        values[i] = value;
                    }
                    i++;
                }
                clade.attributeValues.add(values);

                //progressStream.println(clade + " " + clade.getValuesSize());
                clade.setCount(clade.getCount() + 1);
            }
        }

        public Map getCladeMap() {
            return cladeMap;
        }

        public void calculateCladeCredibilities(int totalTreesUsed) {
            for (Clade clade : cladeMap.values()) {

                if (clade.getCount() > totalTreesUsed) {

                    throw new AssertionError("clade.getCount=(" + clade.getCount() +
                            ") should be <= totalTreesUsed = (" + totalTreesUsed + ")");
                }

                clade.setCredibility(((double) clade.getCount()) / (double) totalTreesUsed);
            }
        }

//        public double getSumCladeCredibility(Tree tree, NodeRef node, BitSet bits) {
//
//            double sum = 0.0;
//
//            if (tree.isExternal(node)) {
//
//                int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
//                bits.set(index);
//            } else {
//
//                BitSet bits2 = new BitSet();
//                for (int i = 0; i < tree.getChildCount(node); i++) {
//
//                    NodeRef node1 = tree.getChild(node, i);
//
//                    sum += getSumCladeCredibility(tree, node1, bits2);
//                }
//
//                sum += getCladeCredibility(bits2);
//
//                if (bits != null) {
//                    bits.or(bits2);
//                }
//            }
//
//            return sum;
//        }

        public double getLogCladeCredibility(Tree tree, NodeRef node, BitSet bits) {

            double logCladeCredibility = 0.0;

            if (tree.isExternal(node)) {

                int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
                bits.set(index);
            } else {

                BitSet bits2 = new BitSet();
                for (int i = 0; i < tree.getChildCount(node); i++) {

                    NodeRef node1 = tree.getChild(node, i);

                    logCladeCredibility += getLogCladeCredibility(tree, node1, bits2);
                }

                logCladeCredibility += Math.log(getCladeCredibility(bits2));

                if (bits != null) {
                    bits.or(bits2);
                }
            }

            return logCladeCredibility;
        }

        private double getCladeCredibility(BitSet bits) {
            Clade clade = cladeMap.get(bits);
            if (clade == null) {
                return 0.0;
            }
            return clade.getCredibility();
        }

        public void annotateTree(MutableTree tree, NodeRef node, BitSet bits, HeightsSummary heightsOption) {

            BitSet bits2 = new BitSet();

            if (tree.isExternal(node)) {

                int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
                bits2.set(index);

                annotateNode(tree, node, bits2, true, heightsOption);
            } else {

                for (int i = 0; i < tree.getChildCount(node); i++) {

                    NodeRef node1 = tree.getChild(node, i);

                    annotateTree(tree, node1, bits2, heightsOption);
                }

                annotateNode(tree, node, bits2, false, heightsOption);
            }

            if (bits != null) {
                bits.or(bits2);
            }
        }

        private void annotateNode(MutableTree tree, NodeRef node, BitSet bits, boolean isTip, HeightsSummary heightsOption) {
            Clade clade = cladeMap.get(bits);
            assert clade != null : "Clade missing?";

            boolean filter = false;
            if (!isTip) {
                final double posterior = clade.getCredibility();
                tree.setNodeAttribute(node, "posterior", posterior);
                if (posterior < posteriorLimit) {
                    filter = true;
                }
            }

            int i = 0;
            for (String attributeName : attributeNames) {

                if (clade.attributeValues != null && clade.attributeValues.size() > 0) {
                    double[] values = new double[clade.attributeValues.size()];

                    HashMap<Object, Integer> hashMap = new HashMap<Object, Integer>();

                    Object[] v = clade.attributeValues.get(0);
                    if (v[i] != null) {

                        final boolean isHeight = attributeName.equals("height");
                        boolean isBoolean = v[i] instanceof Boolean;

                        boolean isDiscrete = v[i] instanceof String;

                        if (forceIntegerToDiscrete && v[i] instanceof Integer) isDiscrete = true;

                        double minValue = Double.MAX_VALUE;
                        double maxValue = -Double.MAX_VALUE;

                        final boolean isArray = v[i] instanceof Object[];
                        boolean isDoubleArray = isArray && ((Object[]) v[i])[0] instanceof Double;
                        // This is Java, friends - first value type does not imply all.
                        if (isDoubleArray) {
                            for (Object n : (Object[]) v[i]) {
                                if (!(n instanceof Double)) {
                                    isDoubleArray = false;
                                    break;
                                }
                            }
                        }
                        // todo Handle other types of arrays

                        double[][] valuesArray = null;
                        double[] minValueArray = null;
                        double[] maxValueArray = null;
                        int lenArray = 0;

                        if (isDoubleArray) {
                            lenArray = ((Object[]) v[i]).length;

                            valuesArray = new double[lenArray][clade.attributeValues.size()];
                            minValueArray = new double[lenArray];
                            maxValueArray = new double[lenArray];

                            for (int k = 0; k < lenArray; k++) {
                                minValueArray[k] = Double.MAX_VALUE;
                                maxValueArray[k] = -Double.MAX_VALUE;
                            }
                        }

                        for (int j = 0; j < clade.attributeValues.size(); j++) {
                            Object value = clade.attributeValues.get(j)[i];
                            if (isDiscrete) {
                                final Object s = value;
                                if (hashMap.containsKey(s)) {
                                    hashMap.put(s, hashMap.get(s) + 1);
                                } else {
                                    hashMap.put(s, 1);
                                }
                            } else if (isBoolean) {
                                values[j] = (((Boolean) value) ? 1.0 : 0.0);
                            } else if (isDoubleArray) {
                                // Forcing to Double[] causes a cast exception. MAS
                                Object[] array = (Object[]) value;
                                for (int k = 0; k < lenArray; k++) {
                                    valuesArray[k][j] = ((Double) array[k]);
                                    if (valuesArray[k][j] < minValueArray[k]) minValueArray[k] = valuesArray[k][j];
                                    if (valuesArray[k][j] > maxValueArray[k]) maxValueArray[k] = valuesArray[k][j];
                                }
                            } else {
                                // Ignore other (unknown) types
                                if (value instanceof Number) {
                                    values[j] = ((Number) value).doubleValue();
                                    if (values[j] < minValue) minValue = values[j];
                                    if (values[j] > maxValue) maxValue = values[j];
                                }
                            }
                        }
                        if (isHeight) {
                            if (heightsOption == HeightsSummary.MEAN_HEIGHTS) {
                                final double mean = DiscreteStatistics.mean(values);
                                tree.setNodeHeight(node, mean);
                            } else if (heightsOption == HeightsSummary.MEDIAN_HEIGHTS) {
                                final double median = DiscreteStatistics.median(values);
                                tree.setNodeHeight(node, median);
                            } else {
                                // keep the existing height
                            }
                        }

                        if (!filter) {
                            boolean processed = false;
                            for (TreeAnnotationPlugin plugin : plugins) {
                                if (plugin.handleAttribute(tree, node, attributeName, values)) {
                                    processed = true;
                                }
                            }

                            if (!processed) {
                                if (!isDiscrete) {
                                    if (!isDoubleArray)
                                        annotateMeanAttribute(tree, node, attributeName, values);
                                    else {
                                        for (int k = 0; k < lenArray; k++) {
                                            annotateMeanAttribute(tree, node, attributeName + (k + 1), valuesArray[k]);
                                        }
                                    }
                                } else {
                                    annotateModeAttribute(tree, node, attributeName, hashMap);
                                    annotateFrequencyAttribute(tree, node, attributeName, hashMap);
                                }
                                if (!isBoolean && minValue < maxValue && !isDiscrete && !isDoubleArray) {
                                    // Basically, if it is a boolean (0, 1) then we don't need the distribution information
                                    // Likewise if it doesn't vary.
                                    annotateMedianAttribute(tree, node, attributeName + "_median", values);
                                    annotateHPDAttribute(tree, node, attributeName + "_95%_HPD", 0.95, values);
                                    annotateRangeAttribute(tree, node, attributeName + "_range", values);
                                }

                                if (isDoubleArray) {
                                    String name = attributeName;
                                    // todo
//                                    if (name.equals(location1Attribute)) {
//                                        name = locationOutputAttribute;
//                                    }
                                    boolean want2d = processBivariateAttributes && lenArray == 2;
                                    if (name.equals("dmv")) {  // terrible hack
                                        want2d = false;
                                    }
                                    for (int k = 0; k < lenArray; k++) {
                                        if (minValueArray[k] < maxValueArray[k]) {
                                            annotateMedianAttribute(tree, node, name + (k + 1) + "_median", valuesArray[k]);
                                            annotateRangeAttribute(tree, node, name + (k + 1) + "_range", valuesArray[k]);
                                            if (!want2d)
                                                annotateHPDAttribute(tree, node, name + (k + 1) + "_95%_HPD", 0.95, valuesArray[k]);
                                        }
                                    }
                                    // 2D contours
                                    if (want2d) {

                                        boolean variationInFirst = (minValueArray[0] < maxValueArray[0]);
                                        boolean variationInSecond = (minValueArray[1] < maxValueArray[1]);

                                        if (variationInFirst && !variationInSecond)
                                            annotateHPDAttribute(tree, node, name + "1" + "_95%_HPD", 0.95, valuesArray[0]);

                                        if (variationInSecond && !variationInFirst)
                                            annotateHPDAttribute(tree, node, name + "2" + "_95%_HPD", 0.95, valuesArray[1]);

                                        if (variationInFirst && variationInSecond){

                                            for (int l = 0; l < hpd2D.length; l++) {

                                                if (hpd2D[l] > 1) {
                                                    System.err.println("no HPD for proportion > 1 (" + hpd2D[l] + ")");
                                                } else if (hpd2D[l] < 0){
                                                    System.err.println("no HPD for proportion < 0 (" + hpd2D[l] + ")");
                                                }  else {
                                                    annotate2DHPDAttribute(tree, node, name, "_" + (int) (100 * hpd2D[l]) + "%HPD", hpd2D[l], valuesArray);
                                                }

                                           }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                i++;
            }
        }

        private void annotateMeanAttribute(MutableTree tree, NodeRef node, String label, double[] values) {
            double mean = DiscreteStatistics.mean(values);
            tree.setNodeAttribute(node, label, mean);
        }

        private void annotateMedianAttribute(MutableTree tree, NodeRef node, String label, double[] values) {
            double median = DiscreteStatistics.median(values);
            tree.setNodeAttribute(node, label, median);

        }

        private void annotateModeAttribute(MutableTree tree, NodeRef node, String label, HashMap<Object, Integer> values) {
            Object mode = null;
            int maxCount = 0;
            int totalCount = 0;
            int countInMode = 1;

            for (Object key : values.keySet()) {
                int thisCount = values.get(key);
                if (thisCount == maxCount) {
                    // I hope this is the intention
                    mode = mode.toString().concat("+" + key);
                    countInMode++;
                } else if (thisCount > maxCount) {
                    mode = key;
                    maxCount = thisCount;
                    countInMode = 1;
                }
                totalCount += thisCount;
            }
            double freq = (double) maxCount / (double) totalCount * countInMode;
            tree.setNodeAttribute(node, label, mode);
            tree.setNodeAttribute(node, label + ".prob", freq);
        }

        private void annotateFrequencyAttribute(MutableTree tree, NodeRef node, String label, HashMap<Object, Integer> values) {
            double totalCount = 0;
            Set keySet = values.keySet();
            int length = keySet.size();
            String[] name = new String[length];
            Double[] freq = new Double[length];
            int index = 0;
            for (Object key : values.keySet()) {
                name[index] = key.toString();
                freq[index] = new Double(values.get(key));
                totalCount += freq[index];
                index++;
            }
            for (int i = 0; i < length; i++)
                freq[i] /= totalCount;

            tree.setNodeAttribute(node, label + ".set", name);
            tree.setNodeAttribute(node, label + ".set.prob", freq);
        }

        private void annotateRangeAttribute(MutableTree tree, NodeRef node, String label, double[] values) {
            double min = DiscreteStatistics.min(values);
            double max = DiscreteStatistics.max(values);
            tree.setNodeAttribute(node, label, new Object[]{min, max});
        }

        private void annotateHPDAttribute(MutableTree tree, NodeRef node, String label, double hpd, double[] values) {
            int[] indices = new int[values.length];
            HeapSort.sort(values, indices);

            double minRange = Double.MAX_VALUE;
            int hpdIndex = 0;

            int diff = (int) Math.round(hpd * (double) values.length);
            for (int i = 0; i <= (values.length - diff); i++) {
                double minValue = values[indices[i]];
                double maxValue = values[indices[i + diff - 1]];
                double range = Math.abs(maxValue - minValue);
                if (range < minRange) {
                    minRange = range;
                    hpdIndex = i;
                }
            }
            double lower = values[indices[hpdIndex]];
            double upper = values[indices[hpdIndex + diff - 1]];
            tree.setNodeAttribute(node, label, new Object[]{lower, upper});
        }

        // todo Move rEngine to outer class; create once.
        Rengine rEngine = null;

        private final String[] rArgs = {"--no-save"};

//	    private int called = 0;

        private final String[] rBootCommands = {
                "library(MASS)",
                "makeContour = function(var1, var2, prob=0.95, n=50, h=c(1,1)) {" +
                        "post1 = kde2d(var1, var2, n = n, h=h); " +    // This had h=h in argument
                        "dx = diff(post1$x[1:2]); " +
                        "dy = diff(post1$y[1:2]); " +
                        "sz = sort(post1$z); " +
                        "c1 = cumsum(sz) * dx * dy; " +
                        "levels = sapply(prob, function(x) { approx(c1, sz, xout = 1 - x)$y }); " +
                        "line = contourLines(post1$x, post1$y, post1$z, level = levels); " +
                        "return(line) }"
        };

        private String makeRString(double[] values) {
            StringBuffer sb = new StringBuffer("c(");
            sb.append(values[0]);
            for (int i = 1; i < values.length; i++) {
                sb.append(",");
                sb.append(values[i]);
            }
            sb.append(")");
            return sb.toString();
        }

        public static final String CORDINATE = "cordinates";

//		private String formattedLocation(double loc1, double loc2) {
//			return formattedLocation(loc1) + "," + formattedLocation(loc2);
//		}

        private String formattedLocation(double x) {
            return String.format("%5.8f", x);
        }

        private void annotate2DHPDAttribute(MutableTree tree, NodeRef node, String preLabel, String postLabel,
                                            double hpd, double[][] values) {
            int N = 50;
            if (USE_R) {

                // Uses R-Java interface, and the HPD routines from 'emdbook' and 'coda'

                if (rEngine == null) {

                    if (!Rengine.versionCheck()) {
                        throw new RuntimeException("JRI library version mismatch");
                    }

                    rEngine = new Rengine(rArgs, false, null);

                    if (!rEngine.waitForR()) {
                        throw new RuntimeException("Cannot load R");
                    }

                    for (String command : rBootCommands) {
                        rEngine.eval(command);
                    }
                }

                // todo Need a good method to pick grid size


                REXP x = rEngine.eval("makeContour(" +
                        makeRString(values[0]) + "," +
                        makeRString(values[1]) + "," +
                        hpd + "," +
                        N + ")");

                RVector contourList = x.asVector();
                int numberContours = contourList.size();

                if (numberContours > 1) {
                    System.err.println("Warning: a node has a disjoint " + 100 * hpd + "% HPD region.  This may be an artifact!");
                    System.err.println("Try decreasing the enclosed mass or increasing the number of samples.");
                }


                tree.setNodeAttribute(node, preLabel + postLabel + "_modality", numberContours);

                StringBuffer output = new StringBuffer();
                for (int i = 0; i < numberContours; i++) {
                    output.append("\n<" + CORDINATE + ">\n");
                    RVector oneContour = contourList.at(i).asVector();
                    double[] xList = oneContour.at(1).asDoubleArray();
                    double[] yList = oneContour.at(2).asDoubleArray();
                    StringBuffer xString = new StringBuffer("{");
                    StringBuffer yString = new StringBuffer("{");
                    for (int k = 0; k < xList.length; k++) {
                        xString.append(formattedLocation(xList[k])).append(",");
                        yString.append(formattedLocation(yList[k])).append(",");
                    }
                    xString.append(formattedLocation(xList[0])).append("}");
                    yString.append(formattedLocation(yList[0])).append("}");

                    tree.setNodeAttribute(node, preLabel + "1" + postLabel + "_" + (i + 1), xString);
                    tree.setNodeAttribute(node, preLabel + "2" + postLabel + "_" + (i + 1), yString);
                }


            } else { // do not use R


//                KernelDensityEstimator2D kde = new KernelDensityEstimator2D(values[0], values[1], N);
                //ContourMaker kde = new ContourWithSynder(values[0], values[1], N);
                boolean bandwidthLimit = false;

                ContourMaker kde = new ContourWithSynder(values[0], values[1], bandwidthLimit);

                ContourPath[] paths = kde.getContourPaths(hpd);

                tree.setNodeAttribute(node, preLabel + postLabel + "_modality", paths.length);

                if (paths.length > 1) {
                    System.err.println("Warning: a node has a disjoint " + 100 * hpd + "% HPD region.  This may be an artifact!");
                    System.err.println("Try decreasing the enclosed mass or increasing the number of samples.");
                }

                StringBuffer output = new StringBuffer();
                int i = 0;
                for (ContourPath p : paths) {
                    output.append("\n<" + CORDINATE + ">\n");
                    double[] xList = p.getAllX();
                    double[] yList = p.getAllY();
                    StringBuffer xString = new StringBuffer("{");
                    StringBuffer yString = new StringBuffer("{");
                    for (int k = 0; k < xList.length; k++) {
                        xString.append(formattedLocation(xList[k])).append(",");
                        yString.append(formattedLocation(yList[k])).append(",");
                    }
                    xString.append(formattedLocation(xList[0])).append("}");
                    yString.append(formattedLocation(yList[0])).append("}");

                    tree.setNodeAttribute(node, preLabel + "1" + postLabel + "_" + (i + 1), xString);
                    tree.setNodeAttribute(node, preLabel + "2" + postLabel + "_" + (i + 1), yString);
                    i++;

                }
            }
        }

        public BitSet removeClades(Tree tree, NodeRef node, boolean includeTips) {

            BitSet bits = new BitSet();

            if (tree.isExternal(node)) {

                int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
                bits.set(index);

                if (includeTips) {
                    removeClade(bits);
                }

            } else {

                for (int i = 0; i < tree.getChildCount(node); i++) {

                    NodeRef node1 = tree.getChild(node, i);

                    bits.or(removeClades(tree, node1, includeTips));
                }

                removeClade(bits);
            }

            return bits;
        }

        private void removeClade(BitSet bits) {
            Clade clade = cladeMap.get(bits);
            if (clade != null) {
                clade.setCount(clade.getCount() - 1);
            }

        }

        // Get tree clades as bitSets on target taxa
        // codes is an array of existing BitSet objects, which are reused

        void getTreeCladeCodes(Tree tree, BitSet[] codes) {
            getTreeCladeCodes(tree, tree.getRoot(), codes);
        }

        int getTreeCladeCodes(Tree tree, NodeRef node, BitSet[] codes) {
            final int inode = node.getNumber();
            codes[inode].clear();
            if (tree.isExternal(node)) {
                int index = taxonList.getTaxonIndex(tree.getNodeTaxon(node).getId());
                codes[inode].set(index);
            } else {
                for (int i = 0; i < tree.getChildCount(node); i++) {
                    final NodeRef child = tree.getChild(node, i);
                    final int childIndex = getTreeCladeCodes(tree, child, codes);

                    codes[inode].or(codes[childIndex]);
                }
            }
            return inode;
        }

        class Clade {
            public Clade(BitSet bits) {
                this.bits = bits;
                count = 0;
                credibility = 0.0;
            }

            public int getCount() {
                return count;
            }

            public void setCount(int count) {
                this.count = count;
            }

            public double getCredibility() {
                return credibility;
            }

            public void setCredibility(double credibility) {
                this.credibility = credibility;
            }

            public boolean equals(Object o) {
                if (this == o) return true;
                if (o == null || getClass() != o.getClass()) return false;

                final Clade clade = (Clade) o;

                return !(bits != null ? !bits.equals(clade.bits) : clade.bits != null);

            }

            public int hashCode() {
                return (bits != null ? bits.hashCode() : 0);
            }

            public String toString() {
                return "clade " + bits.toString();
            }

            int count;
            double credibility;
            BitSet bits;
            List<Object[]> attributeValues = null;
        }

        //
        // Private stuff
        //
        TaxonList taxonList = null;
        Map<BitSet, Clade> cladeMap = new HashMap<BitSet, Clade>();

        Tree targetTree;
    }

    int totalTrees = 0;
    int totalTreesUsed = 0;
    double posteriorLimit = 0.0;
//PL:    double hpd2D = 0.80;
    double[] hpd2D = {0.80};

    private final List<TreeAnnotationPlugin> plugins = new ArrayList<TreeAnnotationPlugin>();

    Set<String> attributeNames = new HashSet<String>();
    TaxonList taxa = null;

    static boolean processBivariateAttributes = false;

    static {
        try {
            System.loadLibrary("jri");
            processBivariateAttributes = true;
            System.err.println("JRI loaded. Will process bivariate attributes");
        } catch (UnsatisfiedLinkError e) {
//            System.err.print("JRI not available. ");
            if (!USE_R) {
                processBivariateAttributes = true;
//                System.err.println("Using Java bivariate attributes");
            } else {
//                System.err.println("Will not process bivariate attributes");
            }
        }
    }

    public static void printTitle() {
        progressStream.println();
        centreLine("TreeAnnotator " + version.getVersionString() + ", " + version.getDateString(), 60);
        centreLine("MCMC Output analysis", 60);
        centreLine("by", 60);
        centreLine("Andrew Rambaut and Alexei J. Drummond", 60);
        progressStream.println();
        centreLine("Institute of Evolutionary Biology", 60);
        centreLine("University of Edinburgh", 60);
        centreLine("a.rambaut@ed.ac.uk", 60);
        progressStream.println();
        centreLine("Department of Computer Science", 60);
        centreLine("University of Auckland", 60);
        centreLine("alexei@cs.auckland.ac.nz", 60);
        progressStream.println();
        progressStream.println();
    }

    public static void centreLine(String line, int pageWidth) {
        int n = pageWidth - line.length();
        int n1 = n / 2;
        for (int i = 0; i < n1; i++) {
            progressStream.print(" ");
        }
        progressStream.println(line);
    }


    public static void printUsage(Arguments arguments) {

        arguments.printUsage("treeannotator", "<input-file-name> [<output-file-name>]");
        progressStream.println();
        progressStream.println("  Example: treeannotator test.trees out.txt");
        progressStream.println("  Example: treeannotator -burnin 100 -heights mean test.trees out.txt");
        progressStream.println("  Example: treeannotator -burnin 100 -target map.tree test.trees out.txt");
        progressStream.println();
    }

    public static double[] parseVariableLengthDoubleArray(String inString) throws Arguments.ArgumentException {

        List<Double> returnList = new ArrayList<Double>();
        StringTokenizer st = new StringTokenizer(inString,",");
        while(st.hasMoreTokens()) {
            try {
                returnList.add(Double.parseDouble(st.nextToken()));
            } catch (NumberFormatException e) {
                throw new Arguments.ArgumentException();
            }

        }

        if (returnList.size()>0) {
            double[] doubleArray = new double[returnList.size()];
            for(int i=0; i<doubleArray.length; i++)
                doubleArray[i] = returnList.get(i);
            return doubleArray;
        }
        return null;
    }

    //Main method
    public static void main(String[] args) throws IOException {

        // There is a major issue with languages that use the comma as a decimal separator.
        // To ensure compatibility between programs in the package, enforce the US locale.
        Locale.setDefault(Locale.US);

        String targetTreeFileName = null;
        String inputFileName = null;
        String outputFileName = null;

        if (args.length == 0) {
            System.setProperty("com.apple.macos.useScreenMenuBar", "true");
            System.setProperty("apple.laf.useScreenMenuBar", "true");
            System.setProperty("apple.awt.showGrowBox", "true");

            java.net.URL url = LogCombiner.class.getResource("/images/utility.png");
            javax.swing.Icon icon = null;

            if (url != null) {
                icon = new javax.swing.ImageIcon(url);
            }

            final String versionString = version.getVersionString();
            String nameString = "TreeAnnotator " + versionString;
            String aboutString = "<html><center><p>" + versionString + ", " + version.getDateString() + "</p>" +
                    "<p>by<br>" +
                    "Andrew Rambaut and Alexei J. Drummond</p>" +
                    "<p>Institute of Evolutionary Biology, University of Edinburgh<br>" +
                    "<a href=\"mailto:a.rambaut@ed.ac.uk\">a.rambaut@ed.ac.uk</a></p>" +
                    "<p>Department of Computer Science, University of Auckland<br>" +
                    "<a href=\"mailto:alexei@cs.auckland.ac.nz\">alexei@cs.auckland.ac.nz</a></p>" +
                    "<p>Part of the BEAST package:<br>" +
                    "<a href=\"http://beast.community\">http://beast.community</a></p>" +
                    "</center></html>";

            new ConsoleApplication(nameString, aboutString, icon, true);

            // The ConsoleApplication will have overridden System.out so set progressStream
            // to capture the output to the window:
            progressStream = System.out;

            printTitle();

            TreeAnnotatorDialog dialog = new TreeAnnotatorDialog(new JFrame());

            if (!dialog.showDialog("TreeAnnotator " + versionString)) {
                return;
            }

            int burninStates = dialog.getBurninStates();
            int burninTrees = dialog.getBurninTrees();
            double posteriorLimit = dialog.getPosteriorLimit();
            double[] hpd2D = {0.80};
            Target targetOption = dialog.getTargetOption();
            HeightsSummary heightsOption = dialog.getHeightsOption();

            targetTreeFileName = dialog.getTargetFileName();
            if (targetOption == Target.USER_TARGET_TREE && targetTreeFileName == null) {
                System.err.println("No target file specified");
                return;
            }

            inputFileName = dialog.getInputFileName();
            if (inputFileName == null) {
                System.err.println("No input file specified");
                return;
            }

            outputFileName = dialog.getOutputFileName();
            if (outputFileName == null) {
                System.err.println("No output file specified");
                return;
            }

            try {
                new TreeAnnotator(
                        burninTrees,
                        burninStates,
                        heightsOption,
                        posteriorLimit,
                        hpd2D,
                        targetOption,
                        targetTreeFileName,
                        inputFileName,
                        outputFileName);

            } catch (Exception ex) {
                System.err.println("Exception: " + ex.getMessage());
            }

            progressStream.println("Finished - Quit program to exit.");
            while (true) {
                try {
                    Thread.sleep(1000);
                } catch (InterruptedException e) {
                    e.printStackTrace();
                }
            }
        }

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        //new Arguments.StringOption("target", new String[] { "maxclade", "maxtree" }, false, "an option of 'maxclade' or 'maxtree'"),
                        new Arguments.StringOption("heights", new String[]{"keep", "median", "mean", "ca"}, false,
                                "an option of 'keep' (default), 'median', 'mean' or 'ca'"),
                        new Arguments.IntegerOption("burnin", "the number of states to be considered as 'burn-in'"),
                        new Arguments.IntegerOption("burninTrees", "the number of trees to be considered as 'burn-in'"),
                        new Arguments.RealOption("limit", "the minimum posterior probability for a node to be annotated"),
                        new Arguments.StringOption("target", "target_file_name", "specifies a user target tree to be annotated"),
                        new Arguments.Option("help", "option to print this message"),
                        new Arguments.Option("forceDiscrete", "forces integer traits to be treated as discrete traits."),
                        new Arguments.StringOption("hpd2D", "the HPD interval to be used for the bivariate traits", "specifies a (vector of comma seperated) HPD proportion(s)")
                });

        try {
            arguments.parseArguments(args);
        } catch (Arguments.ArgumentException ae) {
            progressStream.println(ae);
            printUsage(arguments);
            System.exit(1);
        }

        if (arguments.hasOption("forceDiscrete")) {
            System.out.println("  Forcing integer traits to be treated as discrete traits.");
            forceIntegerToDiscrete = true;
        }

        if (arguments.hasOption("help")) {
            printUsage(arguments);
            System.exit(0);
        }

        HeightsSummary heights = HeightsSummary.KEEP_HEIGHTS;
        if (arguments.hasOption("heights")) {
            String value = arguments.getStringOption("heights");
            if (value.equalsIgnoreCase("mean")) {
                heights = HeightsSummary.MEAN_HEIGHTS;
            } else if (value.equalsIgnoreCase("median")) {
                heights = HeightsSummary.MEDIAN_HEIGHTS;
            } else if (value.equalsIgnoreCase("ca")) {
                heights = HeightsSummary.CA_HEIGHTS;
                System.out.println("Please cite: Heled and Bouckaert: Looking for trees in the forest:\n" +
                        "summary tree from posterior samples. BMC Evolutionary Biology 2013 13:221.");
            }
        }

        int burninStates = -1;
        int burninTrees = -1;
        if (arguments.hasOption("burnin")) {
            burninStates = arguments.getIntegerOption("burnin");
        }
        if (arguments.hasOption("burninTrees")) {
            burninTrees = arguments.getIntegerOption("burninTrees");
        }

        double posteriorLimit = 0.0;
        if (arguments.hasOption("limit")) {
            posteriorLimit = arguments.getRealOption("limit");
        }

        double[] hpd2D = {80};
        if (arguments.hasOption("hpd2D")) {
            try {
                hpd2D = parseVariableLengthDoubleArray(arguments.getStringOption("hpd2D"));
            } catch (Arguments.ArgumentException e) {
                System.err.println("Error reading " + arguments.getStringOption("hpd2D"));
            }
        }

        Target target = Target.MAX_CLADE_CREDIBILITY;
        if (arguments.hasOption("target")) {
            target = Target.USER_TARGET_TREE;
            targetTreeFileName = arguments.getStringOption("target");
        }

        final String[] args2 = arguments.getLeftoverArguments();

        switch (args2.length) {
            case 2:
                outputFileName = args2[1];
                // fall to
            case 1:
                inputFileName = args2[0];
                break;
            default: {
                System.err.println("Unknown option: " + args2[2]);
                System.err.println();
                printUsage(arguments);
                System.exit(1);
            }
        }

        new TreeAnnotator(burninTrees, burninStates, heights, posteriorLimit, hpd2D, target, targetTreeFileName, inputFileName, outputFileName);

        System.exit(0);
    }

    /**
     * @author Andrew Rambaut
     * @version $Id$
     */
    public static interface TreeAnnotationPlugin {
        Set<String> setAttributeNames(Set<String> attributeNames);

        boolean handleAttribute(Tree tree, NodeRef node, String attributeName, double[] values);
    }

    // very inefficient, but Java wonderful bitset has no subset op
    // perhaps using bit iterator would be faster, I can't br bothered.

    static boolean isSubSet(BitSet x, BitSet y) {
        y = (BitSet) y.clone();
        y.and(x);
        return y.equals(x);
    }

    boolean setTreeHeightsByCA(MutableTree targetTree, final String inputFileName, final int burnin)
            throws IOException, Importer.ImportException {
        progressStream.println("Setting node heights...");
        progressStream.println("0              25             50             75            100");
        progressStream.println("|--------------|--------------|--------------|--------------|");

        int reportStepSize = totalTrees / 60;
        if (reportStepSize < 1) reportStepSize = 1;

        final FileReader fileReader = new FileReader(inputFileName);
        final NexusImporter importer = new NexusImporter(fileReader, true);

        // this call increments the clade counts and it shouldn't
        // this is remedied with removeClades call after while loop below
        CladeSystem cladeSystem = new CladeSystem(targetTree);
        final int nClades = cladeSystem.getCladeMap().size();

        // allocate posterior tree nodes order once
        int[] postOrderList = new int[nClades];
        BitSet[] ctarget = new BitSet[nClades];
        BitSet[] ctree = new BitSet[nClades];

        for (int k = 0; k < nClades; ++k) {
            ctarget[k] = new BitSet();
            ctree[k] = new BitSet();
        }

        cladeSystem.getTreeCladeCodes(targetTree, ctarget);

        // temp collecting heights inside loop allocated once
        double[] hs = new double[nClades];

        // heights total sum from posterior trees
        double[] ths = new double[nClades];

        totalTreesUsed = 0;

        int counter = 0;
        while (importer.hasTree()) {
            final Tree tree = importer.importNextTree();

            if (counter >= burnin) {
                TreeUtils.preOrderTraversalList(tree, postOrderList);
                cladeSystem.getTreeCladeCodes(tree, ctree);
                for (int k = 0; k < nClades; ++k) {
                    int j = postOrderList[k];
                    for (int i = 0; i < nClades; ++i) {
                        if( isSubSet(ctarget[i], ctree[j]) ) {
                            hs[i] = tree.getNodeHeight(tree.getNode(j));
                        }
                    }
                }
                for (int k = 0; k < nClades; ++k) {
                    ths[k] += hs[k];
                }
                totalTreesUsed += 1;
            }
            if (counter > 0 && counter % reportStepSize == 0) {
                progressStream.print("*");
                progressStream.flush();
            }
            counter++;

        }
        cladeSystem.removeClades(targetTree, targetTree.getRoot(), true);
        for (int k = 0; k < nClades; ++k) {
            ths[k] /= totalTreesUsed;
            final NodeRef node = targetTree.getNode(k);
            targetTree.setNodeHeight(node, ths[k]);
        }
        fileReader.close();

        progressStream.println();
        progressStream.println();

        return true;
    }
}


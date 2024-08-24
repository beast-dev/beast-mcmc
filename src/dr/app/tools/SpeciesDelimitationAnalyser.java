/*
 * SpeciesDelimitationAnalyser.java
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

package dr.app.tools;

import dr.app.beast.BeastVersion;
import dr.app.util.Arguments;
import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evomodel.alloppnet.speciation.BirthDeathCollapseModel;
import dr.util.Version;

import java.io.*;
import java.util.*;

/**
 * @author Graham Jones
 *         Date: 01/09/2013
 */
public class SpeciesDelimitationAnalyser {

    private final static Version version = new BeastVersion();
    private ArrayList<Cluster> clusters;
    private TaxonList taxonList;

    private int burnin;
    double collapseheight;
    private double  similaritycutoff;
    private String inputFileName;
    private String outputFileName;

    private class Cluster {
        private int count;
        private int [] partition;
        private double totalsimilarity;
        private boolean deleted;

        public Cluster(TaxonList taxonList, Tree tree, double collapseheight) {
            int nshort = 0;
            count = 1;
            partition = new int[taxonList.getTaxonCount()];
            totalsimilarity = 1.0;
            deleted = false;

            int nnodes = tree.getNodeCount();
            for (int n = 0; n < nnodes; n++) {
                NodeRef nr = tree.getNode(n);
                boolean collapse;
                if (tree.isRoot(nr)) {
                    collapse = BirthDeathCollapseModel.belowCollapseHeight(tree.getNodeHeight(nr), collapseheight);
                } else {
                    NodeRef anc  = tree.getParent(nr);
                    collapse = (BirthDeathCollapseModel.belowCollapseHeight(tree.getNodeHeight(nr), collapseheight)
                            && !BirthDeathCollapseModel.belowCollapseHeight(tree.getNodeHeight(anc), collapseheight));
                }
                if (collapse) {
                    nshort++;
                    nodeToClade(taxonList, tree, nr, nshort);
                }
            }
        }


        double distance(Cluster sd) {
            double d = 0.0;
            for (int i = 0; i < partition.length; i++) {
                for (int j = 0; j < i; j ++) {
                    boolean pair = partition[i] == partition[j];
                    boolean sdpair = sd.partition[i] == sd.partition[j];
                    if (pair != sdpair) {
                        d += 1.0;
                    }
                }
            }
            return d / (partition.length * (partition.length-1.0) / 2.0);
        }


        double similarity(Cluster sd) {
            return 1.0 - distance(sd);
        }



        private void nodeToClade(TaxonList taxonlist, Tree tree, NodeRef nr, int idx) {
            if (tree.isExternal(nr)) {
                String id = tree.getNodeTaxon(nr).getId();
                int i = taxonIdToInt(taxonlist, id);
                partition[i] = idx;
            } else {
                assert 2 == tree.getChildCount(nr);
                nodeToClade(taxonlist, tree, tree.getChild(nr, 0), idx);
                nodeToClade(taxonlist, tree, tree.getChild(nr, 1), idx);
            }
        }

        int taxonIdToInt(TaxonList taxonList, String id) {
            int j = -1;
            for (int i = 0; i < taxonList.getTaxonCount(); i++) {
                String tid = taxonList.getTaxonId(i);
                if (id.compareTo(tid) == 0) {
                    assert j < 0;
                    j = i;
                }
            }
            assert j >= 0;
            return j;
        }
    }


    SpeciesDelimitationAnalyser(int burnin, double collapseheight, double  similaritycutoff, String inputFileName, String outputFileName) {
        this.burnin = burnin;
        this.collapseheight = collapseheight;
        this.similaritycutoff = similaritycutoff;
        this.inputFileName = inputFileName;
        this.outputFileName = outputFileName;
    }

    void readtrees() throws IOException {
        clusters = new ArrayList<Cluster>(0);
        int totalTrees;
        FileReader fileReader = new FileReader(inputFileName);
        NexusImporter importer = new NexusImporter(fileReader);
        System.out.println("Reading trees...");
        try {
            taxonList = importer.parseTaxaBlock();
            totalTrees = 0;
            while (importer.hasTree()) {
                Tree tree = importer.importNextTree();

                if (totalTrees >= burnin) {
                    Cluster sd = new Cluster(taxonList, tree, collapseheight);
                    clusters.add(sd);
                }

                if (totalTrees > 0 && (totalTrees % 100 == 0)) {
                    System.out.print("*");
                    System.out.flush();
                }
                if (totalTrees > 0 && (totalTrees % 5000 == 0)) {
                    System.out.println(" " + totalTrees);
                    System.out.flush();
                }
                totalTrees++;
            }
            System.out.println("");
        } catch (Importer.ImportException e) {
            System.err.println("Error Parsing Input Tree: " + e.getMessage());
            return;
        }
        fileReader.close();
        if (totalTrees < 1) {
            System.err.println("No trees");
            return;
        }
        System.out.println("Total trees read: " + totalTrees);
    }



    void countclusterings() {
        if (burnin > 0) {
            System.out.println("Ignoring first " + burnin + " trees.");
        }
        System.out.println("Counting clusterings... ");
        for (int i = 0; i < clusters.size(); i++) {
            if (i > 0 && (i % 100 == 0)) {
                System.out.print("*");
                System.out.flush();
            }
            if (i > 0 && (i % 5000 == 0)) {
                System.out.println(" " + i);
                System.out.flush();
            }

            for (int j = i+1; j < clusters.size(); j++) {
                if (clusters.get(j).count > 0) {
                    double simij = clusters.get(i).similarity(clusters.get(j));
                    if (simij == 1.0) {
                        clusters.get(i).count += 1;
                        clusters.get(j).count = 0;
                        clusters.get(j).deleted = true;
                    }
                    if (simij >= similaritycutoff  ||  simij == 1.0) {
                        clusters.get(i).totalsimilarity += simij;
                    }
                }
            }
        }
        Collections.sort(clusters, CLUSTER_COMPARATOR);
    }



    void writeresults() throws IOException {
        FileWriter fileWriter = new FileWriter(outputFileName);

        fileWriter.write("count      fraction               similarity             nclusters  ");
        for (int i = 0; i < taxonList.getTaxonCount(); i++) {
            String tstr = taxonList.getTaxonId(i);
            while (tstr.length() < 4) tstr += " ";
            fileWriter.write(tstr + " ");
        }
        fileWriter.write("\n");
        for (int i = 0; i < clusters.size(); i++) {
            Cluster sd = clusters.get(i);
            if (sd.count > 0) {
                String countstr = "" + sd.count;
                while (countstr.length() < 10) countstr += " ";
                String fracstr = "" + (double)sd.count / (double) clusters.size();
                while (fracstr.length() < 22) fracstr += " ";
                String simstr = "" + sd.totalsimilarity;
                while (simstr.length() < 22) simstr += " ";
                int clustercount = 0;
                for (int j = 0; j < sd.partition.length; j ++) {
                    if (sd.partition[j] > clustercount) {
                        clustercount = sd.partition[j];
                    }
                }
                String ccstr = "" + clustercount;
                while (ccstr.length() < 10) ccstr += " ";
                fileWriter.write(countstr + " " + fracstr + " " + simstr + " " + ccstr + " ");
                for (int j = 0; j < sd.partition.length; j ++) {
                    String pjstr = "" + sd.partition[j];
                    while (pjstr.length() < 4) { pjstr += " "; }
                    fileWriter.write(pjstr + " ");
                }
                fileWriter.write("\n");
            }
        }
        fileWriter.write("\n");
        fileWriter.close();

    }



    public static void printTitle() {
        System.out.println();
        centreLine("SpeciesDelimitationAnalyser " + version.getVersionString() + ", " + version.getDateString(), 60);
        centreLine("Finds clusterings of individuals into clusters", 60);
        centreLine("(ie possible species) from MCMC tree samples", 60);
        centreLine("by", 60);
        centreLine("Graham Jones", 60);
        centreLine("www.indriid.com", 60);
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

        arguments.printUsage("SpeciesDelimitationAnalyser", "<input-file-name> <output-file-name>");
        System.out.println();
        System.out.println("  Example: SpeciesDelimitationAnalyser treesamples.txt out.txt");
        System.out.println();
    }


    static final Comparator<Cluster> CLUSTER_COMPARATOR = new Comparator<Cluster>() {
        public int compare(Cluster a, Cluster b) {
            if (a.deleted != b.deleted) {
                return a.deleted ? 1 : -1;
            }
            if (b.totalsimilarity != a.totalsimilarity) {
                return (b.totalsimilarity > a.totalsimilarity) ? 1 : -1;
            }
            return 0;
        }
    };




    public static void main(String[] args) throws java.io.IOException {

        Locale.setDefault(Locale.US);

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.IntegerOption("burnin", "the number of states to be considered as 'burn-in' [default = none]"),
                        new Arguments.RealOption("collapseheight", "the height below which nodes get collapsed [default = .001]"),
                        new Arguments.RealOption("simcutoff", "the value above which two clusters are regarded as similar enough to support one another's credibility [default = .9]"),
                        new Arguments.Option("help", "option to print this message")
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

        int burnin = -1;
        if (arguments.hasOption("burnin")) {
            burnin = arguments.getIntegerOption("burnin");
        }
        System.out.println();

        double collapseheight = 0.001;
        if (arguments.hasOption("collapseheight")) {
            collapseheight = arguments.getRealOption("collapseheight");
        }

        double similaritycutoff = 0.9;
        if (arguments.hasOption("simcutoff")) {
            similaritycutoff = arguments.getRealOption("simcutoff");
        }

        System.out.println("burnin " + burnin + " collapseheight " + collapseheight + " simcutoff "+ similaritycutoff);

        String[] args2 = arguments.getLeftoverArguments();
        if (args2.length > 2) {
            System.err.println("Unknown option: " + args2[2]);
            System.err.println();
            printUsage(arguments);
            System.exit(1);
        }
        if (args2.length < 2) {
            System.err.println("Input filename and outputfilename required");
            System.err.println();
            printUsage(arguments);
            System.exit(1);
        }

        String inputFileName = args2[0];
        String outputFileName = args2[1];

        SpeciesDelimitationAnalyser spDA =
                new SpeciesDelimitationAnalyser(burnin, collapseheight, similaritycutoff, inputFileName, outputFileName);
        spDA.readtrees();
        spDA.countclusterings();
        spDA.writeresults();

        System.exit(0);

    }
}

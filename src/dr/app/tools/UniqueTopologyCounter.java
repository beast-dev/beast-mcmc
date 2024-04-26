/*
 * UniqueTopologyCounter.java
 *
 * Copyright (c) 2002-2022 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import dr.app.util.Arguments;
import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.TreeImporter;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Locale;

/**
 * UniqueTopologyCounter processes a BEAST .trees file and counts the number of unique topologies
 * Branch lengths are not taken into account and a number of trees can be discarded as burn-in
 *
 * @author Guy Baele
 */

public class UniqueTopologyCounter {

    private long totalTrees;
    private HashMap<String, Integer> treeCounter;

    public UniqueTopologyCounter(String inputFileName, int burninTrees) {
        this.totalTrees = 0;
        this.treeCounter = new HashMap<String, Integer>();

        int counter = 0;

        FileReader fileReader = null;
        TreeImporter importer = null;
        try {
            fileReader = new FileReader(inputFileName);
            importer = new NexusImporter(fileReader, true);

            while (importer.hasTree()) {
                Tree tree = importer.importNextTree();
                String treeString = TreeUtils.uniqueNewick(tree, tree.getRoot());
                //String treeString = TreeUtils.newickNoLengths(tree);
                counter++;
                if (counter > burninTrees) {
                    if (treeCounter.containsKey(treeString)) {
                        int count = treeCounter.remove(treeString);
                        treeCounter.put(treeString, count + 1);
                    } else {
                        treeCounter.put(treeString, 1);
                    }
                    this.totalTrees++;
                }
            }
            //print out all trees and their counts
            System.out.println(this.totalTrees + " trees read from file " + inputFileName);
            System.out.println(counter + " trees processed after removing " + burninTrees + " as burn-in");
            System.out.println(treeCounter.size() + " unique trees found");
            for (String s : treeCounter.keySet()) {
                System.out.println(s + " : " + treeCounter.get(s));
            }
        } catch (FileNotFoundException fnfe) {
            System.err.println("burninTrees = " + burninTrees);
            fnfe.printStackTrace();
        } catch (IOException io) {
            io.printStackTrace();
        } catch (Importer.ImportException e) {
            System.err.println("Error Parsing Input Tree: " + e.getMessage());
        }
    }

    public static void printUsage(Arguments arguments) {
        arguments.printUsage("UniqueTopologyCounter", "<input-file-name>");
        System.out.println();
        System.out.println("  Example: UniqueTopologyCounter -burninTrees 1000 test.trees");
        System.out.println();
    }

    public static void main(String[] args) {

        // There is a major issue with languages that use the comma as a decimal separator.
        // To ensure compatibility between programs in the package, enforce the US locale.
        Locale.setDefault(Locale.US);

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.IntegerOption("burninTrees", "the number of trees to be considered as 'burn-in'"),
                        new Arguments.Option("help", "option to print this message")
                });

        try {
            arguments.parseArguments(args);
        } catch (Arguments.ArgumentException ae) {
            System.out.println(ae.getMessage());
            printUsage(arguments);
            System.exit(1);
        }

        int burninTrees = -1;
        if (arguments.hasOption("burninTrees")) {
            burninTrees = arguments.getIntegerOption("burninTrees");
        }

        String[] args2 = arguments.getLeftoverArguments();

        if (args2.length == 1) {
            String inputFileName = args2[0];
            new UniqueTopologyCounter(inputFileName, burninTrees);
        } else {
            throw new RuntimeException("Incorrect number of arguments.");
        }

        System.exit(0);

    }

}

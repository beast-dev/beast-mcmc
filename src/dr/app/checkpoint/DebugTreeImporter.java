/*
 * DebugTreeImporter.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.app.checkpoint;

import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.Tree;

import java.io.*;

/**
 * @author Guy Baele
 */
public class DebugTreeImporter {

    private static final boolean DEBUG = true;

    //just trolling
    public static boolean isInteger(String s) {
        try {
            Integer.parseInt(s);
        } catch(NumberFormatException e) {
            return false;
        } catch(NullPointerException e) {
            return false;
        }
        // only got here if we didn't return false
        return true;
    }

    //args[0]: saved dump file
    //args[1]: new XML file
    //args[2]: new tree (Newick format)
    //args[3]: output file
    //e.g.: dump.trump.160 worobey-con-sc-resume-162-taxa.xml newtree.nwk restart.dump
    public static void main (String[] args) {

        if (args.length != 4) {
            throw new RuntimeException("Incorrect number of arguments.");
        }

        try {
            FileReader fileIn = new FileReader(args[2]);
            BufferedReader in = new BufferedReader(fileIn);

            String fullTree = in.readLine();

            NewickImporter importer = new NewickImporter(fullTree);
            Tree tree = importer.importNextTree();

            //assume rootHeight in dumpfile is ALWAYS followed by all the nodeHeights
            System.out.println("root height: " + tree.getNodeHeight(tree.getRoot()));

            FileReader dumpFileIn = new FileReader(args[0]);
            BufferedReader dumpIn = new BufferedReader(dumpFileIn);

            FileWriter dumpFileOut = new FileWriter(args[3]);
            BufferedWriter dumpOut = new BufferedWriter(dumpFileOut);

            String original = dumpIn.readLine();
            String[] fields = original.split("\t");

            //write new rootHeight to new output dump file
            while (!fields[0].equals("treeModel.rootHeight")) {
                if (DEBUG) {
                    System.out.println(original);
                }
                dumpOut.write(original + "\n");

                original = dumpIn.readLine();
                fields = original.split("\t");
            }

            fields[2] = Double.toString(tree.getNodeHeight(tree.getRoot()));
            for (int i = 0; i < fields.length; i++) {
                dumpOut.write(fields[i] + "\t");
            }
            dumpOut.write("\n");

            //write all new node heights to new output dump file
            int nodeCount = tree.getNodeCount();
            if (DEBUG) {
                System.out.println(nodeCount + " nodes found in tree.");
            }
            for (int i = 0; i < (nodeCount-1); i++) {
                if (DEBUG) {
                    System.out.println(tree.getNode(i).getNumber() + "\t" + tree.getNodeHeight(tree.getNode(i)));
                }
                dumpOut.write(tree.getNode(i).getNumber() + "\t1\t" + tree.getNodeHeight(tree.getNode(i)) + "\n");
            }

            //skip all the node heights in the original dump file
            //no clue as to how many there are ...
            //best I can tell only node heights have integers as parameter names
            original = dumpIn.readLine();
            if (DEBUG) {
                System.out.println(original);
            }
            fields = original.split("\t");
            while (isInteger(fields[0])) {
                original = dumpIn.readLine();
                fields = original.split("\t");
            }

            while (!fields[0].equals("treeModel")) {
                dumpOut.write(original + "\n");
                original = dumpIn.readLine();
                if (DEBUG) {
                    System.out.println(original);
                }
                fields = original.split("\t");
            }

            dumpOut.write(fields[0] + "\t" + fullTree);

            dumpOut.flush();
            dumpOut.close();

        } catch (FileNotFoundException fnfe) {
            throw new RuntimeException("Tree file not found.");
        } catch (IOException ioe) {
            throw new RuntimeException("Unable to read file: " + ioe.getMessage());
        } catch (Importer.ImportException ie) {
            throw new RuntimeException("Unable to import tree: " + ie.getMessage());
        }
    }

}

/*
 * SequentialTreeReader.java
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

import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.TreeImporter;
import dr.evolution.tree.Tree;

import java.io.FileReader;
import java.io.IOException;


public class SequentialTreeReader {

    private static final Boolean DEBUG = false;

    private final String inputFileName;
    private final int burnin;
    private static BaseTreeTool.TreeProgressPrinter progressPrinter = new BaseTreeTool.TreeProgressPrinter(System.err);
    private int totalTrees;
    private int totalUsedTrees;

    private Tree currentTree;
    private int currentIndex;

    private final FileReader reader;
    private final TreeImporter importer;


    SequentialTreeReader(String inputFileName, int burnin) throws IOException {
        this.inputFileName = inputFileName;
        this.burnin = burnin;
        this.reader = new FileReader(inputFileName);
        this.importer = new NexusImporter(reader, false);

    }


    public Tree getTree(int index) throws IOException {

        if (DEBUG) {
            System.out.println(index);
        }

        if (index < burnin) {
            throw new RuntimeException("Tree " + index + " was requested, but burnin is set at " + burnin + ".");
        }

        try {
            if (currentTree == null) {
                progressPrinter.printReadingTrees();
                currentTree = importer.importNextTree();
                currentIndex = 0;
                totalTrees = 1;
                totalUsedTrees = 1;
            }

            if (index < currentIndex) {
                throw new RuntimeException("The last tree accessed was " + currentIndex + ", but " + index +
                        " was requested. Cannot go backwards. (for developers: consider using BaseTreeTools.readTrees)");
            } else if (index == currentIndex) {
                return currentTree;
            } else {
                for (int i = currentIndex; i < index; i++) {


                    if (importer.hasTree()) {
                        progressPrinter.printProgress(totalTrees);
                        currentTree = importer.importNextTree();
                        totalTrees++;
                    } else {
                        reader.close();
                        progressPrinter.printSummary(totalTrees, totalUsedTrees, burnin);
                        return null;
                    }
                }

                totalUsedTrees++;
                currentIndex = index;
                return currentTree;
            }

        } catch (Importer.ImportException e) {
            reader.close();
            System.err.println("Error Parsing Input Tree: " + e.getMessage());
            return null;
        }
    }
}

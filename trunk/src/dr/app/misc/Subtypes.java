/*
 * Subtypes.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.app.misc;

import dr.app.util.Utils;
import dr.evolution.alignment.Alignment;
import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.evolution.parsimony.ParsimonyCriterion;
import dr.evolution.parsimony.SankoffParsimony;

import java.io.EOFException;
import java.io.FileReader;
import java.io.IOException;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 *
 * @version $Id: Subtypes.java,v 1.1 2005/06/28 11:51:50 rambaut Exp $
 */
public class Subtypes {

    public Subtypes(String fileName, String treeFileName) throws IOException {

        Alignment alignment = null;
        Tree tree = null;

        try {
            alignment = readNexusFile(fileName);
            tree = readTreeFile(treeFileName);

        } catch (IOException ioex) {
            System.err.println("File I/O Error: " + ioex);
            return;
        } catch (Exception ex) {
            System.err.println("Fatal exception: " + ex);
            return;
        }

        if (alignment == null) {
            System.err.println("No alignment in file: " + fileName);
            return;
        }

        if (tree == null) {
            System.err.println("No tree in file: " + fileName);
            return;
        }

        System.out.println("Alignment read: " + alignment.getSequenceCount() + " taxa, " + alignment.getSiteCount() + " sites.");
        System.out.println("Tree read: " + tree.getTaxonCount() + " taxa.");

        ParsimonyCriterion parsimony = new SankoffParsimony(alignment);
        
    }

    /**
     * @param fileName
     * @throws java.io.IOException
     */
    private Alignment readNexusFile(String fileName) throws IOException {

        Alignment alignment = null;
        TaxonList taxonList = null;

        try {
            FileReader reader = new FileReader(fileName);

            NexusImporter importer = new NexusImporter(reader);

            boolean done = false;

            while (!done) {
                try {

                    NexusImporter.NexusBlock block = importer.findNextBlock();

                    if (block == NexusImporter.TAXA_BLOCK) {

                        if (taxonList != null) {
                            throw new NexusImporter.MissingBlockException("TAXA block already defined");
                        }

                        taxonList = importer.parseTaxaBlock();

                    } else if (block == NexusImporter.DATA_BLOCK) {

                        // A data block doesn't need a taxon block before it
                        // but if one exists then it will use it.
                        alignment = importer.parseDataBlock(taxonList);
                        if (taxonList == null) {
                            taxonList = alignment;
                        }

                    } else if (block == NexusImporter.TREES_BLOCK) {

                        // ignore tree block
                    } else {
                        // Ignore the block..
                    }

                } catch (EOFException ex) {
                    done = true;
                }
            }

        } catch (Importer.ImportException ime) {
            System.err.println("Error reading alignment: " + ime);
        }

        return alignment;
    }

    /**
     * @param fileName
     * @throws java.io.IOException
     */
    private Tree readTreeFile(String fileName) throws IOException {

        Alignment alignment = null;
        Tree[] trees = null;
        TaxonList taxonList = null;

        try {
            FileReader reader = new FileReader(fileName);

            NexusImporter importer = new NexusImporter(reader);

            boolean done = false;

            while (!done) {
                try {

                    NexusImporter.NexusBlock block = importer.findNextBlock();

                    if (block == NexusImporter.TAXA_BLOCK) {

                        if (taxonList != null) {
                            throw new NexusImporter.MissingBlockException("TAXA block already defined");
                        }

                        taxonList = importer.parseTaxaBlock();

                    } else if (block == NexusImporter.DATA_BLOCK) {

                        // A data block doesn't need a taxon block before it
                        // but if one exists then it will use it.
                        alignment = importer.parseDataBlock(taxonList);
                        if (taxonList == null) {
                            taxonList = alignment;
                        }

                    } else if (block == NexusImporter.TREES_BLOCK) {

                        trees = importer.parseTreesBlock(taxonList);
                        if (taxonList == null) {
                            taxonList = alignment;
                        }

                    } else {
                        // Ignore the block..
                    }

                } catch (EOFException ex) {
                    done = true;
                }
            }

        } catch (Importer.ImportException ime) {
            System.err.println("Error reading alignment: " + ime);
        }

        return trees[0];
    }

    // Main entry point
    static public void main(String[] args) throws IOException {
        if (args.length > 2) {
            System.err.println("Unknown option: " + args[2]);
            System.err.println();
            System.out.println("Usage: subtype <alignmentfilename> <treefilename>");
            System.exit(1);
        }

        String inputFileName = null;
        String treeFileName = null;

        if (args.length > 0) {
            inputFileName = args[0];
        }

        if (args.length > 1) {
            treeFileName = args[1];
        }


        if (inputFileName == null) {
            // No input file name was given so throw up a dialog box...
            inputFileName = Utils.getLoadFileName("Select NEXUS file to analyse");
        }

        if (treeFileName == null) {
            // No input file name was given so throw up a dialog box...
            treeFileName = Utils.getLoadFileName("Select tree file to analyse");
        }

        new Subtypes(inputFileName, treeFileName);
    }

}

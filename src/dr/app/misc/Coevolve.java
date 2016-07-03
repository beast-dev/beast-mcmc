/*
 * Coevolve.java
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

package dr.app.misc;

import dr.app.util.Utils;
import dr.evolution.alignment.*;
import dr.evolution.datatype.DataType;
import dr.evolution.datatype.GeneralDataType;
import dr.evolution.io.Importer;
import dr.evolution.io.NexusImporter;
import dr.evolution.sequence.Sequence;
import dr.evolution.tree.Tree;
import dr.evolution.util.TaxonList;
import dr.oldevomodel.sitemodel.GammaSiteModel;
import dr.oldevomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.*;
import dr.evomodel.tree.TreeModel;
import dr.oldevomodel.substmodel.FrequencyModel;
import dr.oldevomodel.treelikelihood.TreeLikelihood;
import dr.inference.model.Parameter;
import dr.math.*;
import dr.oldevomodel.substmodel.GeneralSubstitutionModel;
import dr.oldevomodel.substmodel.SubstitutionModel;

import java.io.*;

/**
 * @author Alexei Drummond
 * @author Andrew Rambaut
 *
 * @version $Id: Coevolve.java,v 1.3 2006/01/09 15:12:17 rambaut Exp $
 */
public class Coevolve {

    public Coevolve(String fileName, String treeFileName) throws java.io.IOException {

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

        PrintWriter pw = new PrintWriter(new FileWriter(fileName + ".out"));
        for (int dist = 3; dist < Math.min(alignment.getSiteCount()-200,2000); dist += 3) {

            System.out.print(dist+"\t");
            pw.write(dist+"\t");
            // create a pair alignment of only third positions
            Alignment pairAlignment = makePairAlignment(alignment, dist,2,3);
            SitePatterns sitePatterns = new SitePatterns(pairAlignment);

            Parameter beta = new Parameter.Default(2, 1.0);
            beta.setId("betakappa");
            SubstitutionModel pairModel = makePairSubstModel(alignment, pairAlignment, beta);
            runModel2(sitePatterns, pw, tree, pairModel, beta);
        }
    }

    private SubstitutionModel makePairSubstModel(Alignment alignment, Alignment pairAlignment, Parameter beta) {

        DataType newDataType = makePairDataType(alignment);

        double[] frequencies = pairAlignment.getStateFrequencies();
        Parameter frequencyParameter = new Parameter.Default(frequencies);
        FrequencyModel freqModel = new FrequencyModel(newDataType, frequencyParameter);

        return new SitePairSubstitutionModel(newDataType, freqModel, beta);
    }

    class SitePairSubstitutionModel extends GeneralSubstitutionModel {

        public SitePairSubstitutionModel(DataType dataType, FrequencyModel freqModel, Parameter beta) {
            super(dataType, freqModel, beta, 0);
        }

        protected void setupRelativeRates() {

            // WARNING!! Assumes nucleotides!

            int stateCount = dataType.getStateCount();

            int k = 0;
            for (int i = 0; i < stateCount; i++) {
                for (int j = i+1; j < stateCount; j++) {
                    relativeRates[k] = 1.0;
                    int x1 = i % 4;
                    int y1 = i / 4;
                    int x2 = j % 4;
                    int y2 = j / 4;
                    // if both nucleotides different in the two pairs then use beta
                    if (x1 != x2 && y1 != y2) {
                        relativeRates[k] = ratesParameter.getParameterValue(0);
                    }

                    // if transition at first position then multiply rate by kappa
                    if (x1 != x2 && (x1 % 2 == x2 % 2)) {
                        relativeRates[k] *= ratesParameter.getParameterValue(1);
                    }
                    // if transition at second position then multiply rate by kappa
                    if (y1 != y2 && (y1 % 2 == y2 % 2)) {
                        relativeRates[k] *= ratesParameter.getParameterValue(1);
                    }

                    k += 1;
                }
            }
        }
    }

    private DataType makePairDataType(Alignment alignment) {

        DataType dataType = alignment.getDataType();
        int stateCount = dataType.getStateCount();

        String[] stateChars = new String[stateCount*stateCount];
        for (int i = 0; i < Math.min(10, stateChars.length); i++) { stateChars[i] = String.valueOf((char)('0' + i)); }
        if (stateCount > 10) {
            for (int i = 10; i < stateChars.length; i++) { stateChars[i] = String.valueOf((char)('A' + i-10)); }
        }

        GeneralDataType newDataType = new GeneralDataType(stateChars);
//        int[] ambiguities = new int[stateChars.length];
//        for (int i = 0; i < ambiguities.length; i++) {
//            ambiguities[i] = i;
//        }
//
//        // question mark for unknown states
//        newDataType.addAmbiguity("?", ambiguities);

        return newDataType;
    }

    /**
     * @param alignment
     * @param dist the distance |i-j| between the site pairs (j = i + dist).
     * @param start the first i value considered
     * @param step the distance between successive i values (i += step)
     * @return an alignment of sites, where each site corresponds to a pair (i, j) in the old alignment.
     */
    private Alignment makePairAlignment(Alignment alignment, int dist, int start, int step) {

        DataType dataType = alignment.getDataType();
        int stateCount = dataType.getStateCount();

        DataType newDataType = makePairDataType(alignment);

        SimpleAlignment pairAlignment = new SimpleAlignment();
        pairAlignment.setDataType(newDataType);
        for (int k = 0; k < alignment.getSequenceCount(); k++) {
            StringBuffer sequence = new StringBuffer();
            for (int i = start; (i+dist) < alignment.getSiteCount(); i+= step) {
                int j = i + dist;

                int state = alignment.getState(k,i) * stateCount + alignment.getState(k, j);
                if (state < newDataType.getStateCount()) {
                    sequence.append(newDataType.getChar(state));
                } else {
                    sequence.append('?');
                }
            }
            pairAlignment.addSequence(new Sequence(alignment.getTaxon(k),sequence.toString()));
        }
        return pairAlignment;
    }

    private void runModel2(PatternList patternList, PrintWriter pw, Tree tree, SubstitutionModel substModel, final Parameter betaParameter) {

        final Parameter muParameter = new Parameter.Default(1.0);
        muParameter.setId("mu");

        SiteModel siteModel = new GammaSiteModel(substModel, muParameter, null, 1, null);
        TreeModel treeModel = new TreeModel(tree);
        final TreeLikelihood treeLikelihood = new TreeLikelihood(patternList, treeModel, siteModel, null, null, false, false, true, false, false);
        treeLikelihood.setId("likelihood");

        MultivariateFunction function = new MultivariateFunction() {

            public double evaluate(double[] argument) {
                betaParameter.setParameterValue(0, argument[0]);
                betaParameter.setParameterValue(1, argument[1]);
                muParameter.setParameterValue(0, argument[2]);
                double lnL = -treeLikelihood.getLogLikelihood();

                // System.err.println("" + argument[0] + "\t" + argument[1] + "\t" + argument[2] + "\t" + lnL);
                return lnL;
            }

            public int getNumArguments() {
                return 3;
            }

            public double getLowerBound(int n) {
                return 0.0;
            }

            public double getUpperBound(int n) {
                return 100.0;
            }
        };

        MultivariateMinimum optimizer = new ConjugateGradientSearch();
        double lnL = optimizer.findMinimum(function, new double[] { 1.0, 1.0, 1.0}, 6, 6);
        pw.write(betaParameter.getParameterValue(0)+"\t");
        pw.write(betaParameter.getParameterValue(1)+"\t");
        pw.write(muParameter.getParameterValue(0)+"\t");
        pw.write(lnL+"\n");
        pw.flush();
        System.out.println("" + betaParameter.getParameterValue(0) + "\t" + betaParameter.getParameterValue(1)+"\t" + muParameter.getParameterValue(0)+"\t" + lnL);
    }

    /**
     * @param fileName
     * @throws java.io.IOException
     */
    private Alignment readNexusFile(String fileName) throws java.io.IOException {

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
    private Tree readTreeFile(String fileName) throws java.io.IOException {

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
    static public void main(String[] args) throws java.io.IOException {
        if (args.length > 2) {
            System.err.println("Unknown option: " + args[3]);
            System.err.println();
            System.out.println("Usage: misc <alignmentfilename> <treefilename> <outputFileName>");
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

        new Coevolve(inputFileName, treeFileName);
    }

}

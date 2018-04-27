/*
 * CheckPointUpdater.java
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

package dr.app.realtime;

import dr.app.beast.BeastParser;
import dr.app.checkpoint.BeastCheckpointer;
import dr.app.util.Arguments;
import dr.evolution.alignment.PatternList;
import dr.evolution.distance.DistanceMatrix;
import dr.evolution.distance.F84DistanceMatrix;
import dr.evolution.distance.JukesCantorDistanceMatrix;
import dr.evolution.util.Taxon;
import dr.inference.loggers.Logger;
import dr.inference.loggers.MCLogger;
import dr.inference.markovchain.MarkovChain;
import dr.inference.mcmc.MCMC;
import dr.inference.model.Likelihood;
import dr.util.Transform;
import dr.xml.XMLParseException;
import dr.xml.XMLParser;
import org.xml.sax.SAXException;

import javax.xml.parsers.ParserConfigurationException;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Locale;
import java.util.Set;

/**
 * @author Guy Baele
 */
public class CheckPointUpdaterApp {

    private final boolean ADD_TAXA = true;

    private final boolean VERBOSE = true;
    private final boolean PARSER_WARNINGS = true;
    private final boolean STRICT_XML = false;

    public enum UpdateChoice {
        JC69DISTANCE("JC69Distance", new JukesCantorDistanceMatrix()), F84DISTANCE("F84Distance", new F84DistanceMatrix()), SIMPLE("Simple", new SimpleDistanceMatrix());

        private String name;
        private DistanceMatrix matrix;

        UpdateChoice(String name, DistanceMatrix matrix) {
            this.name = name;
            this.matrix = matrix;
        }

        public void setPatterns(PatternList patterns) {
            this.matrix.setPatterns(patterns);
        }

        public Taxon getClosestTaxon(Taxon taxon) {
            if (matrix == null) {
                throw new RuntimeException("Patterns need to be set first.");
            }
            int taxonIndex = matrix.getTaxonIndex(taxon);
            int closestIndex = 0;
            double minimumDistance = Double.MAX_VALUE;
            for (int i = 0; i < matrix.getColumnCount(); i++) {
                if (i != taxonIndex) {
                    if (matrix.getElement(taxonIndex, i) < minimumDistance) {
                        minimumDistance = matrix.getElement(taxonIndex, i);
                        closestIndex = i;
                    }
                }
            }
            return matrix.getTaxon(closestIndex);
        }

        public Taxon getClosestTaxon(Taxon taxon, ArrayList<Taxon> taxa) {
            if (matrix == null) {
                throw new RuntimeException("Patterns need to be set first.");
            }
            int taxonIndex = matrix.getTaxonIndex(taxon);
            int closestIndex = 0;
            double minimumDistance = Double.MAX_VALUE;
            for (int i = 0; i < matrix.getColumnCount(); i++) {
                if (i != taxonIndex) {
                    if ((matrix.getElement(taxonIndex, i) < minimumDistance) && (taxa.contains(matrix.getTaxon(i)))) {
                        minimumDistance = matrix.getElement(taxonIndex, i);
                        closestIndex = i;
                    }
                }
            }
            return matrix.getTaxon(closestIndex);
        }

        public double getDistance(Taxon taxonOne, Taxon taxonTwo) {
            System.out.println("taxon 1: " + taxonOne + " (" + matrix.getTaxonIndex(taxonOne) + ")");
            System.out.println("taxon 2: " + taxonTwo + " (" + matrix.getTaxonIndex(taxonTwo) + ")");
            return matrix.getElement(matrix.getTaxonIndex(taxonOne), matrix.getTaxonIndex(taxonTwo));
        }

        public String getName() {
            return this.name;
        }

        public DistanceMatrix getMatrix() {
            return this.matrix;
        }

        public String toString() {
            return this.name;
        }
    }

    /**
     * The goal is to modify an existing checkpoint file with additional information and generate a novel checkpoint file.
     * Running the MCMC chain after parsing the file(s) should not happen.
     * @param beastXMLFileName
     */
    public CheckPointUpdaterApp(String beastXMLFileName, String debugStateFile, UpdateChoice choice) {
        //no additional parsers, we don't need BEAGLE at the moment just yet
        XMLParser parser = new BeastParser(new String[]{beastXMLFileName}, null, VERBOSE, PARSER_WARNINGS, STRICT_XML, null);
        try {
            FileReader fileReader = new FileReader(beastXMLFileName);

            //Don't run the analysis, so set the argument to false
            //parser.parse(fileReader, false);

            //get the MCMC object
            MCMC mcmc = (MCMC) parser.parse(fileReader, MCMC.class);
            MarkovChain mc = mcmc.getMarkovChain();

            //make sure that no output files (i.e. logs) are being created
            Logger[] loggers = mcmc.getLoggers();
            for (int j = 0; j < loggers.length; j++) {
                ((MCLogger) loggers[j]).setFormatters(Collections.EMPTY_LIST);
            }

            // Install the checkpointer. This creates a factory that returns
            // appropriate savers and loaders according to the user's options.
            // BeastCheckpointer checkpoint = new BeastCheckpointer();

            CheckPointModifier checkpoint = new CheckPointModifier();

            //load the stored checkpoint file
            //this will/should also copy any trait information present in the checkpoint file
            //TODO Check if this works for multiple trees (e.g. for multiple partitions)
            long state = checkpoint.loadState(mc, new double[]{Double.NaN});

            //TODO Check if this can be uncommented again
            //System.out.println("Pre-checking likelihood values ...");
            //probably don't need this but it's good to check
            //double logL = mc.evaluate();
            //System.out.println("likelihood = " + logL);
            //mc.getLikelihood().makeDirty();
            //logL = mc.evaluate();
            //System.out.println("likelihood = " + logL);

            double logL;

            if (ADD_TAXA) {

                checkpoint.extendLoadState(choice);

                mc.getLikelihood().makeDirty();
                logL = mc.evaluate();
                System.out.println("likelihood = " + logL);
                mc.getLikelihood().makeDirty();
                logL = mc.evaluate();
                System.out.println("likelihood = " + logL);

                //TODO Print full compoundLikelihood evaluation
                Set<Likelihood> likelihoodSet = mc.getLikelihood().getLikelihoodSet();
                for (Likelihood l : likelihoodSet) {
                    System.out.println("  " + l.getLogLikelihood());
                }

            }

            checkpoint.saveState(mc, state, logL);

            //TODO .log and .trees files are being created; not necessary here as we're not running an analysis

            fileReader.close();

        } catch (FileNotFoundException fnf) {
            System.out.println(fnf);
        } catch (IOException io) {
            System.out.println(io);
            io.printStackTrace();
        } catch (SAXException sax) {
            System.out.println(sax);
        } catch (XMLParseException xml) {
            System.out.println(xml);
        } catch (ParserConfigurationException pce) {
            System.out.println(pce);
        }
    }

    public static void main(String[] args) throws java.io.IOException {

        // There is a major issue with languages that use the comma as a decimal separator.
        // To ensure compatibility between programs in the package, enforce the US locale.
        Locale.setDefault(Locale.US);

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.StringOption("BEAST_XML", "FILENAME", "Specify a BEAST XML file"),
                        new Arguments.StringOption("load_dump", "FILENAME", "Specify a filename to load a dumped state from"),
                        new Arguments.StringOption("output_file", "FILENAME", "Specify a filename for the output file"),
                        new Arguments.StringOption("update_choice", "UPDATECHOICE", "Specify a function by which to update the tree"),
                        new Arguments.Option("help", "Print this information and stop")
                });

        try {
            arguments.parseArguments(args);
        } catch (Arguments.ArgumentException ae) {
            System.out.println();
            System.out.println(ae.getMessage());
            System.out.println();
            //printUsage(arguments);
            System.exit(1);
        }

        String inputFile = null;

        if (arguments.hasOption("BEAST_XML")) {
            inputFile = arguments.getStringOption("BEAST_XML");
        } else {
            throw new RuntimeException("No BEAST XML file specified.");
        }

        String debugStateFile;
        if (arguments.hasOption("load_dump")) {
            debugStateFile = arguments.getStringOption("load_dump");
            //pass on as argument
            System.setProperty(BeastCheckpointer.LOAD_STATE_FILE, debugStateFile);
        } else {
            throw new RuntimeException("No dump file specified.");
        }

        String choice = "";
        if (arguments.hasOption("update_choice")) {
            choice = arguments.getStringOption("update_choice");
        } else {
            throw new RuntimeException("Update mechanism needs to be specified.");
        }
        UpdateChoice chosen = null;
        for (UpdateChoice ch : UpdateChoice.values()) {
            if (choice.equals(ch.getName())) {
                chosen = ch;
                break;
            }
        }
        if (chosen == null) {
            throw new RuntimeException("Incorrect update mechanism specified.");
        }

        if (arguments.hasOption("output_file")) {
            String outputStateFile = arguments.getStringOption("output_file");
            //pass on as argument
            System.setProperty(BeastCheckpointer.SAVE_STATE_FILE, outputStateFile);
        } else {
            throw new RuntimeException("No output file specified.");
        }

        new CheckPointUpdaterApp(inputFile, debugStateFile, chosen);

        System.exit(0);

    }

}



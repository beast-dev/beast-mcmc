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
import dr.evolution.alignment.Patterns;
import dr.evolution.distance.DistanceMatrix;
import dr.evolution.distance.F84DistanceMatrix;
import dr.evolution.distance.JukesCantorDistanceMatrix;
import dr.evolution.tree.BranchRates;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;
import dr.evomodel.tree.TreeModel;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.MultiPartitionDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.markovchain.MarkovChain;
import dr.inference.mcmc.MCMC;
import dr.inference.model.Likelihood;
import dr.inference.model.Model;
import dr.xml.XMLParseException;
import dr.xml.XMLParser;
import org.xml.sax.SAXException;

import javax.xml.parsers.ParserConfigurationException;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Locale;

/**
 * @author Guy Baele
 */
public class CheckPointUpdaterApp {

    private final boolean VERBOSE = true;
    private final boolean PARSER_WARNINGS = true;
    private final boolean STRICT_XML = false;

    public enum UpdateChoice {
        JC69DISTANCE("JC69Distance", new JukesCantorDistanceMatrix()), F84DISTANCE("F84Distance", new F84DistanceMatrix());

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

        public String getName() {
            return this.name;
        }

        public String toString() {
            return this.name;
        }
    }

    /**
     * The goal is to modify and existing checkpoint file with additional information and generate a novel checkpoint file.
     * Running the MCMC chain after parsing the file(s) should not happen.
     * @param beastXMLFileName
     */
    public CheckPointUpdaterApp(String beastXMLFileName, String debugStateFile, UpdateChoice choice, FileWriter fw) {
        //no additional parsers, we don't need BEAGLE at the moment just yet
        XMLParser parser = new BeastParser(new String[]{beastXMLFileName}, null, VERBOSE, PARSER_WARNINGS, STRICT_XML);
        try {
            FileReader fileReader = new FileReader(beastXMLFileName);

            //Don't run the analysis, so set the argument to false
            //parser.parse(fileReader, false);

            //get the MCMC object
            MCMC mcmc = (MCMC) parser.parse(fileReader, MCMC.class);
            MarkovChain mc = mcmc.getMarkovChain();

            // Install the checkpointer. This creates a factory that returns
            // appropriate savers and loaders according to the user's options.
            // BeastCheckpointer checkpoint = new BeastCheckpointer();

            CheckPointModifier checkpoint = new CheckPointModifier();

            //load the stored checkpoint file
            long state = checkpoint.loadState(mc, new double[]{Double.NaN});

            //probably don't need this but it's good to check
            double lnL = mc.evaluate();
            System.out.println("likelihood = " + lnL);
            mc.getLikelihood().makeDirty();
            lnL = mc.evaluate();
            System.out.println("likelihood = " + lnL);

            checkpoint.extendLoadState(choice);

            lnL = mc.evaluate();
            System.out.println("likelihood = " + lnL);
            mc.getLikelihood().makeDirty();
            lnL = mc.evaluate();
            System.out.println("likelihood = " + lnL);





            System.exit(0);






            //TODO Move most of the code below to the CheckPointModifier class and keep the app code simple?
            //check the Tree(Data)Likelihoods in the connected set of likelihoods
            //focus on TreeDataLikelihood, which has getTree() to get the tree for each likelihood
            //also get the DataLikelihoodDelegate from TreeDataLikelihood
            ArrayList<TreeDataLikelihood> likelihoods = new ArrayList<TreeDataLikelihood>();
            ArrayList<Tree> trees = new ArrayList<Tree>();
            ArrayList<DataLikelihoodDelegate> delegates = new ArrayList<DataLikelihoodDelegate>();
            for (Likelihood likelihood : Likelihood.CONNECTED_LIKELIHOOD_SET) {
                if (likelihood instanceof TreeDataLikelihood) {
                    likelihoods.add((TreeDataLikelihood)likelihood);
                    trees.add(((TreeDataLikelihood) likelihood).getTree());
                    delegates.add(((TreeDataLikelihood) likelihood).getDataLikelihoodDelegate());
                }
            }

            BranchRates rateModel = null;
            for (Model model : Model.CONNECTED_MODEL_SET) {
                if (model instanceof BranchRates) {
                    rateModel = (BranchRates) model;
                }
            }

            //suggested to go through TreeDataLikelihoodParser and give it an extra option to create a HashMap
            //keyed by the tree; am currently not overly fond of this approach
            ArrayList<PatternList> patternLists = new ArrayList<PatternList>();
            for (DataLikelihoodDelegate del : delegates) {
                if (del instanceof BeagleDataLikelihoodDelegate) {
                    patternLists.add(((BeagleDataLikelihoodDelegate) del).getPatternList());
                } else if (del instanceof MultiPartitionDataLikelihoodDelegate) {
                    //TODO complete code
                }
            }

            //aggregate all patterns to create distance matrix
            //TODO What about different trees for different partitions?
            Patterns patterns = new Patterns(patternLists.get(0));
            if (patternLists.size() > 1) {
                for (int i = 1; i < patternLists.size(); i++) {
                    patterns.addPatterns(patternLists.get(i));
                }
            }

            //identify the additional taxa/sequences
            int taxonCount = patterns.getTaxonCount();

            //TODO Could be multiple trees (for multiple partitions)
            Tree currentTree = trees.get(0);
            int treeTaxa = currentTree.getExternalNodeCount();

            if (taxonCount > treeTaxa) {
                System.out.println("Additional taxa found:");
                //list to contain the taxa being added
                ArrayList<Taxon> additionalTaxa = new ArrayList<Taxon>();

                Iterator<Taxon> iterator = patterns.iterator();
                while (iterator.hasNext()) {
                    Taxon taxon = iterator.next();
                    boolean taxonFound = false;
                    for (int i = 0; i < currentTree.getExternalNodeCount(); i++) {
                        if (currentTree.getNodeTaxon(currentTree.getExternalNode(i)) == taxon) {
                            taxonFound = true;
                        }
                    }
                    if (!taxonFound) {
                        additionalTaxa.add(taxon);
                    }
                }

                for (Taxon tax : additionalTaxa) {
                    System.out.println(tax);
                }

                TreeModel newTreeModel = null;
                if (choice.equals(UpdateChoice.JC69DISTANCE)) {
                    //build a distance matrix according to JC69
                    JukesCantorDistanceMatrix jcDistanceMatrix = new JukesCantorDistanceMatrix(patterns);
                    //newTreeModel = new GeneticDistanceTree(currentTree, rateModel).addTaxa(additionalTaxa,jcDistanceMatrix);
                } else if (choice.equals(UpdateChoice.F84DISTANCE)) {
                    //build a distance matrix according to F84
                    F84DistanceMatrix f84DistanceMatrix = new F84DistanceMatrix(patterns);
                    //newTreeModel = new GeneticDistanceTree(currentTree, rateModel).addTaxa(additionalTaxa,f84DistanceMatrix);
                } else {
                    throw new RuntimeException("Invalid update option provided.");
                }

            } else {
                throw new RuntimeException("Removing taxa from previous analysis currently not supported.");
            }

            mc.getLikelihood().makeDirty();
            double logL = mc.evaluate();

            System.out.println("new logLikelihood value = " + logL);

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

        FileWriter fw = null;
        if (arguments.hasOption("output_file")) {
            String output = arguments.getStringOption("output_file");
            fw = new FileWriter(output);
        } else {
            throw new RuntimeException("No output file specified.");
        }

        new CheckPointUpdaterApp(inputFile, debugStateFile, chosen, fw);

        fw.flush();
        fw.close();

        System.exit(0);

    }

}



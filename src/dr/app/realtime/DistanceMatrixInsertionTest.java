/*
 * TaxonInsertionTest.java
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
import dr.app.util.Arguments;
import dr.evolution.alignment.PatternList;
import dr.evolution.alignment.Patterns;
import dr.evolution.distance.DistanceMatrix;
import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.BeagleDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.MultiPartitionDataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.inference.markovchain.MarkovChain;
import dr.inference.mcmc.MCMC;
import dr.inference.model.Likelihood;
import dr.xml.XMLParseException;
import dr.xml.XMLParser;
import org.xml.sax.SAXException;

import javax.swing.*;
import javax.xml.parsers.ParserConfigurationException;
import java.awt.*;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

/**
 * @author Guy Baele
 */
public class DistanceMatrixInsertionTest extends JPanel {

    private final boolean VERBOSE = true;
    private final boolean PARSER_WARNINGS = true;
    private final boolean STRICT_XML = false;

    public DistanceMatrixInsertionTest(String beastXMLFileName, CheckPointUpdaterApp.UpdateChoice choice) {

        super(new GridLayout(1,0));

        //no additional parsers, we don't need BEAGLE at the moment just yet
        XMLParser parser = new BeastParser(new String[]{beastXMLFileName}, null, VERBOSE, PARSER_WARNINGS, STRICT_XML, null);
        try {
            FileReader fileReader = new FileReader(beastXMLFileName);

            //get the MCMC object
            MCMC mcmc = (MCMC) parser.parse(fileReader, MCMC.class);
            MarkovChain mc = mcmc.getMarkovChain();

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

            //suggested to go through TreeDataLikelihoodParser and give it an extra option to create a HashMap
            //keyed by the tree; am currently not overly fond of this approach
            ArrayList<PatternList> patternLists = new ArrayList<PatternList>();
            for (DataLikelihoodDelegate del : delegates) {
                if (del instanceof BeagleDataLikelihoodDelegate) {
                    patternLists.add(((BeagleDataLikelihoodDelegate) del).getPatternList());
                } else if (del instanceof MultiPartitionDataLikelihoodDelegate) {
                    MultiPartitionDataLikelihoodDelegate mpdld = (MultiPartitionDataLikelihoodDelegate)del;
                    List<PatternList> list = mpdld.getPatternLists();
                    for (PatternList pList : list) {
                        patternLists.add(pList);
                    }
                }
            }

            if (patternLists.size() == 0) {
                throw new RuntimeException("No patterns detected. Please make sure the XML file is BEAST 1.9 compatible.");
            }

            //aggregate all patterns to create distance matrix
            Patterns patterns = new Patterns(patternLists.get(0));
            if (patternLists.size() > 1) {
                for (int i = 1; i < patternLists.size(); i++) {
                    patterns.addPatterns(patternLists.get(i));
                }
            }

            String[] columnNames = new String[patterns.getTaxonCount()+1];
            columnNames[0] = "";
            for (int i = 0; i < patterns.getTaxonCount(); i++) {
                columnNames[i+1] = patterns.getTaxonId(i);
            }

            //set the patterns for the distance matrix
            DistanceMatrix matrix = choice.getMatrix();

            Object[][] data = new Object[patterns.getTaxonCount()][patterns.getTaxonCount()+1];

            //perform actual tests on the constructed distance matrix
            for (int i = 0; i < matrix.getRowCount(); i++) {
                data[i][0] = new String(columnNames[i+1]);
                for (int j = 1; j < matrix.getColumnCount()+1; j++) {
                    data[i][j] = new Double(matrix.getElement(i, j-1));
                }
            }

            final JTable table = new JTable(data, columnNames);
            //table.setPreferredScrollableViewportSize(new Dimension(800, 600));
            table.setFillsViewportHeight(true);
            table.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);

            for (int i = 0; i < patterns.getTaxonCount(); i++) {
                table.getColumnModel().getColumn(i).setPreferredWidth(100);
            }

            //Create the scroll pane and add the table to it.
            JScrollPane scrollPane = new JScrollPane(table, JScrollPane.VERTICAL_SCROLLBAR_ALWAYS, JScrollPane.HORIZONTAL_SCROLLBAR_ALWAYS);

            //Add the scroll pane to this panel.
            add(scrollPane);

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

    public static void main(String[] args) {

        // There is a major issue with languages that use the comma as a decimal separator.
        // To ensure compatibility between programs in the package, enforce the US locale.
        Locale.setDefault(Locale.US);

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.StringOption("BEAST_XML", "FILENAME", "Specify a BEAST XML file"),
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

        String choice = "";
        if (arguments.hasOption("update_choice")) {
            choice = arguments.getStringOption("update_choice");
        } else {
            throw new RuntimeException("Update mechanism needs to be specified.");
        }
        CheckPointUpdaterApp.UpdateChoice chosen = null;
        for (CheckPointUpdaterApp.UpdateChoice ch : CheckPointUpdaterApp.UpdateChoice.values()) {
            if (choice.equals(ch.getName())) {
                chosen = ch;
                break;
            }
        }
        if (chosen == null) {
            throw new RuntimeException("Incorrect update mechanism specified.");
        }

        //Create and set up the window.
        JFrame frame = new JFrame("Test distance matrix");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);

        DistanceMatrixInsertionTest test = new DistanceMatrixInsertionTest(inputFile, chosen);
        test.setOpaque(true);
        frame.setContentPane(test);

        //Display the window.
        frame.pack();
        frame.setVisible(true);

        //System.exit(0);

    }

}

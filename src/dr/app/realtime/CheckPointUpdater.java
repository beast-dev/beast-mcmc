/*
 * TopologyTracer.java
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
import dr.inference.markovchain.MarkovChain;
import dr.inference.mcmc.MCMC;
import dr.xml.XMLParseException;
import dr.xml.XMLParser;
import org.xml.sax.SAXException;

import javax.xml.parsers.ParserConfigurationException;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Locale;

/**
 * @author Guy Baele
 */
public class CheckPointUpdater {

    private final boolean VERBOSE = true;
    private final boolean PARSER_WARNINGS = true;
    private final boolean STRICT_XML = false;

    /**
     * The goal is to modify and existing checkpoint file with additional information and generate a novel checkpoint file.
     * Running the MCMC chain after parsing the file(s) should not happen.
     * @param beastXMLFileName
     */
    public CheckPointUpdater(String beastXMLFileName, String debugStateFile, FileWriter fw) {
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
            BeastCheckpointer checkpoint = new BeastCheckpointer();

            //load the stored checkpoint file
            checkpoint.loadState(mc, new double[]{Double.NaN});




            /*Iterator iter = parser.getParsers();
            while (iter.hasNext()) {
                Object nextParser = iter.next();
                if (nextParser instanceof AlignmentParser) {
                    System.out.println("Alignment parser found!");
                    AlignmentParser alParser = (AlignmentParser) nextParser;

                    Map<String, XMLObject> map = alParser.getStore();
                    System.out.println(map.size());
                }
            }*/





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

        FileWriter fw = null;
        if (arguments.hasOption("output_file")) {
            String output = arguments.getStringOption("output_file");
            fw = new FileWriter(output);
        } else {
            throw new RuntimeException("No output file specified.");
        }

        new CheckPointUpdater(inputFile, debugStateFile, fw);

        fw.flush();
        fw.close();

        System.exit(0);

    }

}



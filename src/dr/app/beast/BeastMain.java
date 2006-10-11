/*
 * BeastMain.java
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

package dr.app.beast;

import dr.app.util.Arguments;
import dr.app.util.Utils;
import dr.util.Version;
import dr.util.MessageLogHandler;
import dr.xml.XMLParser;
import dr.math.MathUtils;

import java.io.*;
import java.util.Iterator;
import java.util.logging.*;

import org.virion.jam.util.IconUtils;

public class BeastMain {

    private final static Version version = new BeastVersion();

    static class BeastConsoleApp extends org.virion.jam.console.ConsoleApplication {
        XMLParser parser = null;

        public BeastConsoleApp(String nameString, String aboutString, javax.swing.Icon icon) throws IOException {
            super(nameString, aboutString, icon);
        }

        public void doStop() {
            Iterator iter = parser.getThreads();
            while (iter.hasNext()) {
                Thread thread = (Thread)iter.next();
                thread.stop();
            }
        }
    };

    public BeastMain(File inputFile, BeastConsoleApp consoleApp, boolean verbose) {

        if (inputFile == null) {
            System.err.println();
            System.err.println("Error: no input file specified");
            return;
        }

        try {
            String fileName = inputFile.getName();

            FileReader fileReader = new FileReader(inputFile);
            System.out.println(fileReader.getEncoding());

            XMLParser parser = new BeastParser(new String[] {fileName}, verbose);

            if (consoleApp != null) {
                consoleApp.parser = parser;
            }

            // Add a handler to handle warnings and errors. This is a ConsoleHandler
            // so the messages will go to StdOut..
            Logger logger = Logger.getLogger("dr");
            Handler handler = new MessageLogHandler();
            handler.setFilter(new Filter() {
                public boolean isLoggable(LogRecord record) {
                    return record.getLevel().intValue() < Level.WARNING.intValue();
                }
            });
            logger.addHandler(handler);

            // Add a handler to handle warnings and errors. This is a ConsoleHandler
            // so the messages will go to StdErr..
            handler = new ConsoleHandler();
            handler.setFilter(new Filter() {
                public boolean isLoggable(LogRecord record) {
                    return record.getLevel().intValue() >= Level.WARNING.intValue();
                }
            });
            logger.addHandler(handler);

            logger.setUseParentHandlers(false);

            Logger.getLogger("dr.apps.beast").info("Parsing XML file: " + fileName);

            parser.parse(fileReader, true);

        } catch (java.io.IOException ioe) {
            Logger.getLogger("dr.apps.beast").severe("File error: " + ioe.getMessage());
        } catch (org.xml.sax.SAXParseException spe) {
            if (spe.getMessage() != null && spe.getMessage().equals("Content is not allowed in prolog")) {
                Logger.getLogger("dr.apps.beast").severe("Parsing error - the input file is not a valid XML file.");
            } else {
                Logger.getLogger("dr.apps.beast").severe("Parsing error - poorly formed XML (possibly not an XML file):\n" +
                        spe.getMessage());
            }
        } catch (org.w3c.dom.DOMException dome) {
            Logger.getLogger("dr.apps.beast").severe("Parsing error - poorly formed XML:\n" +
                    dome.getMessage());
        } catch (dr.xml.XMLParseException pxe) {
            if (pxe.getMessage() != null && pxe.getMessage().equals("Unknown root document element, beauti")) {
                Logger.getLogger("dr.apps.beast").severe("The file you just tried to run in BEAST is actually a BEAUti document.\n" +
                        "Although this uses XML, it is not a format that BEAST understands.\n" +
                        "These files are used by BEAUti to save and load your settings so that\n" +
                        "you can go back and alter them. To generate a BEAST file you must\n" +
                        "select the 'Generate BEAST File' option, either from the File menu or\n" +
                        "the button at the bottom right of the window.");

            } else {
                Logger.getLogger("dr.apps.beast").severe("Parsing error - poorly formed BEAST file:\n" +
                        pxe.getMessage());
            }

        } catch (RuntimeException rex) {
            if (rex.getMessage() != null && rex.getMessage().startsWith("The initial model is invalid")) {
                Logger.getLogger("dr.apps.beast").severe("The initial model is invalid because state has a zero likelihood.\n" +
                        "This may be because the initial, random tree is so large that it\n" +
                        "has an extremely bad likelihood which is being rounded to zero.\n" +
                        "Alternatively, it may be that the product of starting mutation rate\n" +
                        "and tree height is extremely small or extremely large. Try to set\n" +
                        "initial values such that the product is similar to the average\n" +
                        "pairwise genetic distance between the sequences.\n" +
                        "For more information go to <http://evolve.zoo.ox.ac.uk/beast/help/>.");

            } else {
                Logger.getLogger("dr.apps.beast").severe("Fatal exception (email the authors)");
                System.err.println("Fatal exception (email the authors)");
                rex.printStackTrace(System.err);
            }

        } catch (Exception ex) {
            Logger.getLogger("dr.apps.beast").severe("Fatal exception (email the authors)");
            System.err.println("Fatal exception (email the authors)");
            ex.printStackTrace(System.err);
        }
    }

    public static void printTitle() {

        System.out.println("+-----------------------------------------------\\");

        String versionString = "BEAST " + version.getVersionString() + " 2002-2006";
        System.out.print("|");
        int n = 47 - versionString.length();
        int n1 = n / 2;
        int n2 = n1 + (n % 2);
        for (int i = 0; i < n1; i++) { System.out.print(" "); }
        System.out.print(versionString);
        for (int i = 0; i < n2; i++) { System.out.print(" "); }
        System.out.println("|\\");

        System.out.println("| Bayesian Evolutionary Analysis Sampling Trees ||");

        String buildString = "BEAST Library: " + version.getBuildString();
        System.out.print("|");
        n = 47 - buildString.length();
        n1 = n / 2;
        n2 = n1 + (n % 2);
        for (int i = 0; i < n1; i++) { System.out.print(" "); }
        System.out.print(buildString);
        for (int i = 0; i < n2; i++) { System.out.print(" "); }
        System.out.println("||");

        System.out.println("|       Alexei Drummond and Andrew Rambaut      ||");
        System.out.println("|              University of Oxford             ||");
        System.out.println("|       http://evolve.zoo.ox.ac.uk/Beast/       ||");
        System.out.println("\\-----------------------------------------------\\|");
        System.out.println(" \\-----------------------------------------------\\");
    }

    public static void printHeader() {
        System.out.println(" +-----------------------------------------------+");
        System.out.println(" | Components created by:                        |");
        System.out.println(" |       Alexei Drummond                         |");
        System.out.println(" |       Roald Forsberg                          |");
        System.out.println(" |       Gerton Lunter                           |");
        System.out.println(" |       Oliver Pybus                            |");
        System.out.println(" |       Andrew Rambaut                          |");
        System.out.println(" | Thanks to (for use of their code):            |");
        System.out.println(" |       Korbinian Strimmer                      |");
        System.out.println(" |       Oliver Pybus                            |");
        System.out.println(" +-----------------------------------------------+");
        System.out.println();

    }

    public static void printUsage(Arguments arguments) {

        arguments.printUsage("beast", "[<input-file-name>]");
        System.out.println();
        System.out.println("  Example: beast test.xml");
        System.out.println("  Example: beast -window test.xml");
        System.out.println();
    }

    //Main method
    public static void main(String[] args) throws java.io.IOException {


        Arguments arguments = new Arguments(
            new Arguments.Option[] {
                new Arguments.Option("verbose", "verbose XML parsing messages"),
                new Arguments.Option("window", "provide a console window"),
                new Arguments.Option("working", "change working directory to input file's directory"),
                        new Arguments.LongOption("seed", "specify a random number generator seed"),
                new Arguments.Option("help", "option to print this message")
            });

        try {
            arguments.parseArguments(args);
        } catch (Arguments.ArgumentException ae) {
            System.out.println();
            System.out.println(ae.getMessage());
            System.out.println();
            printTitle();
            printUsage(arguments);
            System.exit(1);
        }

        if (arguments.hasOption("help")) {
            printTitle();
            printUsage(arguments);
            System.exit(0);
        }

        boolean verbose = arguments.hasOption("verbose");
        boolean window = arguments.hasOption("window");
        boolean working = arguments.hasOption("working");

        long seed = MathUtils.getSeed();
        if (arguments.hasOption("seed")) {
            seed = arguments.getLongOption("seed");
            MathUtils.setSeed(seed);
        }

//		if (System.getProperty("dr.app.beast.main.window", "false").toLowerCase().equals("true")) {
//			window = true;
//		}

        BeastConsoleApp consoleApp = null;

        if (window) {
            System.setProperty("com.apple.macos.useScreenMenuBar","true");
            System.setProperty("apple.laf.useScreenMenuBar","true");
            System.setProperty("apple.awt.showGrowBox","true");

            javax.swing.Icon icon = IconUtils.getIcon(BeastMain.class, "images/beast.png");

            String nameString = "BEAST " + version.getVersionString();
            String aboutString = "Bayesian Evolutionary Analysis Sampling Trees\n" +
                                    version.getVersionString() + " ©2002-2006 Alexei Drummond & Andrew Rambaut\n" +
                                    "University of Oxford";

            consoleApp = new BeastConsoleApp(nameString, aboutString, icon);
        }

        String inputFileName = null;

        String[] args2 = arguments.getLeftoverArguments();

        if (args2.length > 1) {
            System.err.println("Unknown option: " + args2[1]);
            System.err.println();
            printTitle();
            printUsage(arguments);
            System.exit(1);
        }

        File inputFile = null;

        if (args2.length > 0) {
            inputFileName = args2[0];
            inputFile = new File(inputFileName);
        }

        if (inputFileName == null) {
            // No input file name was given so throw up a dialog box...
            inputFile = Utils.getLoadFile("BEAST "+version.getVersionString()+" - Select XML input file");
        }

        if (inputFile != null && working) {
            System.setProperty("user.dir", inputFile.getParent());
        }

        printTitle();
        printHeader();

        System.out.println();
        System.out.println("Random number seed: " + seed);
        System.out.println();

        new BeastMain(inputFile, consoleApp, verbose);
    }
}


/*
 * BeastMC3.java
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
import dr.inference.mcmc.MCMC;
import dr.inference.mcmcmc.MCMCMC;
import dr.inference.mcmcmc.MCMCMCOptions;
import dr.util.MessageLogHandler;
import dr.util.Version;
import dr.xml.XMLParser;
import jam.console.ConsoleApplication;
import jam.util.IconUtils;

import javax.swing.*;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Iterator;
import java.util.Locale;
import java.util.logging.*;

public class BeastMC3 {

    public final static int HOT_CHAIN_COUNT = 2;
    public static final double DELTA = 1.0;
    public static final int SWAP_CHAIN_EVERY = 100;

    private final static Version version = new dr.app.beast.BeastVersion();

    static class BeastConsoleApp extends ConsoleApplication {
        XMLParser parser = null;

        public BeastConsoleApp(String nameString, String aboutString, javax.swing.Icon icon) throws IOException {
            super(nameString, aboutString, icon, false);
        }

        public void doStop() {
            Iterator iter = parser.getThreads();
            while (iter.hasNext()) {
                Thread thread = (Thread) iter.next();
                // @todo should never use this method...
                thread.stop();
                // @todo who cares? - it is bad in complex multithreading applications but this is a really simple case.
                // At the moment it works and when it doesn't we can implement a polling alternative.
            }
        }
    }

    public BeastMC3(double[] chainTemperatures, int swapChainsEvery, File inputFile, BeastConsoleApp consoleApp,
                    boolean verbose, boolean parserWarning, boolean strictXML) {

        if (inputFile == null) {
            System.err.println();
            System.err.println("Error: no input file specified");
            return;
        }

        int chainCount = chainTemperatures.length;
        MCMC[] chains = new MCMC[chainCount];
        MCMCMCOptions options = new MCMCMCOptions();
        options.setChainTemperatures(chainTemperatures);
        options.setSwapChainsEvery(swapChainsEvery);

        // Add a handler to handle warnings and errors. This is a ConsoleHandler
        // so the messages will go to StdOut..
        Logger logger = Logger.getLogger("dr");
        Handler messageHandler = new MessageLogHandler();
        messageHandler.setFilter(new Filter() {
            public boolean isLoggable(LogRecord record) {
                return record.getLevel().intValue() < Level.WARNING.intValue();
            }
        });
        logger.addHandler(messageHandler);

        // Add a handler to handle warnings and errors. This is a ConsoleHandler
        // so the messages will go to StdErr..
        Handler errorHandler = new ConsoleHandler();
        errorHandler.setLevel(Level.WARNING);
        logger.addHandler(errorHandler);

        logger.setUseParentHandlers(false);

        try {
            String fileName = inputFile.getName();

            FileReader fileReader = new FileReader(inputFile);

            XMLParser parser = new BeastParser(new String[]{fileName}, null, verbose, parserWarning, strictXML);

            if (consoleApp != null) {
                consoleApp.parser = parser;
            }

            Logger.getLogger("dr.apps.beast").info("Starting cold chain plus hot chains with temperatures: ");
            for (int i = 1; i < chainTemperatures.length; i++) {
                Logger.getLogger("dr.apps.beast").info("Hot Chain " + i + ": " + chainTemperatures[i]);
            }

            Logger.getLogger("dr.apps.beast").info("Parsing XML file: " + fileName);

            chains[0] = (MCMC) parser.parse(fileReader, MCMC.class);
            if (chains[0] == null) {
                throw new dr.xml.XMLParseException("BEAST XML file is missing an MCMC element");
            }
            fileReader.close();

            chainTemperatures[0] = 1.0;

            for (int i = 1; i < chainCount; i++) {
                fileReader = new FileReader(inputFile);

                messageHandler.setLevel(Level.OFF);
                parser = new BeastParser(new String[]{fileName}, null, verbose, parserWarning, strictXML);

                chains[i] = (MCMC) parser.parse(fileReader, MCMC.class);
                if (chains[i] == null) {
                    throw new dr.xml.XMLParseException("BEAST XML file is missing an MCMC element");
                }
                fileReader.close();
            }
            messageHandler.setLevel(Level.ALL);

        } catch (IOException ioe) {
            System.err.println();
            System.err.println("File error:");
            System.err.println(ioe.getMessage());
        } catch (org.xml.sax.SAXParseException spe) {
            System.err.println();
            if (spe.getMessage() != null && spe.getMessage().equals("Content is not allowed in prolog")) {
                System.err.println("Parsing error - the input file is not a valid XML file.");
            } else {
                System.err.println("Parsing error - poorly formed XML (possibly not an XML file):");
                System.err.println(spe.getMessage());
            }
        } catch (org.w3c.dom.DOMException dome) {
            System.err.println();
            System.err.println("Parsing error - poorly formed XML:");
            System.err.println(dome.getMessage());
        } catch (dr.xml.XMLParseException pxe) {
            if (pxe.getMessage() != null && pxe.getMessage().equals("Unknown root document element, beauti")) {
                System.err.println();
                System.err.println("The file you just tried to run in BEAST is actually a BEAUti document.");
                System.err.println("Although this uses XML, it is not a format that BEAST understands. ");
                System.err.println("These files are used by BEAUti to save and load your settings so that ");
                System.err.println("you can go back and alter them. To generate a BEAST file you must ");
                System.err.println("select the 'Generate BEAST File' option, either from the File menu or ");
                System.err.println("the button at the bottom right of the window.");

            } else {

                System.err.println();
                System.err.println("Parsing error - poorly formed BEAST file:");
                System.err.println(pxe.getMessage());
            }
        } catch (RuntimeException rex) {
            if (rex.getMessage() != null && rex.getMessage().startsWith("The initial model is invalid")) {
                System.err.println();
                System.err.println("The initial model is invalid because state has a zero likelihood.");
                System.err.println("This may be because the initial, random tree is so large that it");
                System.err.println("has an extremely bad likelihood which is being rounded to zero.");
                System.err.println("Alternatively, it may be that the product of starting mutation rate");
                System.err.println("and tree height is extremely small or extremely large. Try to set");
                System.err.println("initial values such that the product is similar to the average");
                System.err.println("pairwise genetic distance between the sequences.");
                System.err.println("For more information go to <http://evolve.zoo.ox.ac.uk/beast/help/>.");

            } else {

                System.err.println();
                System.err.println("Fatal exception (email the authors)");
                rex.printStackTrace(System.out);
            }
        } catch (Exception ex) {
            System.err.println();
            System.err.println("Fatal exception (email the authors)");
            ex.printStackTrace(System.err);
        }

        MCMCMC mc3 = new MCMCMC(chains, options);
        Thread thread = new Thread(mc3);
        thread.start();

    }

    public static void printTitle() {

        System.out.println();
        BeastMain.centreLine("BEAST " + version.getVersionString() + ", " + version.getDateString(), 60);
        BeastMain.centreLine("Bayesian Evolutionary Analysis Sampling Trees", 60);
        for (String creditLine : version.getCredits()) {
            BeastMain.centreLine(creditLine, 60);
        }
        System.out.println();
    }

    public static void printUsage(Arguments arguments) {

        arguments.printUsage("beastmc3", " [-chains <chain-count>] [-delta <temperature-delta>|-temperatures <temperature-list>] [-swap <swap-every>] [-verbose] [-window] [-working] [-help] [<input-file-name>]");
        System.out.println();
        System.out.println("  Example: beastmc3 -chains 3 -delta 1.0 -swap 100 test.xml");
        System.out.println("  Example: beastmc3 -temperatures 0.2,0.4 -swap 200 -window test.xml");
        System.out.println();
    }

    //Main method
    public static void main(String[] args) throws IOException {

        // There is a major issue with languages that use the comma as a decimal separator.
        // To ensure compatibility between programs in the package, enforce the US locale.
        Locale.setDefault(Locale.US);

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.IntegerOption("chains", 2, Integer.MAX_VALUE, "number of chains"),
                        new Arguments.RealOption("delta", 0.0, Double.MAX_VALUE, "temperature increment parameter"),
                        new Arguments.RealArrayOption("temperatures", -1, "a comma-separated list of the hot chain temperatures"),
                        new Arguments.IntegerOption("swap", 1, Integer.MAX_VALUE, "frequency at which chains temperatures will be swapped"),
                        new Arguments.Option("verbose", "verbose XML parsing messages"),
                        new Arguments.Option("strict", "Fail on non conforming BEAST XML file"),
                        new Arguments.Option("window", "provide a console window"),
                        new Arguments.Option("working", "change working directory to input file's directory"),
                        new Arguments.Option("overwrite", "Allow overwriting of log files"),
                        new Arguments.Option("options", "Display an options dialog"),
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

        int chainCount = HOT_CHAIN_COUNT + 1;
        if (arguments.hasOption("chains")) {
            chainCount = arguments.getIntegerOption("chains");
        } else if (arguments.hasOption("temperatures")) {
            chainCount = 1 + arguments.getRealArrayOption("temperatures").length;
        }

        double delta = DELTA;
        if (arguments.hasOption("delta")) {
            if (arguments.hasOption("temperatures")) {
                System.err.println("Either the -delta or the -temperatures option should be used, not both");
                System.err.println();
                printUsage(arguments);
                System.exit(1);
            }
            delta = arguments.getRealOption("delta");
        }

        double[] chainTemperatures = new double[chainCount];
        chainTemperatures[0] = 1.0;
        if (arguments.hasOption("temperatures")) {
            double[] hotChainTemperatures = arguments.getRealArrayOption("temperatures");
            assert hotChainTemperatures.length == chainCount - 1;

            System.arraycopy(hotChainTemperatures, 0, chainTemperatures, 1, chainCount - 1);
        } else {
            for (int i = 1; i < chainCount; i++) {
                chainTemperatures[i] = 1.0 / (1.0 + (delta * i));
            }
        }

        int swapChainsEvery = SWAP_CHAIN_EVERY;
        if (arguments.hasOption("swap")) {
            swapChainsEvery = arguments.getIntegerOption("swap");
        }

        final boolean strictXML = arguments.hasOption("strict");
        final boolean verbose = arguments.hasOption("verbose");
        final boolean parserWarning = arguments.hasOption("pwarning"); // if dev, then auto turn on, otherwise default to turn off
        final boolean window = arguments.hasOption("window");
        final boolean working = arguments.hasOption("working");
        final boolean options = arguments.hasOption("options");
        boolean allowOverwrite = arguments.hasOption("overwrite");

//		if (System.getProperty("dr.app.beast.main.window", "false").toLowerCase().equals("true")) {
//			window = true;
//		}

        BeastConsoleApp consoleApp = null;

        String nameString = "BEAST " + version.getVersionString();

        if (window) {
            System.setProperty("com.apple.macos.useScreenMenuBar", "true");
            System.setProperty("apple.laf.useScreenMenuBar", "true");
            System.setProperty("apple.awt.showGrowBox", "true");

            java.net.URL url = BeastMC3.class.getResource("/images/beast.png");
            javax.swing.Icon icon = null;

            if (url != null) {
                icon = new javax.swing.ImageIcon(url);
            }

            String aboutString = "<html><div style=\"font-family:sans-serif;\"><center>" +
                    "<div style=\"font-size:12;\"><p>Bayesian Evolutionary Analysis Sampling Trees<br>" +
                    "Version " + version.getVersionString() + ", " + version.getDateString() + "</p>" +
                    version.getHTMLCredits() +
                    "</div></center></div></html>";

            consoleApp = new BeastConsoleApp(nameString, aboutString, icon);
        }

        File inputFile = null;

        if (options) {

            String titleString = "<html><center><p>Bayesian Evolutionary Analysis Sampling Trees<br>" +
                    "Metropolis-coupled version<br>" +
                    "Version " + version.getVersionString() + ", " + version.getDateString() + "</p></center></html>";
            javax.swing.Icon icon = IconUtils.getIcon(BeastMC3.class, "images/beast.png");

            BeastMC3Dialog dialog = new BeastMC3Dialog(new JFrame(), titleString, icon);

            if (!dialog.showDialog(nameString)) {
                return;
            }

            if (dialog.allowOverwrite()) {
                allowOverwrite = true;
            }

            chainCount = dialog.getChains();

            delta = dialog.getDelta();

            chainTemperatures = dialog.getTemperaturesArray();

            swapChainsEvery = dialog.getSwap();

            inputFile = dialog.getInputFile();

            if (inputFile == null) {
                System.err.println("No input file specified");
                return;
            }

        }


        if (allowOverwrite) {
            System.setProperty("log.allow.overwrite", "true");
        }

        if (inputFile == null) {
            String[] args2 = arguments.getLeftoverArguments();

            if (args2.length > 1) {
                System.err.println("Unknown option: " + args2[1]);
                System.err.println();
                printUsage(arguments);
                System.exit(1);
            }

            String inputFileName = null;

            if (args2.length > 0) {
                inputFileName = args2[0];
                inputFile = new File(inputFileName);
            }
            if (inputFileName == null) {
                // No input file name was given so throw up a dialog box...
                inputFile = Utils.getLoadFile("BEAST " + version.getVersionString() + " - Select XML input file");
            }
        }


        if (inputFile != null && working) {
            System.setProperty("user.dir", inputFile.getParent());
        }

        printTitle();

        try {
            new BeastMC3(chainTemperatures, swapChainsEvery, inputFile, consoleApp, verbose, parserWarning, strictXML);
        } catch (RuntimeException rte) {
            if (window) {
                System.out.println();
                System.out.println("BEAST has terminated with an error. Please select QUIT from the menu.");
                // logger.severe will throw a RTE but we want to keep the console visible
            } else {
                System.exit(1);
            }
        }

        if (!window) {
            System.exit(0);
        }
    }
}


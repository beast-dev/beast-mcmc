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
import dr.inference.mcmc.MCMC;
import dr.math.MathUtils;
import dr.util.ErrorLogHandler;
import dr.util.MessageLogHandler;
import dr.util.Version;
import dr.xml.XMLParser;
import org.virion.jam.util.IconUtils;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Iterator;
import java.util.logging.*;

public class BeastMain {

    private final static Version version = new BeastVersion();

    static class BeastConsoleApp extends org.virion.jam.console.ConsoleApplication {
        XMLParser parser = null;

        public BeastConsoleApp(String nameString, String aboutString, javax.swing.Icon icon) throws IOException {
            super(nameString, aboutString, icon, false);
        }

        public void doStop() {
            Iterator iter = parser.getThreads();
            while (iter.hasNext()) {
                Thread thread = (Thread) iter.next();
                thread.stop();
            }
        }
    }

    public BeastMain(File inputFile, BeastConsoleApp consoleApp, int maxErrorCount, boolean verbose, boolean strictXML) {

        if (inputFile == null) {
            System.err.println();
            System.err.println("Error: no input file specified");
            return;
        }

        String fileName = inputFile.getName();

        final Logger infoLogger = Logger.getLogger("dr.apps.beast");
        try {

            FileReader fileReader = new FileReader(inputFile);

            XMLParser parser = new BeastParser("release_parsers.properties", new String[]{fileName}, verbose, strictXML);

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
//            logger.addHandler(handler);

            logger.setUseParentHandlers(false);

            infoLogger.info("Parsing XML file: " + fileName);
            infoLogger.info("  File encoding: " + fileReader.getEncoding());

            // This is a special logger that is for logging numerical and statistical errors
            // during the MCMC run. It will tolerate up to maxErrorCount before throwing a
            // RuntimeException to shut down the run.
            //Logger errorLogger = Logger.getLogger("error");
            handler = new ErrorLogHandler(maxErrorCount);
            handler.setLevel(Level.WARNING);
            logger.addHandler(handler);

            parser.parse(fileReader, true);

        } catch (java.io.IOException ioe) {
            infoLogger.severe("File error: " + ioe.getMessage());
        } catch (org.xml.sax.SAXParseException spe) {
            if (spe.getMessage() != null && spe.getMessage().equals("Content is not allowed in prolog")) {
                infoLogger.severe("Parsing error - the input file, " + fileName + ", is not a valid XML file.");
            } else {
                infoLogger.severe("Error running file: " + fileName);
                infoLogger.severe("Parsing error - poorly formed XML (possibly not an XML file):\n" +
                        spe.getMessage());
            }
        } catch (org.w3c.dom.DOMException dome) {
            infoLogger.severe("Error running file: " + fileName);
            infoLogger.severe("Parsing error - poorly formed XML:\n" +
                    dome.getMessage());
        } catch (dr.xml.XMLParseException pxe) {
            if (pxe.getMessage() != null && pxe.getMessage().equals("Unknown root document element, beauti")) {
                infoLogger.severe("Error running file: " + fileName);
                infoLogger.severe(
                        "The file you just tried to run in BEAST is actually a BEAUti document.\n" +
                                "Although this uses XML, it is not a format that BEAST understands.\n" +
                                "These files are used by BEAUti to save and load your settings so that\n" +
                                "you can go back and alter them. To generate a BEAST file you must\n" +
                                "select the 'Generate BEAST File' option, either from the File menu or\n" +
                                "the button at the bottom right of the window.");

            } else {
                infoLogger.severe("Parsing error - poorly formed BEAST file, " + fileName + ":\n" +
                        pxe.getMessage());
            }

        } catch (RuntimeException rex) {
            if (rex.getMessage() != null && rex.getMessage().startsWith("The initial posterior is zero")) {
                infoLogger.warning("Error running file: " + fileName);
                infoLogger.severe(
                        "The initial model is invalid because state has a zero probability.\n\n" +
                                "If the log likelihood of the tree is -Inf, his may be because the\n" +
                                "initial, random tree is so large that it has an extremely bad\n" +
                                "likelihood which is being rounded to zero.\n\n" +
                                "Alternatively, it may be that the product of starting mutation rate\n" +
                                "and tree height is extremely small or extremely large. \n\n" +
                                "Finally, it may be that the initial state is incompatible with\n" +
                                "one or more 'hard' constraints (on monophyly or bounds on parameter\n" +
                                "values. This will result in Priors with zero probability.\n\n" +
                                "The individual components of the posterior are as follows:\n" +
                                rex.getMessage() + "\n" +
                                "For more information go to <http://beast.bio.ed.ac.uk/>.");

            } else {
                // This call never returns as another RuntimeException exception is raised by
                // the error log handler???
                infoLogger.warning("Error running file: " + fileName);
                System.err.println("Fatal exception: " + rex.getMessage());
                rex.printStackTrace(System.err);
            }

        } catch (Exception ex) {
            infoLogger.warning("Error running file: " + fileName);
            infoLogger.severe("Fatal exception: " + ex.getMessage());
            System.err.println("Fatal exception: " + ex.getMessage());
            ex.printStackTrace(System.err);
        }
    }

    public static void centreLine(String line, int pageWidth) {
        int n = pageWidth - line.length();
        int n1 = n / 2;
        for (int i = 0; i < n1; i++) {
            System.out.print(" ");
        }
        System.out.println(line);
    }

    public static void printTitle() {
        System.out.println();
        centreLine("BEAST " + version.getVersionString() + ", 2002-2008", 60);
        centreLine("Bayesian Evolutionary Analysis Sampling Trees", 60);
        centreLine("by", 60);
        centreLine("Alexei J. Drummond and Andrew Rambaut", 60);
        System.out.println();
        centreLine("Department of Computer Science", 60);
        centreLine("University of Auckland", 60);
        centreLine("alexei@cs.auckland.ac.nz", 60);
        System.out.println();
        centreLine("Institute of Evolutionary Biology", 60);
        centreLine("University of Edinburgh", 60);
        centreLine("a.rambaut@ed.ac.uk", 60);
        System.out.println();

    }

    public static void printHeader() {
        System.out.println("Downloads, Help & Resources:\n" +
                "\thttp://beast.bio.ed.ac.uk/\n" +
                "\n" +
                "Source code distributed under the GNU Lesser General Public License:\n" +
                "\thttp://code.google.com/p/beast-mcmc/\n" +
                "\n" +
                "Additional programming & components created by:\n" +
                "\tRoald Forsberg\n" +
                "\tGerton Lunter\n" +
                "\tSidney Markowitz\n" +
                "\tOliver Pybus\n" +
                "\n" +
                "Thanks to (for use of their code):\n" +
                "\tKorbinian Strimmer");
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

        printTitle();
        Arguments arguments = new Arguments(
                new Arguments.Option[]{

                        new Arguments.Option("verbose", "verbose XML parsing messages"),
                        new Arguments.Option("strict", "Fail on non conforming BEAST XML file"),
                        new Arguments.Option("window", "provide a console window"),
                        new Arguments.Option("working", "change working directory to input file's directory"),
                        new Arguments.LongOption("seed", "specify a random number generator seed"),
                        new Arguments.IntegerOption("errors", "maximum number of numerical errors before stopping"),
                        // new Arguments.Option("logops", "hack: log ops to stderr"),
                        new Arguments.IntegerOption("otfops", "experimental: on the fly op weigths. recompute frequency" +
                                "in number of states."),
                        new Arguments.Option("help", "option to print this message"),
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

        final boolean verbose = arguments.hasOption("verbose");
        final boolean strictXML = arguments.hasOption("strict");
        final boolean window = arguments.hasOption("window");
        final boolean working = arguments.hasOption("working");

        // (HACK)
        //MCMC.logOps =  arguments.hasOption("logops");
        MCMC.ontheflyFreq = arguments.hasOption("otfops") ? arguments.getIntegerOption("otfops") : 0;

        long seed = MathUtils.getSeed();
        if (arguments.hasOption("seed")) {
            seed = arguments.getLongOption("seed");
            if (seed <= 0) {
                printTitle();
                System.err.println("The random number seed should be > 0");
                System.exit(1);
            }
            MathUtils.setSeed(seed);
        }

        int maxErrorCount = 0;
        if (arguments.hasOption("errors")) {
            maxErrorCount = arguments.getIntegerOption("errors");
            if (maxErrorCount < 0) {
                maxErrorCount = 0;
            }
        }

//		if (System.getProperty("dr.app.beast.main.window", "false").toLowerCase().equals("true")) {
//			window = true;
//		}

        BeastConsoleApp consoleApp = null;

        if (window) {
            System.setProperty("com.apple.macos.useScreenMenuBar", "true");
            System.setProperty("apple.laf.useScreenMenuBar", "true");
            System.setProperty("apple.awt.showGrowBox", "true");

            javax.swing.Icon icon = IconUtils.getIcon(BeastMain.class, "images/beast.png");

            String nameString = "BEAST " + version.getVersionString();
            String aboutString = "<html><center><p>Bayesian Evolutionary Analysis Sampling Trees<br>" +
                    "Version " + version.getVersionString() + ", 2002-2008</p>" +
                    "<p>by<br>" +
                    "Alexei J. Drummond and Andrew Rambaut</p>" +
                    "<p>Department of Computer Science, University of Auckland<br>" +
                    "<a href=\"mailto:alexei@cs.auckland.ac.nz\">alexei@cs.auckland.ac.nz</a></p>" +
                    "<p>Institute of Evolutionary Biology, University of Edinburgh<br>" +
                    "<a href=\"mailto:a.rambaut@ed.ac.uk\">a.rambaut@ed.ac.uk</a></p>" +
                    "<p><a href=\"http://beast.bio.ed.ac.uk/\">http://beast.bio.ed.ac.uk/</a></p>" +
                    "<p>Source code distributed under the GNU LGPL:<br>" +
                    "<a href=\"http://code.google.com/p/beast-mcmc/\">http://code.google.com/p/beast-mcmc/</a></p>" +
                    "<p>Model contributions from:<br>" +
                    "Erik Bloomquist, Roald Forsberg, Joseph Heled, Gerton Lunter, Sidney Markowitz, " +
                    "Vladimir Minin, Oliver Pybus, Marc Suchard, Jen Tom</p>" +
                    "<p>Thanks to Korbinian Strimmer for use of his code</p>" +
                    "</center></html>";

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
            inputFile = Utils.getLoadFile("BEAST " + version.getVersionString() + " - Select XML input file");
        }

        if (inputFile != null && working) {
            System.setProperty("user.dir", inputFile.getParent());
        }

        printTitle();
        printHeader();

        System.out.println();
        System.out.println("Random number seed: " + seed);
        System.out.println();

        new BeastMain(inputFile, consoleApp, maxErrorCount, verbose, strictXML);
    }
}


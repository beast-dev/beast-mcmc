/*
 * BeastMain.java
 *
 * Copyright (c) 2002-2012 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

import beagle.BeagleFlag;
import beagle.BeagleInfo;
import dr.app.plugin.Plugin;
import dr.app.plugin.PluginLoader;
import dr.app.util.Arguments;
import dr.app.util.Utils;
import dr.math.MathUtils;
import dr.util.ErrorLogHandler;
import dr.util.MessageLogHandler;
import dr.util.Version;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParser;
import jam.util.IconUtils;

import javax.swing.*;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
import java.util.logging.*;

public class BeastMain {

    private final static Version version = new BeastVersion();

    static class BeastConsoleApp extends jam.console.ConsoleApplication {
        XMLParser parser = null;

        public BeastConsoleApp(String nameString, String aboutString, javax.swing.Icon icon) throws IOException {
            super(nameString, aboutString, icon, false);
            getDefaultFrame().setDefaultCloseOperation(WindowConstants.DO_NOTHING_ON_CLOSE);
        }

        public void doStop() {
            Iterator iter = parser.getThreads();
            while (iter.hasNext()) {
                Thread thread = (Thread) iter.next();
                thread.stop(); // http://java.sun.com/j2se/1.5.0/docs/guide/misc/threadPrimitiveDeprecation.html
            }
        }

        public void setTitle(String title) {
            getDefaultFrame().setTitle(title);
        }
    }

    public BeastMain(File inputFile, BeastConsoleApp consoleApp, int maxErrorCount, final boolean verbose,
                     boolean parserWarning, boolean strictXML, List<String> additionalParsers) {

        if (inputFile == null) {
            System.err.println();
            System.err.println("Error: no input file specified");
            System.exit(1);
        }

        String fileName = inputFile.getName();

        final Logger infoLogger = Logger.getLogger("dr.app.beast");
        try {

            FileReader fileReader = new FileReader(inputFile);

            XMLParser parser = new BeastParser(new String[]{fileName}, additionalParsers, verbose, parserWarning, strictXML);

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

//            // Add a handler to handle warnings and errors. This is a ConsoleHandler
//            // so the messages will go to StdErr..
//            handler = new ConsoleHandler();
//            handler.setFilter(new Filter() {
//                public boolean isLoggable(LogRecord record) {
//                    if (verbose) {
//                        return record.getLevel().intValue() >= Level.WARNING.intValue();
//                    } else {
//                        return record.getLevel().intValue() >= Level.SEVERE.intValue();
//                    }
//                }
//            });
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

            for (String pluginName : PluginLoader.getAvailablePlugins()) {
                Plugin plugin = PluginLoader.loadPlugin(pluginName);
                if (plugin != null) {
                    Set<XMLObjectParser> parserSet = plugin.getParsers();
                    for (XMLObjectParser pluginParser : parserSet) {
                        parser.addXMLObjectParser(pluginParser);
                    }
                }
            }

            parser.parse(fileReader, true);

        } catch (java.io.IOException ioe) {
            infoLogger.severe("File error: " + ioe.getMessage());
            System.exit(1);
        } catch (org.xml.sax.SAXParseException spe) {
            if (spe.getMessage() != null && spe.getMessage().equals("Content is not allowed in prolog")) {
                infoLogger.severe("Parsing error - the input file, " + fileName + ", is not a valid XML file.");
            } else {
                infoLogger.severe("Error running file: " + fileName);
                infoLogger.severe("Parsing error - poorly formed XML (possibly not an XML file):\n" +
                        spe.getMessage());
            }
            System.exit(1);
        } catch (org.w3c.dom.DOMException dome) {
            infoLogger.severe("Error running file: " + fileName);
            infoLogger.severe("Parsing error - poorly formed XML:\n" +
                    dome.getMessage());
            System.exit(1);
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
            System.exit(1);
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
            System.exit(1);
        } catch (Exception ex) {
            infoLogger.warning("Error running file: " + fileName);
            infoLogger.severe("Fatal exception: " + ex.getMessage());
            System.err.println("Fatal exception: " + ex.getMessage());
            ex.printStackTrace(System.err);
            System.exit(1);
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
        centreLine("BEAST " + version.getVersionString() + ", " + version.getDateString(), 60);
        centreLine("Bayesian Evolutionary Analysis Sampling Trees", 60);
        for (String creditLine : version.getCredits()) {
            centreLine(creditLine, 60);
        }
        System.out.println();

    }

    public static void printUsage(Arguments arguments) {

        arguments.printUsage("beast", "[<input-file-name>]");
        System.out.println();
        System.out.println("  Example: beast test.xml");
        System.out.println("  Example: beast -window test.xml");
        System.out.println("  Example: beast -help");
        System.out.println();
    }

    //Main method
    public static void main(String[] args) throws java.io.IOException {

        // There is a major issue with languages that use the comma as a decimal separator.
        // To ensure compatibility between programs in the package, enforce the US locale.
        Locale.setDefault(Locale.US);

        Arguments arguments = new Arguments(
                new Arguments.Option[]{

                        new Arguments.Option("verbose", "Give verbose XML parsing messages"),
                        new Arguments.Option("warnings", "Show warning messages about BEAST XML file"),
                        new Arguments.Option("strict", "Fail on non-conforming BEAST XML file"),
                        new Arguments.Option("window", "Provide a console window"),
                        new Arguments.Option("options", "Display an options dialog"),
                        new Arguments.Option("working", "Change working directory to input file's directory"),
                        new Arguments.LongOption("seed", "Specify a random number generator seed"),
                        new Arguments.StringOption("prefix", "PREFIX", "Specify a prefix for all output log filenames"),
                        new Arguments.Option("overwrite", "Allow overwriting of log files"),
                        new Arguments.IntegerOption("errors", "Specify maximum number of numerical errors before stopping"),
                        // new Arguments.Option("logops", "hack: log ops to stderr"),
//                        new Arguments.IntegerOption("otfops", "experimental: on the fly op weigths. recompute frequency" +
//                                "in number of states."),
                        new Arguments.IntegerOption("threads", "The number of computational threads to use (default auto)"),
                        new Arguments.Option("java", "Use Java only, no native implementations"),
                        new Arguments.Option("beagle", "Use beagle library if available"),
                        new Arguments.Option("beagle_info", "BEAGLE: show information on available resources"),
                        new Arguments.StringOption("beagle_order", "order", "BEAGLE: set order of resource use"),
                        new Arguments.IntegerOption("beagle_instances", "BEAGLE: divide site patterns amongst instances"),
                        new Arguments.Option("beagle_CPU", "BEAGLE: use CPU instance"),
                        new Arguments.Option("beagle_GPU", "BEAGLE: use GPU instance if available"),
                        new Arguments.Option("beagle_SSE", "BEAGLE: use SSE extensions if available"),
                        new Arguments.Option("beagle_single", "BEAGLE: use single precision if available"),
                        new Arguments.Option("beagle_double", "BEAGLE: use double precision if available"),
                        new Arguments.StringOption("beagle_scaling", new String[]{"default", "dynamic", "delayed", "always", "none"},
                                false, "BEAGLE: specify scaling scheme to use"),
                        new Arguments.IntegerOption("beagle_rescale", "BEAGLE: frequency of rescaling (dynamic scaling only)"),
                        new Arguments.Option("help", "Print this information and stop"),
                });

        int argumentCount = 0;

        try {
            argumentCount = arguments.parseArguments(args);
        } catch (Arguments.ArgumentException ae) {
            System.out.println();
            System.out.println(ae.getMessage());
            System.out.println();
            printUsage(arguments);
            System.exit(1);
        }

        if (arguments.hasOption("help")) {
            printUsage(arguments);
            System.exit(0);
        }

        List<String> additionalParsers = new ArrayList<String>();

        final boolean verbose = arguments.hasOption("verbose");
        final boolean parserWarning = arguments.hasOption("warnings"); // if dev, then auto turn on, otherwise default to turn off
        final boolean strictXML = arguments.hasOption("strict");
        final boolean window = arguments.hasOption("window");
        final boolean options = arguments.hasOption("options") || (argumentCount == 0);
        final boolean working = arguments.hasOption("working");
        String fileNamePrefix = null;
        boolean allowOverwrite = arguments.hasOption("overwrite");

        long seed = MathUtils.getSeed();
        boolean useJava = false;

        int threadCount = -1;

        if (arguments.hasOption("java")) {
            useJava = true;
        }

        if (arguments.hasOption("prefix")) {
            fileNamePrefix = arguments.getStringOption("prefix");
        }

        long beagleFlags = 0;

        // if any beagle flag is specified then use beagle...
        boolean useBeagle = arguments.hasOption("beagle") ||
                arguments.hasOption("beagle_CPU") ||
                arguments.hasOption("beagle_GPU") ||
                arguments.hasOption("beagle_SSE") ||
                arguments.hasOption("beagle_double") ||
                arguments.hasOption("beagle_single") ||
                arguments.hasOption("beagle_order") ||
                arguments.hasOption("beagle_scaling") ||
                arguments.hasOption("beagle_rescale") ||
                arguments.hasOption("beagle_instances");

        boolean beagleShowInfo = arguments.hasOption("beagle_info");

        if (arguments.hasOption("beagle_CPU")) {
            beagleFlags |= BeagleFlag.PROCESSOR_CPU.getMask();
        }
        if (arguments.hasOption("beagle_GPU")) {
            beagleFlags |= BeagleFlag.PROCESSOR_GPU.getMask();
        }
        if (arguments.hasOption("beagle_SSE")) {
            beagleFlags |= BeagleFlag.PROCESSOR_CPU.getMask();
            beagleFlags |= BeagleFlag.VECTOR_SSE.getMask();
        }
        if (arguments.hasOption("beagle_double")) {
            beagleFlags |= BeagleFlag.PRECISION_DOUBLE.getMask();
        }
        if (arguments.hasOption("beagle_single")) {
            beagleFlags |= BeagleFlag.PRECISION_SINGLE.getMask();
        }

        if (arguments.hasOption("beagle_order")) {
            System.setProperty("beagle.resource.order", arguments.getStringOption("beagle_order"));
        }

        if (arguments.hasOption("beagle_instances")) {
            System.setProperty("beagle.instance.count", Integer.toString(arguments.getIntegerOption("beagle_instances")));
        }

        if (arguments.hasOption("beagle_scaling")) {
            System.setProperty("beagle.scaling", arguments.getStringOption("beagle_scaling"));
        }

        if (arguments.hasOption("beagle_rescale")) {
            System.setProperty("beagle.rescale", Integer.toString(arguments.getIntegerOption("beagle_rescale")));
        }

        if (arguments.hasOption("threads")) {
            // threadCount defaults to -1 unless the user specifies an option
            threadCount = arguments.getIntegerOption("threads");
            if (threadCount < 0) {
                printTitle();
                System.err.println("The the number of threads should be >= 0");
                System.exit(1);
            }
        }

        if (arguments.hasOption("seed")) {
            seed = arguments.getLongOption("seed");
            if (seed <= 0) {
                printTitle();
                System.err.println("The random number seed should be > 0");
                System.exit(1);
            }
        }

        int maxErrorCount = 0;
        if (arguments.hasOption("errors")) {
            maxErrorCount = arguments.getIntegerOption("errors");
            if (maxErrorCount < 0) {
                maxErrorCount = 0;
            }
        }

        BeastConsoleApp consoleApp = null;

        String nameString = "BEAST " + version.getVersionString();

        if (window) {
            System.setProperty("com.apple.macos.useScreenMenuBar", "true");
            System.setProperty("apple.laf.useScreenMenuBar", "true");
            System.setProperty("apple.awt.showGrowBox", "true");

            javax.swing.Icon icon = IconUtils.getIcon(BeastMain.class, "images/beast.png");

            String aboutString = "<html><div style=\"font-family:sans-serif;\"><center>" +
                    "<div style=\"font-size:12;\"><p>Bayesian Evolutionary Analysis Sampling Trees<br>" +
                    "Version " + version.getVersionString() + ", " + version.getDateString() + "</p>" +
                    version.getHTMLCredits() +
                    "</div></center></div></html>";

            consoleApp = new BeastConsoleApp(nameString, aboutString, icon);
        }

        printTitle();

        File inputFile = null;

        if (options) {

            String titleString = "<html><center><p>Bayesian Evolutionary Analysis Sampling Trees<br>" +
                    "Version " + version.getVersionString() + ", " + version.getDateString() + "</p></center></html>";
            javax.swing.Icon icon = IconUtils.getIcon(BeastMain.class, "images/beast.png");

            BeastDialog dialog = new BeastDialog(new JFrame(), titleString, icon);

            dialog.setAllowOverwrite(allowOverwrite);
            dialog.setSeed(seed);

            dialog.setUseBeagle(useBeagle);

            if (BeagleFlag.PROCESSOR_GPU.isSet(beagleFlags)) {
                dialog.setPreferBeagleGPU();
            }

            dialog.setPreferBeagleSSE(BeagleFlag.VECTOR_SSE.isSet(beagleFlags));

            if (BeagleFlag.PRECISION_SINGLE.isSet(beagleFlags)) {
                dialog.setPreferBeagleSingle();
            }

            if (!dialog.showDialog(nameString)) {
                return;
            }

            if (dialog.allowOverwrite()) {
                allowOverwrite = true;
            }

            seed = dialog.getSeed();
            threadCount = dialog.getThreadPoolSize();

            useBeagle = dialog.useBeagle();
            if (useBeagle) {
                beagleShowInfo = dialog.showBeagleInfo();
                if (dialog.preferBeagleCPU()) {
                    beagleFlags |= BeagleFlag.PROCESSOR_CPU.getMask();
                }
                if (dialog.preferBeagleSSE()) {
                    beagleFlags |= BeagleFlag.VECTOR_SSE.getMask();
                }
                if (dialog.preferBeagleGPU()) {
                    beagleFlags |= BeagleFlag.PROCESSOR_GPU.getMask();
                }
                if (dialog.preferBeagleDouble()) {
                    beagleFlags |= BeagleFlag.PRECISION_DOUBLE.getMask();
                }
                if (dialog.preferBeagleSingle()) {
                    beagleFlags |= BeagleFlag.PRECISION_SINGLE.getMask();
                }
                System.setProperty("beagle.scaling", dialog.scalingScheme());
            }

            inputFile = dialog.getInputFile();
            if (!beagleShowInfo && inputFile == null) {
                System.err.println("No input file specified");
                return;
            }

        }

        if (beagleShowInfo) {
            BeagleInfo.printResourceList();
            return;
        }

        if (inputFile == null) {

            String[] args2 = arguments.getLeftoverArguments();

            if (args2.length > 1) {
                System.err.println("Unknown option: " + args2[1]);
                System.err.println();
                printUsage(arguments);
                return;
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

        if (inputFile != null && inputFile.getParent() != null && working) {
            System.setProperty("user.dir", inputFile.getParent());
        }

        if (window) {
            if (inputFile == null) {
                consoleApp.setTitle("null");
            } else {
                consoleApp.setTitle(inputFile.getName());
            }
        }

        if (useJava) {
            System.setProperty("java.only", "true");
        }

        if (fileNamePrefix != null && fileNamePrefix.trim().length() > 0) {
            System.setProperty("file.name.prefix", fileNamePrefix.trim());
        }

        if (allowOverwrite) {
            System.setProperty("log.allow.overwrite", "true");
        }

        if (useBeagle) {
            additionalParsers.add("beagle");
        }

        if (beagleFlags != 0) {
            System.setProperty("beagle.preferred.flags", Long.toString(beagleFlags));

        }

        if (threadCount >= 0) {
            System.setProperty("thread.count", String.valueOf(threadCount));
        }

        MathUtils.setSeed(seed);

        System.out.println();
        System.out.println("Random number seed: " + seed);
        System.out.println();

        try {
            new BeastMain(inputFile, consoleApp, maxErrorCount, verbose, parserWarning, strictXML, additionalParsers);
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



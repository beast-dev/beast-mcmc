package dr.app.beast;

import dr.app.util.Arguments;
import dr.math.MathUtils;
import mpi.MPI;

import java.io.File;

/**
 * Created by IntelliJ IDEA.
 * User: msuchard
 * Date: Jan 23, 2007
 * Time: 6:17:34 PM
 * To change this template use File | Settings | File Templates.
 */
public class BeastRemote extends BeastMain {

    public BeastRemote(File inputFile, BeastConsoleApp consoleApp, boolean verbose) {
        super(inputFile, consoleApp, verbose);
    }

    public static void main(String[] oldArgs) throws java.io.IOException {

        // First populate args from MPI.WORLD

        // String[] args = null;
        MPI.Init(oldArgs);
        int argLength = oldArgs.length - 3;
        String[] args = new String[argLength];
        for (int i = 0; i < argLength; i++)
            args[i] = oldArgs[i + 3];

        /*  for (String str : args)
            System.err.println(str);
        System.exit(1);*/


        Arguments arguments = new Arguments(
                new Arguments.Option[]{
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

/*		if (window) {
			System.setProperty("com.apple.macos.useScreenMenuBar", "true");
			System.setProperty("apple.laf.useScreenMenuBar", "true");
			System.setProperty("apple.awt.showGrowBox", "true");

			javax.swing.Icon icon = IconUtils.getIcon(BeastMain.class, "images/beast.png");

			String nameString = "BEAST " + version.getVersionString();
			String aboutString = "<html><center><p>Bayesian Evolutionary Analysis Sampling Trees<br>" +
					"Version " + version.getVersionString() + ", 2002-2006</p>" +
					"<p>by<br>" +
					"Alexei J. Drummond and Andrew Rambaut</p>" +
					"<p>Department of Computer Science, University of Auckland<br>" +
					"<a href=\"mailto:alexei@cs.auckland.ac.nz\">alexei@cs.auckland.ac.nz</a></p>" +
					"<p>Institute of Evolutionary Biology, University of Edinburgh<br>" +
					"<a href=\"mailto:a.rambaut@ed.ac.uk\">a.rambaut@ed.ac.uk</a></p>" +
					"<p><a href=\"http://evolve.zoo.ox.ac.uk/beast/\">http://evolve.zoo.ox.ac.uk/beast/</a></p>" +
					"<p>Source code distributed under the GNU LGPL:<br>" +
					"<a href=\"http://code.google.com/p/beast-mcmc/\">http://code.google.com/p/beast-mcmc/</a></p>" +
					"<p>Additional programming by:<br>" +
					"Roald Forsberg, Gerton Lunter, Sidney Markowitz, Oliver Pybus</p>" +
					"<p>Thanks to Korbinian Strimmer for use of his code</p>" +
					"</center></html>";

			consoleApp = new BeastConsoleApp(nameString, aboutString, icon);
		}*/      // Remote can never be interactive

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

        /*	if (inputFileName == null) {
              // No input file name was given so throw up a dialog box...
              inputFile = Utils.getLoadFile("BEAST " + version.getVersionString() + " - Select XML input file");
          }*/

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

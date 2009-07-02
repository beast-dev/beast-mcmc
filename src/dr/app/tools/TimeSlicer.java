package dr.app.tools;

import dr.app.beast.BeastVersion;
import dr.app.util.Arguments;
import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.TreeImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.util.Version;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 */

public class TimeSlicer {

    public static final String sep = "\t";

    public TimeSlicer(String treeFileName, int burnin, String[] traits, double[] slices) {

        List<Tree> trees = null;

        try {
            trees = importTrees(treeFileName, burnin);
        } catch (IOException e) {
            System.err.println("Error reading file: " + treeFileName);
            System.exit(-1);
        } catch (Importer.ImportException e) {
            System.err.println("Error parsing trees in file: " + treeFileName);
            System.exit(-1);
        }
        if (trees == null || trees.size() == 0) {
            System.err.println("No trees read from file: " + treeFileName);
            System.exit(-1);
        }

        run(trees, traits, slices);

        outputHeader(traits);

        if (slices == null)
            outputSlice(0,1);
        else {
            for(int i=0; i<slices.length; i++)
                outputSlice(i,slices[i]);
        
        }

    }

    private void outputHeader(String[] traits) {
        StringBuffer sb = new StringBuffer("slice");
        for(int i=0; i<traits.length; i++) {
            sb.append(sep);
            sb.append(traits[i].toString());
        }
        sb.append("\n");
        resultsStream.print(sb);
    }

    private List<Tree> importTrees(String treeFileName, int burnin) throws IOException, Importer.ImportException {

        int totalTrees = 10000;
//        int totalTreesUsed = 0;

        progressStream.println("Reading trees (bar assumes 10,000 trees)...");
        progressStream.println("0              25             50             75            100");
        progressStream.println("|--------------|--------------|--------------|--------------|");

        int stepSize = totalTrees / 60;
        if (stepSize < 1) stepSize = 1;


        List<Tree> treeList = new ArrayList<Tree>();
        BufferedReader reader1 = new BufferedReader(new FileReader(treeFileName));

        String line1 = reader1.readLine();
        TreeImporter importer1 = null;
        if (line1.toUpperCase().startsWith("#NEXUS")) {
            importer1 = new NexusImporter(new FileReader(treeFileName));
        } else {
            importer1 = new NewickImporter(new FileReader(treeFileName));
        }
        totalTrees = 0;
        while (importer1.hasTree()) {
            Tree treeTime = importer1.importNextTree();

            if (totalTrees > burnin)
                treeList.add(treeTime);

            if (totalTrees > 0 && totalTrees % stepSize == 0) {
                progressStream.print("*");
                progressStream.flush();
            }
            totalTrees++;
        }
        return treeList;
    }

    class Trait {

        Trait(Object obj) {
            this.obj = obj;
        }

        private Object obj;

        public String toString() { return obj.toString(); }

    }

    private List<List<List<Trait>>> values;

    private void outputSlice(int slice, double sliceValue) {

        List<List<Trait>> thisSlice = values.get(slice);
        int traitCount = thisSlice.size();
        int valueCount = thisSlice.get(0).size();

        StringBuffer sb = new StringBuffer();

        for(int v=0; v<valueCount; v++) {
            sb.append(sliceValue);
            for(int t=0; t<traitCount; t++) {
                sb.append(sep);
                sb.append(thisSlice.get(t).get(v));
            }
            sb.append("\n");
        }

       resultsStream.print(sb);            
    }

    private void run(List<Tree> trees, String[] traits, double[] slices) {

        int traitCount = traits.length;
        int sliceCount = 1;
        boolean doSlices = false;

        if (slices != null) {
            sliceCount = slices.length;
            doSlices = true;
        }

        values = new ArrayList<List<List<Trait>>>(sliceCount);
        for (int i = 0; i < sliceCount; i++) {
            List<List<Trait>> thisSlice = new ArrayList<List<Trait>>(traitCount);
            values.add(thisSlice);
            for (int j = 0; j < traitCount; j++) {
                List<Trait> thisTraitSlice = new ArrayList<Trait>();
                thisSlice.add(thisTraitSlice);
            }
        }

        for (Tree treeTime : trees) {

            for (int x = 0; x < treeTime.getNodeCount(); x++) {

                NodeRef node = treeTime.getNode(x);

                if (!(treeTime.isRoot(node))) {

                    double nodeHeight = treeTime.getNodeHeight(node);
                    double parentHeight = treeTime.getNodeHeight(treeTime.getParent(node));

                    for (int i = 0; i < sliceCount; i++) {

                        if (!doSlices ||
                                (slices[i] > nodeHeight && slices[i] < parentHeight)
                                ) {

                            List<List<Trait>> thisSlice = values.get(i);
                            for (int j = 0; j < traitCount; j++) {
                                
                                List<Trait> thisTraitSlice = thisSlice.get(j);
                                Object trait = treeTime.getNodeAttribute(node, traits[j]);
                                if (trait == null) {
                                    System.err.println("Trait '"+traits[j]+"' not found on branch.");
                                    System.exit(-1);
                                }
                                thisTraitSlice.add(new Trait(trait));
                            }
                        }
                    }
                }
            }
        }
    }

    // Messages to stderr, output to stdout
    private static PrintStream progressStream = System.err;
    private static PrintStream resultsStream = System.out;


    private final static Version version = new BeastVersion();

    private static final String commandName = "treeslicer";


    public static void printUsage(Arguments arguments) {

        arguments.printUsage(commandName, "<input-file-name> [<output-file-name>]");
        progressStream.println();
        progressStream.println("  Example: " + commandName + " test.trees out.txt");
        progressStream.println();
    }

    public static void centreLine(String line, int pageWidth) {
        int n = pageWidth - line.length();
        int n1 = n / 2;
        for (int i = 0; i < n1; i++) {
            progressStream.print(" ");
        }
        progressStream.println(line);
    }

    public static void printTitle() {
        progressStream.println();
        centreLine("TimeSlicer " + version.getVersionString() + ", " + version.getDateString(), 60);
        centreLine("MCMC Output analysis", 60);
        centreLine("by", 60);
        centreLine("Marc A. Suchard, Philippe Lemey,", 60);
        centreLine("Alexei J. Drummond and Andrew Rambaut", 60);
        progressStream.println();
        centreLine("Department of Biomathematics", 60);
        centreLine("University of California, Los Angeles", 60);
        centreLine("msuchard@ucla.edu", 60);
        progressStream.println();
        centreLine("DEPARTMENT", 60);
        centreLine("UNIVERSITY", 60);
        centreLine("EMAIL", 60);
        progressStream.println();
        centreLine("Department of Computer Science", 60);
        centreLine("University of Auckland", 60);
        centreLine("alexei@cs.auckland.ac.nz", 60);
        progressStream.println();
        centreLine("Institute of Evolutionary Biology", 60);
        centreLine("University of Edinburgh", 60);
        centreLine("a.rambaut@ed.ac.uk", 60);
        progressStream.println();
        progressStream.println();
    }


    //Main method
    public static void main(String[] args) throws IOException {


        String inputFileName = null;
        String outputFileName = null;

//        if (args.length == 0) {
//          // TODO Make flash GUI
//        }

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.IntegerOption("burnin", "the number of states to be considered as 'burn-in' [default = 0]"),
                        new Arguments.StringOption("trait", "trait_name", "specifies an attribute to use to create a density map [default = location.angle]"),
                        new Arguments.StringOption("trait2", "trait2_name", "specifies an attribute to use to create a density map [default = null]"),
                        new Arguments.RealArrayOption("slice", 0, Double.MAX_VALUE, "list of times to perform summary"), new Arguments.Option("help", "option to print this message")
                });

        try {
            arguments.parseArguments(args);
        } catch (Arguments.ArgumentException ae) {
            progressStream.println(ae);
            printUsage(arguments);
            System.exit(1);
        }

        if (arguments.hasOption("help")) {
            printUsage(arguments);
            System.exit(0);
        }

        int burnin = -1;
        if (arguments.hasOption("burnin")) {
            burnin = arguments.getIntegerOption("burnin");
        }

//        String[] traitNames = arguments.getStringArrayOption("trait");

//        System.err.println("length = "+traitNames.length);

        String[] traitNames;
        String traitName1 = arguments.getStringOption("trait2");
        if (traitName1 == null)
            traitNames = new String[1];
        else {
            traitNames = new String[2];
            traitNames[1] = traitName1;
        }

        traitNames[0] = arguments.getStringOption("trait"); // TODO Read in a list of traits
        if(traitNames[0] == null)
            traitNames[0] = "location.rate";

        double[] slices = arguments.getRealArrayOption("slice");

        final String[] args2 = arguments.getLeftoverArguments();

        switch (args2.length) {
            case 0:
                printUsage(arguments);
                System.exit(1);
            case 2:
                outputFileName = args2[1];
                // fall to
            case 1:
                inputFileName = args2[0];
                break;
            default: {
                System.err.println("Unknown option: " + args2[2]);
                System.err.println();
                printUsage(arguments);
                System.exit(1);
            }
        }

        new TimeSlicer(inputFileName, burnin, traitNames, slices);

        System.exit(0);
    }
}

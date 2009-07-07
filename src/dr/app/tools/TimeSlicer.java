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
import dr.math.distributions.MultivariateNormalDistribution;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 */

public class TimeSlicer {

    public static final String sep = "\t";

    public static final String PRECISION_STRING = "precision";
    public static final String RATE_STRING = "rate";

    public TimeSlicer(String treeFileName, String outFileName, int burnin, String[] traits, double[] slices, boolean impute, boolean trueNoise) {

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

        run(trees, traits, slices, impute, trueNoise);

        resultsStream = System.out;

        if (outFileName != null) {
            try {
            resultsStream = new PrintStream(new File(outFileName));
            } catch (IOException e) {
                System.err.println("Error opening file: "+outFileName);
                System.exit(-1);
            }
        }

        outputHeader(traits);

        if (slices == null)
            outputSlice(0,Double.NaN);
        else {
            for(int i=0; i<slices.length; i++)
                outputSlice(i,slices[i]);
        
        }

    }

    private void outputHeader(String[] traits) {
        StringBuffer sb = new StringBuffer("slice");
        for(int i=0; i<traits.length; i++) {
            // Load first value to check dimensionality
            Trait trait = values.get(0).get(i).get(0);
            if (trait.isMultivariate()) {
                int dim = trait.getDim();
                for(int j=1; j<=dim; j++)
                    sb.append(sep).append(traits[i]).append(j);
            } else
                sb.append(sep).append(traits[i]);
        }
        sb.append("\n");
        resultsStream.print(sb);
    }

    private List<Tree> importTrees(String treeFileName, int burnin) throws IOException, Importer.ImportException {

        int totalTrees = 10000;

        progressStream.println("Reading trees (bar assumes 10,000 trees)...");
        progressStream.println("0              25             50             75            100");
        progressStream.println("|--------------|--------------|--------------|--------------|");

        int stepSize = totalTrees / 60;
        if (stepSize < 1) stepSize = 1;

        List<Tree> treeList = new ArrayList<Tree>();
        BufferedReader reader1 = new BufferedReader(new FileReader(treeFileName));

        String line1 = reader1.readLine();
        TreeImporter importer1;
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
            if (obj instanceof Object[]) {
                isMultivariate = true;
                array = (Object[])obj;
            }
        }

        public boolean isMultivariate() { return isMultivariate; }

        public boolean isNumber() {
            if (!isMultivariate)
                return (obj instanceof Double);
            return (array[0] instanceof Double);
        }

        public int getDim() {
            if (isMultivariate) {
                return array.length;
            }
            return 1;
        }

        public double[] getValue() {
            int dim = getDim();
            double[] result = new double[dim];
            for(int i=0; i<dim; i++)
                result[i] = (Double)array[i];
            return result;
        }

        private Object obj;
        private Object[] array;
        private boolean isMultivariate = false;

        public String toString() {
            if (!isMultivariate)
                return obj.toString();
            StringBuffer sb = new StringBuffer(array[0].toString());
            for(int i=1; i<array.length; i++)
                sb.append(sep).append(array[i]);
            return sb.toString();
        }
    }

    private List<List<List<Trait>>> values;

    private void outputSlice(int slice, double sliceValue) {

        List<List<Trait>> thisSlice = values.get(slice);
        int traitCount = thisSlice.size();
        int valueCount = thisSlice.get(0).size();

        StringBuffer sb = new StringBuffer();

        for(int v=0; v<valueCount; v++) {
            if (Double.isNaN(sliceValue))
                sb.append("All");
            else
                sb.append(sliceValue);
            for(int t=0; t<traitCount; t++) {
                sb.append(sep);
                sb.append(thisSlice.get(t).get(v));
            }
            sb.append("\n");
        }

       resultsStream.print(sb);            
    }

    private void run(List<Tree> trees, String[] traits, double[] slices, boolean impute, boolean trueNoise) {

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

            double[][] precision = null;

            if (impute) {
                Object o = treeTime.getAttribute("precision");
                if (o != null) {
                    Object[] array = (Object[])o;
                    int dim = (int) Math.sqrt(1+8*array.length) / 2;
                    precision = new double[dim][dim];
                    int c = 0;
                    for(int i=0; i<dim; i++) {
                        for(int j=i; j<dim; j++) {
                            precision[j][i] = precision[i][j] = (Double)array[c++];
                        }
                    }
                }              
            }

            for (int x = 0; x < treeTime.getNodeCount(); x++) {

                NodeRef node = treeTime.getNode(x);

                if (!(treeTime.isRoot(node))) {

                    double nodeHeight = treeTime.getNodeHeight(node);
                    double parentHeight = treeTime.getNodeHeight(treeTime.getParent(node));

                    for (int i = 0; i < sliceCount; i++) {

                        if (!doSlices ||
                                (slices[i] >= nodeHeight && slices[i] < parentHeight)
                                ) {

                            List<List<Trait>> thisSlice = values.get(i);
                            for (int j = 0; j < traitCount; j++) {
                                
                                List<Trait> thisTraitSlice = thisSlice.get(j);
                                Object tmpTrait = treeTime.getNodeAttribute(node, traits[j]);
                                if (tmpTrait == null) {
                                    System.err.println("Trait '"+traits[j]+"' not found on branch.");
                                    System.exit(-1);
                                }
                                Trait trait = new Trait(tmpTrait);
                                if (impute) {
                                    Double rateAttribute = (Double) treeTime.getNodeAttribute(node,RATE_STRING);
                                    double rate = 1.0;
                                    if (rateAttribute != null) {
                                        rate = rateAttribute;
                                        if (outputRateWarning) {
                                            progressStream.println("Warning: using a rate attribute during imputation!");
                                            outputRateWarning = false;
                                        }
                                    }
                                    trait = imputeValue(trait,new Trait(treeTime.getNodeAttribute(treeTime.getParent(node),traits[j])),
                                            slices[i],nodeHeight, parentHeight,precision, rate, trueNoise);
                                }
                                thisTraitSlice.add(trait);
                            }
                        }
                    }
                }
            }
        }
    }

    private boolean outputRateWarning = true;

    private Trait imputeValue(Trait nodeTrait, Trait parentTrait, double time, double nodeHeight, double parentHeight, double[][] precision, double rate, boolean trueNoise) {
        if (!nodeTrait.isNumber()) {
            System.err.println("Can only impute numbers!");
            System.exit(-1);
        }

        int dim = precision.length;
        double[] nodeValue = nodeTrait.getValue();
        double[] parentValue = parentTrait.getValue();

        final double timeTotal = parentHeight - nodeHeight;
        final double timeChild = (time - nodeHeight);
        final double timeParent = (parentHeight - time);
        final double weightTotal = 1.0 / timeChild + 1.0 / timeParent;

        if (timeChild == 0)
            return nodeTrait;

        if (timeParent == 0)
            return parentTrait;

        // Find mean value, weighted average
        double[] mean = new double[dim];
        double[][] scaledPrecision = new double[dim][dim];

        for(int i=0; i<dim; i++) {
            mean[i] = (nodeValue[i] / timeChild + parentValue[i] / timeParent) / weightTotal;
            if (trueNoise) {
                for(int j=i; j<dim; j++)
                    scaledPrecision[j][i] = scaledPrecision[i][j] = precision[i][j] / timeTotal / rate;
            }
        }

        if (trueNoise)
            mean = MultivariateNormalDistribution.nextMultivariateNormalPrecision(mean, precision);
        Object[] result = new Object[dim];
        for(int i=0; i<dim; i++)
            result[i] = mean[i];
        return new Trait(result);
    }

    // Messages to stderr, output to stdout
    private static PrintStream progressStream = System.err;
    private PrintStream resultsStream;


    private final static Version version = new BeastVersion();

    private static final String commandName = "timeslicer";


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
        centreLine("Rega Institute for Medical Research", 60);
        centreLine("Katholieke Unversiteit Leuven", 60);
        centreLine("philippe.lemey@gmail.com", 60);
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

    private static double[] parseVariableLengthDoubleArray(String inString) throws Arguments.ArgumentException {

        List<Double> returnList = new ArrayList<Double>();
        StringTokenizer st = new StringTokenizer(inString,",");
        while(st.hasMoreTokens()) {
            try {
                returnList.add(Double.parseDouble(st.nextToken()));
            } catch (NumberFormatException e) {
                throw new Arguments.ArgumentException();
            }

        }

        if (returnList.size()>0) {
            double[] doubleArray = new double[returnList.size()];
            for(int i=0; i<doubleArray.length; i++)
                doubleArray[i] = returnList.get(i);
            return doubleArray;
        }
        return null;
    }

     private static String[] parseVariableLengthStringArray(String inString) {

        List<String> returnList = new ArrayList<String>();
        StringTokenizer st = new StringTokenizer(inString,",");
        while(st.hasMoreTokens()) {
            returnList.add(st.nextToken());           
        }

        if (returnList.size()>0) {
            String[] stringArray = new String[returnList.size()];
            stringArray = returnList.toArray(stringArray);
            return stringArray;
        }
        return null;
    }

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
                        new Arguments.StringOption("trait", "trait_name", "specifies an attribute-list to use to create a density map [default = location.angle]"),
                        new Arguments.StringOption("slice","time","specifies an slice time-list [default=none]"),
                        new Arguments.Option("help", "option to print this message"),
                        new Arguments.StringOption("noise", new String[]{"false","true"}, false,
                                "add true noise [default = true])"),
                        new Arguments.StringOption("impute", new String[]{"false", "true"}, false,
                                "impute trait at time-slice [default = false])"),
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

        String[] traitNames = null;
        double[] sliceTimes = null;

        try {

            String traitString = arguments.getStringOption("trait");
            if (traitString != null) {
                traitNames = parseVariableLengthStringArray(traitString);
            }
            if (traitNames == null) {
                traitNames = new String[1];
                traitNames[0] = "location.rate";
            }

            String sliceString = arguments.getStringOption("slice");
            if (sliceString != null) {
                sliceTimes = parseVariableLengthDoubleArray(sliceString);
            }
                       
        } catch (Arguments.ArgumentException e) {
            progressStream.println(e);
            printUsage(arguments);
            System.exit(-1);
        }

        String imputeString = arguments.getStringOption("impute");
        boolean impute = false;
        if (imputeString != null && imputeString.compareToIgnoreCase("true") == 0)
            impute = true;

        String noiseString = arguments.getStringOption("noise");
        boolean trueNoise = true;
        if (noiseString != null && noiseString.compareToIgnoreCase("false") == 0)
            trueNoise = false;

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

        new TimeSlicer(inputFileName, outputFileName, burnin, traitNames, sliceTimes, impute, trueNoise);

        System.exit(0);
    }
}

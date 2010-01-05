package dr.app.tools;


import dr.app.util.Arguments;
import dr.geo.math.SphericalPolarCoordinates;
import dr.stats.DiscreteStatistics;
import org.jdom.Element;
import org.jdom.output.Format;
import org.jdom.output.XMLOutputter;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

/**
 * Created by IntelliJ IDEA.
 * User: phil
 * Date: Sep 25, 2009
 * Time: 7:53:00 AM
 * To change this template use File | Settings | File Templates.
 */
public class DiscreteRatePriorGenerator {

    public static final String HELP = "help";
    public static final String COORDINATES = "coordinates";
    public static final String DENSITIES = "densities";
    public static final String GENERICS = "generics";
    public static final String FORMAT = "format";

    public DiscreteRatePriorGenerator(String[] locations, Double[] latitudes, Double[] longitudes, Double[] densities) {

        if (locations == null) {
            //System.out.println(locations[0]);
            System.err.println("no locations specified!");
        } else {
            //this.locations = locations;
        }

        //this.latitudes = latitudes;
        //this.longitudes = longitudes;
        this.densities = densities;

        if ((latitudes == null)||(longitudes == null)) {
            progressStream.println("no latitudes or longitudes specified!");
        } else {
            distances = getUpperTriangleDistanceMatrix(latitudes,longitudes);
        }

    }
    // for the time being, locations, latitudes and longitudes are not equired
    //private String[] locations;
    //private Double[] latitudes;
    //private Double[] longitudes;
    private final Double[] densities;
    private double[] distances;

    private double[] getUpperTriangleDistanceMatrix(Double[] latitudes, Double[] longitudes) {

        double[] distances = new double[(latitudes.length * (latitudes.length - 1))/2];

        int distanceCounter = 0;
        int pairwiseCounter1 = 0;
        for (int c = 0; c < latitudes.length; c++) {
            pairwiseCounter1 ++;
            for (int d = pairwiseCounter1; d < latitudes.length; d++) {
                distances[distanceCounter] = getKilometerGreatCircleDistance(latitudes[c],longitudes[c],latitudes[d],longitudes[d]);
                distanceCounter++;
            }
        }
        return distances;
    }

    private static double getKilometerGreatCircleDistance(double lat1, double long1, double lat2, double long2) {
         SphericalPolarCoordinates coord1 = new SphericalPolarCoordinates(lat1, long1);
         SphericalPolarCoordinates coord2 = new SphericalPolarCoordinates(lat2, long2);
         return (coord1.distance(coord2));
    }

    enum OutputFormat {
        TAB,
        SPACE,
        XML
    }

    public void output(String outputFileName, OutputFormat outputFormat) {

        resultsStream = System.out;

        if (outputFileName != null) {
            try {
            resultsStream = new PrintStream(new File(outputFileName));
            } catch (IOException e) {
                System.err.println("Error opening file: "+outputFileName);
                System.exit(1);
            }
        }
        if (distances != null) {
            outputStringLine("reversible priors: distances",outputFormat);
            resultsStream.print("\r");
            outputArray(distances, outputFormat, false);
            resultsStream.print("\r");
            outputStringLine("reversible priors: normalized inverse distances",outputFormat);
            resultsStream.print("\r");
            outputArray(transform(distances,true, true, false), outputFormat, false);
            resultsStream.print("\r");
            outputStringLine("reversible priors: normalized inverse log distances",outputFormat);
            resultsStream.print("\r");
            outputArray(transform(distances,true, true, true), outputFormat, false);
            resultsStream.print("\r");
        }
        if (densities != null) {
            outputStringLine("reversible priors: densities",outputFormat);
            resultsStream.print("\r");
            outputArray(transform(getUpperTrianglePairwiseProductMatrix(densities),false,false,false), outputFormat, false);
            resultsStream.print("\r");
            outputStringLine("reversible priors: normalized densities",outputFormat);
            resultsStream.print("\r");
            outputArray(transform(getUpperTrianglePairwiseProductMatrix(densities),false,true,false), outputFormat, false);
            resultsStream.print("\r");
            outputStringLine("reversible priors: normalized log densities",outputFormat);
            resultsStream.print("\r");
            outputArray(transform(getUpperTrianglePairwiseProductMatrix(densities),false,true,true), outputFormat, false);
            resultsStream.print("\r");

            if (distances != null) {
                outputStringLine("reversible priors: product of normalized densities divided by normalized distances",outputFormat);
                resultsStream.print("\r");
                outputArray(productOfArrays(transform(getUpperTrianglePairwiseProductMatrix(densities),false,true,false),transform(distances,false, true, false),true), outputFormat, false);
                resultsStream.print("\r");
            }

        }
        if (distances != null) {
            outputStringLine("nonreversible priors: distances",outputFormat);
            resultsStream.print("\r");
            outputArray(distances, outputFormat, true);
            resultsStream.print("\r");
            outputStringLine("nonreversible priors: normalized inverse distances",outputFormat);
            resultsStream.print("\r");
            outputArray(transform(distances,true, true, false), outputFormat, true);
            resultsStream.print("\r");
            outputStringLine("nonreversible priors: normalized inverse log distances",outputFormat);
            resultsStream.print("\r");
            outputArray(transform(distances,true, true, true), outputFormat, true);
            resultsStream.print("\r");
        }
        if (densities != null) {
            outputStringLine("nonreversible priors: densities",outputFormat);
            resultsStream.print("\r");
            outputArray(transform(getUpperTrianglePairwiseProductMatrix(densities),false,false,false), outputFormat, true);
            resultsStream.print("\r");
            outputStringLine("nonreversible priors: normalized densities",outputFormat);
            resultsStream.print("\r");
            outputArray(transform(getUpperTrianglePairwiseProductMatrix(densities),false,true,false), outputFormat, true);
            resultsStream.print("\r");
            outputStringLine("nonreversible priors: normalized log densities",outputFormat);
            resultsStream.print("\r");
            outputArray(transform(getUpperTrianglePairwiseProductMatrix(densities),false,true,true), outputFormat, true);
            resultsStream.print("\r");

            if (distances != null) {
                outputStringLine("nonreversible priors: product of normalized densities divided by normalized distances",outputFormat);
                resultsStream.print("\r");
                outputArray(productOfArrays(transform(getUpperTrianglePairwiseProductMatrix(densities),false,true,false),transform(distances,false, true, false),true), outputFormat, true);
                resultsStream.print("\r");
            }

        }
        

    }

    private double[] getUpperTrianglePairwiseProductMatrix(Double[] matrixValues) {

        double[] pairwiseProducts = new double[(matrixValues.length * (matrixValues.length - 1))/2];

        int counter = 0;
        int pairwiseCounter1 = 0;
        for(Double matrixValue : matrixValues) {
            pairwiseCounter1++;
            for(int d = pairwiseCounter1; d < matrixValues.length; d++) {
                pairwiseProducts[counter] = matrixValue * matrixValues[d];
                counter++;
            }
        }
        return pairwiseProducts;
    }

    private double[] productOfArrays(double[] inputArray1, double[] inputArray2, boolean devision) {

        if (inputArray1.length != inputArray2.length) {
            System.err.println("trying to get a product of arrays of unequals size!");
            System.exit(1);
        }

        double[] productOfArray = new double[inputArray1.length];

        for (int i = 0; i < inputArray1.length; i++) {
            if (devision) {
                productOfArray[i] = inputArray1[i]/inputArray2[i];
            } else {
            }
            productOfArray[i] = inputArray1[i]*inputArray2[i];
        }

        return productOfArray;
    }


    private double[] transform(double[] inputArray, boolean inverse, boolean normalize, boolean log) {

        double[] transformedDistances = new double[inputArray.length];

        for(int u=0; u<inputArray.length; u++) {

            double distance = inputArray[u];

            if (log) {
                distance = Math.log(distance);
            }
            if (inverse) {
                distance = 1/distance;
            }

            transformedDistances[u] = distance;
        }

        double meanDistance = DiscreteStatistics.mean(transformedDistances);
        if (!normalize) meanDistance = 1;

        for(int v=0; v<inputArray.length; v++) {
            transformedDistances[v] =  (transformedDistances[v]/meanDistance);
        }
        return transformedDistances;
    }

    //TODO: check if no other way of converting Double[] to double[]
    private double[] transform(Double[] inputArray, boolean inverse, boolean normalize, boolean log) {

        double[] transformedDistances = new double[inputArray.length];

        for(int u=0; u<inputArray.length; u++) {

            double distance = inputArray[u];

            if (log) {
                distance = Math.log(distance);
            }
            if (inverse) {
                distance = 1/distance;
            }

            transformedDistances[u] = distance;
        }

        double meanDistance = DiscreteStatistics.mean(transformedDistances);
        if (!normalize) meanDistance = 1;

        for(int v=0; v<inputArray.length; v++) {
            transformedDistances[v] =  (transformedDistances[v]/meanDistance);
        }
        return transformedDistances;
    }

    private void outputArray(double[] array,  OutputFormat outputFormat, boolean nonreversible) {
        StringBuffer sb1 = new StringBuffer();
        StringBuffer sb2 = new StringBuffer();
        String sep;
        if (outputFormat == OutputFormat.TAB ) {
            sep = "\t";
        } else {
            sep = " ";
        }
        for(double anArray : array) {
            sb1.append(anArray + sep);
            sb2.append(1 + sep);
        }
        if (outputFormat == OutputFormat.XML ) {
            Element priorElement = new Element("multivariateGammaPrior");
            Element data = new Element("data");
            Element parameter1 = new Element("parameter");
            parameter1.setAttribute("idref","rates");
            data.addContent(parameter1);
            Element meanParameter = new Element("meanParameter");
            Element parameter2 = new Element("parameter");
            if (nonreversible) {
                parameter2.setAttribute("value",(sb1.toString()+sb1.toString()));
            } else {
                parameter2.setAttribute("value",sb1.toString());
            }
            meanParameter.addContent(parameter2);
            Element coefficientOfVariation = new Element("coefficientOfVariation");
            Element parameter3 = new Element("parameter");
            if (nonreversible) {
                parameter3.setAttribute("value",sb2.toString()+sb2.toString());
            } else {
                parameter3.setAttribute("value",sb2.toString());
            }
            coefficientOfVariation.addContent(parameter3);
            priorElement.addContent(data);
            priorElement.addContent(meanParameter);
            priorElement.addContent(coefficientOfVariation);
            XMLOutputter xmlOutputter = new XMLOutputter(Format.getPrettyFormat().setTextMode(Format.TextMode.PRESERVE));
            try {
                xmlOutputter.output(priorElement,resultsStream);
            } catch (IOException e) {
                System.err.println("IO Exception encountered: "+e.getMessage());
                System.exit(-1);
            }
            resultsStream.print("\r");
        }  else {
            if (nonreversible) {
                resultsStream.print(sb1.toString()+sb1.toString()+"\r");
            } else {
                resultsStream.print(sb1+"\r");
            }

        }

    }

    private void outputStringLine(String outputString,  OutputFormat outputFormat) {
        if (outputFormat == OutputFormat.XML ) {
           resultsStream.print("<!-- ");
        }
        resultsStream.print(outputString);
        if (outputFormat == OutputFormat.XML ) {
           resultsStream.print(" -->");
        } else {
           resultsStream.print("\r");
        }
    }

    // Messages to stderr, output to stdout
    private static final PrintStream progressStream = System.err;
    private PrintStream resultsStream;
    private static final String commandName = "discreteRatePriorGenerator";

    public static void printUsage(Arguments arguments) {

        arguments.printUsage(commandName, "[<output-file-name>]");
        progressStream.println();
        progressStream.println("  Example: " + commandName + " coordinates.txt ratePriors.txt");
        progressStream.println();
    }

    private static ArrayList parseCoordinatesFile(String inputFile, String[] locations, Double[] latitudes, Double[] longitudes) {
        ArrayList<Object[]> returnList = new ArrayList<Object[]>();
        List<String> countList = new ArrayList<String>();
        try{
            BufferedReader reader = new BufferedReader(new FileReader(inputFile));
            String line = reader.readLine();
            while (line != null && !line.equals("")) {
                countList.add(line);
                line = reader.readLine();
            }
        } catch (IOException e) {
            System.err.println("Error reading " + inputFile);
            System.exit(1);
        }

        if (countList.size()>0) {
            locations = new String[countList.size()];
            latitudes = new Double[countList.size()];
            longitudes = new Double[countList.size()];
            for(int i=0; i<countList.size(); i++) {
                StringTokenizer tokens = new StringTokenizer(countList.get(i));
                locations[i] = tokens.nextToken("\t");
                latitudes[i] = Double.parseDouble(tokens.nextToken("\t"));
                longitudes[i] = Double.parseDouble(tokens.nextToken("\t"));
                //System.out.println(locations[i]+"\t"+latitudes[i]+"\t"+longitudes[i]);
            }
        }
        returnList.add(locations);
        returnList.add(latitudes);
        returnList.add(longitudes);
        return returnList;

   }

    private static ArrayList parseSingleMeasureFile(String inputFile, String[] locations, Double[] densities) {

        ArrayList<Object[]> returnList = new ArrayList<Object[]>();

        boolean locationsSpecified = true;
        if (locations == null) {
            locationsSpecified = false;
        }

        List<String> countList = new ArrayList<String>();
        try{
            BufferedReader reader = new BufferedReader(new FileReader(inputFile));
            String line = reader.readLine();
            while (line != null && !line.equals("")) {
                countList.add(line);
                line = reader.readLine();
            }
        } catch (IOException e) {
            System.err.println("Error reading " + inputFile);
            System.exit(1);
        }

        if (countList.size()>0) {
            if (!locationsSpecified) {
                locations = new String[countList.size()];
            }
            densities = new Double[countList.size()];
            for(int i=0; i<countList.size(); i++) {
                StringTokenizer tokens = new StringTokenizer(countList.get(i));
                String location = tokens.nextToken("\t");
                if (!locationsSpecified) {
                    locations[i] = location;
                } else {
                    if (!(locations[i].equals(location))) {
                        System.err.println("Error in location specification in different files: " + locations[i] + "is not = " + location);
                    }
                }
                densities[i] = Double.parseDouble(tokens.nextToken("\t"));
            }
        }
        returnList.add(locations);
        returnList.add(densities);
        return returnList;
   }

    public static void main(String[] args) throws IOException {

        String outputFileName = null;
        String[] locations = null;
        Double[] latitudes = null;
        Double[] longitudes = null;
        Double[] densities = null;
        OutputFormat outputFormat = OutputFormat.XML;

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.StringOption(COORDINATES, "coordinate file", "specifies a tab-delimited file with coordinates for the locations"),
                        new Arguments.StringOption(DENSITIES, "density file", "specifies a tab-delimited file with densities for the locations"),
                        new Arguments.StringOption(GENERICS, "generics file", "specifies a tab-delimited file-list to use as measures for the locations"),
                        new Arguments.StringOption(FORMAT, TimeSlicer.enumNamesToStringArray(OutputFormat.values()),false,
                                "prior output format [default = XML]"),
                        //example: new Arguments.StringOption(TRAIT, "trait_name", "specifies an attribute-list to use to create a density map [default = location.rate]"),
                        //example: new Arguments.RealOption(MRSD,"specifies the most recent sampling data in fractional years to rescale time [default=0]"),
                        new Arguments.Option(HELP, "option to print this message"),
                });

        try {
            arguments.parseArguments(args);
        } catch (Arguments.ArgumentException ae) {
            progressStream.println(ae);
            printUsage(arguments);
            System.exit(1);
        }

        if (arguments.hasOption(HELP)) {
            printUsage(arguments);
            System.exit(0);
        }

        String coordinatesFileString = arguments.getStringOption(COORDINATES);
        if (coordinatesFileString != null) {
            ArrayList LocsLatsLongs = parseCoordinatesFile(coordinatesFileString,locations,latitudes,longitudes);
            locations = (String[]) LocsLatsLongs.get(0);
            latitudes = (Double[]) LocsLatsLongs.get(1);
            longitudes = (Double[]) LocsLatsLongs.get(2);
        }

        String densitiesFileString = arguments.getStringOption(DENSITIES);
        if (densitiesFileString != null) {
            ArrayList LocsDens = parseSingleMeasureFile(densitiesFileString,locations,densities);
            locations = (String[]) LocsDens.get(0);
            densities = (Double[]) LocsDens.get(1);
        }

        //TODO: support reading any measure (GENERICS)

        String summaryFormat = arguments.getStringOption(FORMAT);
        if (summaryFormat != null) {
            outputFormat = OutputFormat.valueOf(summaryFormat.toUpperCase());
        }

        final String[] args2 = arguments.getLeftoverArguments();

        switch (args2.length) {
            case 0:
                printUsage(arguments);
                System.exit(1);
            case 1:
                outputFileName = args2[0];
                break;
            default: {
                System.err.println("Unknown option: " + args2[1]);
                System.err.println();
                printUsage(arguments);
                System.exit(1);
            }
        }
        DiscreteRatePriorGenerator priors = new DiscreteRatePriorGenerator(locations, latitudes, longitudes, densities);
        priors.output(outputFileName, outputFormat);

        System.exit(0);

    }
    private static void printArray(Double[] array, String name) {
         try {
             PrintWriter outFile = new PrintWriter(new FileWriter(name), true);

             for(Double anArray : array) {
                 outFile.print(anArray + "\t");
             }
             outFile.close();
         } catch(IOException io) {
            System.err.print("Error writing to file: " + name);
         }
     }
    private static void printArray(String[] array, String name) {
         try {
             PrintWriter outFile = new PrintWriter(new FileWriter(name), true);

             for(String anArray : array) {
                 outFile.print(anArray + "\t");
             }
             outFile.close();
         } catch(IOException io) {
            System.err.print("Error writing to file: " + name);
         }
     }

}

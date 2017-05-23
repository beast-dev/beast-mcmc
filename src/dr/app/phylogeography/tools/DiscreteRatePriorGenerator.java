/*
 * DiscreteRatePriorGenerator.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.app.phylogeography.tools;


import dr.app.tools.TimeSlicer;
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
 * @author Philippe Lemey
 * @author Andrew Rambaut
 * @author Marc A. Suchard
 */
public class DiscreteRatePriorGenerator {

    public static final String HELP = "help";
    public static final String COORDINATES = "coordinates";
    public static final String DENSITIES = "densities";
    public static final String GENERICS = "generics";
    public static final String FORMAT = "format";
    public static final String MODEL = "model";

    public DiscreteRatePriorGenerator(String[] locations, Double[] latitudes, Double[] longitudes, Double[] densities) {

        if (locations == null) {
            //System.out.println(locations[0]);
            System.err.println("no locations specified!");
        } else {
            this.locations = locations;
        }

        this.latitudes = latitudes;
        this.longitudes = longitudes;
        this.densities = densities;

        if ((latitudes == null)||(longitudes == null)) {
            progressStream.println("no latitudes or longitudes specified!");
        } else {
            distances = getUpperTriangleDistanceMatrix(latitudes,longitudes);
        }

        if (densities != null) {
            densityDonorMatrix = getDensityMatrix(densities, true);
            densityRecipientMatrix = getDensityMatrix(densities, false);
        }

    }
    // for the time being, locations, latitudes and longitudes are not required
    private String[] locations;
    private Double[] latitudes;
    private Double[] longitudes;
    private final Double[] densities;
    private double[] distances;
    private double[] densityDonorMatrix;
    private double[] densityRecipientMatrix;

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

    private double[] getFullDistanceMatrix(Double[] latitudes, Double[] longitudes) {

        double[] distances = new double[latitudes.length * latitudes.length];

        int distanceCounter = 0;
        for (int a = 0; a < latitudes.length; a++) {
            //resultsStream.print(locations[a]+"\t");
            for (int b = 0; b < latitudes.length; b++) {
                distances[distanceCounter] = getKilometerGreatCircleDistance(latitudes[a],longitudes[a],latitudes[b],longitudes[b]);
                distanceCounter++;
            }
        }
        return distances;
    }

    private double[] getFullCoordDiffMatrix(Double[] coordinates, boolean positive, boolean negative) {

        double[] differences = new double[coordinates.length * coordinates.length];

        int differenceCounter = 0;
        for (int a = 0; a < coordinates.length; a++) {
            for (int b = 0; b < coordinates.length; b++) {
                double difference = coordinates[b]-coordinates[a];
                if (difference > 0 && positive) {
                    differences[differenceCounter] = difference;
                }  else if (difference < 0 && negative) {
                    differences[differenceCounter] = difference;
                } else {
                    differences[differenceCounter] = 0;
                }
                differenceCounter++;
            }
        }
        return differences;
    }

    private static double getKilometerGreatCircleDistance(double lat1, double long1, double lat2, double long2) {
         SphericalPolarCoordinates coord1 = new SphericalPolarCoordinates(lat1, long1);
         SphericalPolarCoordinates coord2 = new SphericalPolarCoordinates(lat2, long2);
         return (coord1.distance(coord2));
    }

    private double[] getDensityMatrix(Double[] densities, boolean donor) {

        double[] returnMatrix = new double[densities.length * (densities.length - 1)];
        int distanceCounter = 0;
        int matrixEntry;
        for (int c = 0; c < densities.length; c++) {
            for (int d = 0; d < densities.length; d++) {
                if (c == d) {
                    continue;
                }                
                if (donor) {
                    matrixEntry = c;
                } else {
                    matrixEntry = d;
                }
                returnMatrix[distanceCounter] = densities[matrixEntry];
                distanceCounter++;
            }
        }
        return returnMatrix;
    }

    private void printLocations(String[] locs, boolean upper) {

        int pairwiseCounter1 = 0;
        for (int c = 0; c < locs.length; c++) {
            pairwiseCounter1 ++;
            for (int d = pairwiseCounter1; d < locs.length; d++) {
                if (upper) {
                    System.out.println(locs[c]+"\t"+locs[d]);                    
                }  else {
                    System.out.println(locs[d]+"\t"+locs[c]);
                }
            }
        }
    }

    private void printLocationsDensitiesDistances (String[] locs, boolean upper) {

        double[] firstDensity = getSingleElementFromFullMatrix(densities,false);
        double[] secondDensity = getSingleElementFromFullMatrix(densities,true);

        System.out.println("location1\tlocation2\tdistance\tdensity1\tdensity2");
        int pairwiseCounter1 = 0;
        int arrayCounter = 0;
        for (int c = 0; c < locs.length; c++) {
            pairwiseCounter1 ++;
            for (int d = pairwiseCounter1; d < locs.length; d++) {
                if (upper) {
                    System.out.println(locs[c]+"\t"+locs[d]+"\t"+distances[arrayCounter]+"\t"+firstDensity[arrayCounter]+"\t"+secondDensity[arrayCounter]);
                    arrayCounter ++;
                }  else {
                    System.out.println(locs[d]+"\t"+locs[c]+"\t"+distances[arrayCounter]+"\t"+firstDensity[(firstDensity.length/2)+arrayCounter]+"\t"+secondDensity[(secondDensity.length/2)+arrayCounter]);
                    arrayCounter ++;
                    }
            }
        }
    }

    private void printFullDistanceMatrix(String name, boolean locationNames) {
        try {
            PrintWriter outFile = new PrintWriter(new FileWriter(name), true);

            outFile.print("location");
            for (int a = 0; a < locations.length; a++) {
                outFile.print(locations[a]+"\t");
                for (int b = 0; b < locations.length; b++) {
                    double lat1 = latitudes[a];
                    double lat2 = latitudes[b];
                    double long1 = longitudes[a];
                    double long2 = longitudes[b];
                    double distance = getKilometerGreatCircleDistance(lat1,long1,lat2,long2);
                    outFile.print(distance+"\t");
                }
                outFile.print("\r");
            }
            outFile.close();

        } catch(IOException io) {
           System.err.print("Error writing to file: " + name);
        }

    }

    enum OutputFormat {
        TAB,
        SPACE,
        XML
    }

    enum Model {
        PRIOR,
        GLM,
        JUMP
    }

    public void output(String outputFileName, OutputFormat outputFormat, Model model) {

        //printLocations(locations,false);
        //printLocationsDensitiesDistances(locations,true);
        //printLocationsDensitiesDistances(locations,false);
        //printFullDistanceMatrix("test.txt", false);
        resultsStream = System.out;

        if (outputFileName != null) {
            try {
            resultsStream = new PrintStream(new File(outputFileName));
            } catch (IOException e) {
                System.err.println("Error opening file: "+outputFileName);
                System.exit(1);
            }
        }

        if (model == model.PRIOR) {

            int predictor = 0;

            if (distances != null) {
                outputStringLine("reversible priors: distances",outputFormat);
                resultsStream.print("\r");
                outputArray(distances, outputFormat, model, predictor, false);
                resultsStream.print("\r");
                outputStringLine("reversible priors: normalized inverse distances",outputFormat);
                resultsStream.print("\r");
                outputArray(transform(distances,true, true, false, false), outputFormat, model, predictor, false);
                resultsStream.print("\r");
                outputStringLine("reversible priors: normalized inverse log distances",outputFormat);
                resultsStream.print("\r");
                outputArray(transform(distances,true, true, false, true), outputFormat, model, predictor, false);
                resultsStream.print("\r");
            }
            if (densities != null) {
                outputStringLine("reversible priors: densities",outputFormat);
                resultsStream.print("\r");
                outputArray(transform(getUpperTrianglePairwiseProductMatrix(densities),false,false,false,false), outputFormat, model, predictor, false);
                resultsStream.print("\r");
                outputStringLine("reversible priors: normalized densities",outputFormat);
                resultsStream.print("\r");
                outputArray(transform(getUpperTrianglePairwiseProductMatrix(densities),false,true,false,false), outputFormat, model, predictor, false);
                resultsStream.print("\r");
                outputStringLine("reversible priors: normalized log densities",outputFormat);
                resultsStream.print("\r");
                outputArray(transform(getUpperTrianglePairwiseProductMatrix(densities),false,true,false,true), outputFormat, model, predictor, false);
                resultsStream.print("\r");

                if (distances != null) {
                    outputStringLine("reversible priors: product of normalized densities divided by normalized distances",outputFormat);
                    resultsStream.print("\r");
                    outputArray(productOfArrays(transform(getUpperTrianglePairwiseProductMatrix(densities),false,true,false,false),transform(distances,false, true,false,false),true), outputFormat, model, predictor, false);
                    resultsStream.print("\r");
                    outputStringLine("reversible priors: normalized product of densities divided by distances",outputFormat);
                    resultsStream.print("\r");
                    outputArray(transform(productOfArrays(getUpperTrianglePairwiseProductMatrix(densities),distances, true),false,true,false,false), outputFormat, model, predictor, false);
                    resultsStream.print("\r");
                    outputStringLine("reversible priors: normalized log product of densities divided by distances",outputFormat);
                    resultsStream.print("\r");
                    outputArray(transform(productOfArrays(getUpperTrianglePairwiseProductMatrix(densities),distances, true),false,true,false,true), outputFormat, model, predictor, false);
                    resultsStream.print("\r");
                }

            }
            if (distances != null) {
                outputStringLine("nonreversible priors: distances",outputFormat);
                resultsStream.print("\r");
                outputArray(distances, outputFormat, model, predictor, true);
                resultsStream.print("\r");
                outputStringLine("nonreversible priors: normalized inverse distances",outputFormat);
                resultsStream.print("\r");
                outputArray(transform(distances,true, true, false, false), outputFormat, model, predictor, true);
                resultsStream.print("\r");
                outputStringLine("nonreversible priors: normalized inverse log distances",outputFormat);
                resultsStream.print("\r");
                outputArray(transform(distances,true, true, false, true), outputFormat, model,predictor, true);
                resultsStream.print("\r");
            }
            if (densities != null) {
                outputStringLine("nonreversible priors: densities",outputFormat);
                resultsStream.print("\r");
                outputArray(transform(getUpperTrianglePairwiseProductMatrix(densities),false,false,false,false), outputFormat, model, predictor, true);
                resultsStream.print("\r");
                outputStringLine("nonreversible priors: normalized densities",outputFormat);
                resultsStream.print("\r");
                outputArray(transform(getUpperTrianglePairwiseProductMatrix(densities),false,true,false,false), outputFormat, model, predictor, true);
                resultsStream.print("\r");
                outputStringLine("nonreversible priors: normalized log densities",outputFormat);
                resultsStream.print("\r");
                outputArray(transform(getUpperTrianglePairwiseProductMatrix(densities),false,true,false,true), outputFormat, model, predictor, true);
                resultsStream.print("\r");

                if (distances != null) {
                    outputStringLine("nonreversible priors: product of normalized densities divided by normalized distances",outputFormat);
                    resultsStream.print("\r");
                    outputArray(productOfArrays(transform(getUpperTrianglePairwiseProductMatrix(densities),false,true,false,false),transform(distances,false, true, false, false),true), outputFormat, model, predictor, true);
                    resultsStream.print("\r");
                    outputStringLine("nonreversible priors: normalized product of densities divided by distances",outputFormat);
                    resultsStream.print("\r");
                    outputArray(transform(productOfArrays(getUpperTrianglePairwiseProductMatrix(densities),distances, true),false,true,false,false), outputFormat, model, predictor, true);
                    resultsStream.print("\r");
                    outputStringLine("nonreversible priors: normalized log product of densities divided by distances",outputFormat);
                    resultsStream.print("\r");
                    outputArray(transform(productOfArrays(getUpperTrianglePairwiseProductMatrix(densities),distances, true),false,true,false,true), outputFormat, model, predictor, true);
                    resultsStream.print("\r");
                }

            }
        } else if (model == model.GLM) {

            int predictor = 1;

            //TODO: fully implement glm output
            if (distances != null) {
                //printLocations(locations, true);
                //printLocations(locations, false);
                outputStringLine("predictor: distances",outputFormat);
                resultsStream.print("\r");
                outputArray(transform(distances,false, false, true, false), outputFormat, model, predictor, true);
                predictor++;
                outputStringLine("predictor: standardized log distances",outputFormat);
                resultsStream.print("\r");
                outputArray(transform(distances,false, false, true, true), outputFormat, model, predictor, true);
            }
            if (densities != null) {
                predictor++;
                outputStringLine("predictor: standardized donor densities",outputFormat);
                resultsStream.print("\r");
                outputArray(transform(getSingleElementFromFullMatrix(densities, false),false,false,true,true), outputFormat, model, predictor, false);
                predictor++;
                outputStringLine("predictor: standardized recipient densities",outputFormat);
                resultsStream.print("\r");
                outputArray(transform(getSingleElementFromFullMatrix(densities, true),false,false,true,true), outputFormat, model, predictor, false);

  // products not necessary, and shouldn't be normalized but standardized, also log
  //              outputStringLine("predictor: normalized product of densities divided by distances",outputFormat);
  //              resultsStream.print("\r");
  //              outputArray(transform(productOfArrays(getUpperTrianglePairwiseProductMatrix(densities),distances, true),false,true,false,false), outputFormat, predictor, true);
  //              resultsStream.print("\r");
  //              predictor++;
  //              outputStringLine("predictor: normalized log product of densities divided by distances",outputFormat);
  //              resultsStream.print("\r");
  //              outputArray(transform(productOfArrays(getUpperTrianglePairwiseProductMatrix(densities),distances, true),false,true,false,true), outputFormat, predictor, true);
  //              resultsStream.print("\r");
            }


        } else if (model == model.JUMP) {

            int predictor = 0;

            outputStringLine("great circle distance jump matrix",outputFormat);
            resultsStream.print("\r");
            outputArray(transform(getFullDistanceMatrix(latitudes, longitudes),false,false,false,false), outputFormat, model, predictor, false);
            resultsStream.print("\r");
            predictor ++;

            outputStringLine("latitude jump matrix",outputFormat);
            resultsStream.print("\r");
            outputArray(transform(getFullCoordDiffMatrix(latitudes, true, true),false,false,false,false), outputFormat, model, predictor, false);
            resultsStream.print("\r");
            predictor ++;

            outputStringLine("longitude jump matrix",outputFormat);
            resultsStream.print("\r");
            outputArray(transform(getFullCoordDiffMatrix(longitudes, true, true),false,false,false,false), outputFormat, model, predictor, false);
            resultsStream.print("\r");
            predictor ++;

            outputStringLine("westward jump matrix",outputFormat);
            resultsStream.print("\r");
            outputArray(transform(getFullCoordDiffMatrix(longitudes, false, true),false,false,false,false), outputFormat, model, predictor, false);
            resultsStream.print("\r");
            predictor ++;

            outputStringLine("eastward jump matrix",outputFormat);
            resultsStream.print("\r");
            outputArray(transform(getFullCoordDiffMatrix(longitudes, true, false),false,false,false,false), outputFormat, model, predictor, false);
            resultsStream.print("\r");
            predictor ++;

            outputStringLine("northward jump matrix",outputFormat);
            resultsStream.print("\r");
            outputArray(transform(getFullCoordDiffMatrix(latitudes, true, false),false,false,false,false), outputFormat, model, predictor, false);
            resultsStream.print("\r");
            predictor ++;

            outputStringLine("southward jump matrix",outputFormat);
            resultsStream.print("\r");
            outputArray(transform(getFullCoordDiffMatrix(latitudes, false, true),false,false,false,false), outputFormat, model, predictor, false);
            resultsStream.print("\r");
            predictor ++;

//            outputStringLine("latitude reward parameter",outputFormat);
//            resultsStream.print("\r");
//            outputArray(transform(objectToPrimitiveArray(latitudes),false,false,false,false), outputFormat, model, predictor, false);
//            resultsStream.print("\r");
//            predictor ++;

//            outputStringLine("longitude reward parameter",outputFormat);
//            resultsStream.print("\r");
//            outputArray(transform(objectToPrimitiveArray(longitudes),false,false,false,false), outputFormat, model, predictor, false);
//            resultsStream.print("\r");
//            predictor ++;

            for (int i = 0; i < locations.length; i++) {

                double[] locationReward = getLocationReward(i,locations.length);
                outputStringLine(locations[i]+" reward parameter",outputFormat);
                resultsStream.print("\r");
                outputArray(transform(locationReward,false,false,false,false), outputFormat, model, predictor, false);
                resultsStream.print("\r");
                predictor ++;

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

    private double[] getSingleElementFromFullMatrix(Double[] matrixValues, boolean secondElement) {

        double[] singleElements = new double[(matrixValues.length * (matrixValues.length - 1))];
        //System.out.println("matrixsize "+(matrixValues.length * (matrixValues.length - 1)));

        //get upper matrix
        int counter = 0;
        int pairwiseCounter1 = 0;
        for(Double matrixValue : matrixValues) {
            pairwiseCounter1++;
            for(int d = pairwiseCounter1; d < matrixValues.length; d++) {
                if (secondElement) {
                    singleElements[counter] = matrixValues[d];
                } else {
                    singleElements[counter] = matrixValue;

                }
                counter++;
            }
        }
        //get lower matrix
        int pairwiseCounter2 = 0;
        for(Double matrixValue : matrixValues) {
            pairwiseCounter2++;
            for(int d = pairwiseCounter2; d < matrixValues.length; d++) {
                if (secondElement) {
                    singleElements[counter] = matrixValue;
                } else {
                    singleElements[counter] = matrixValues[d];
                }
                counter++;
            }
        }

        //System.out.println("counter "+counter);

        return singleElements;
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

    // normalize is really rescaling the vector to have a mean = 1 here, standardize is rescaling it to have mean = 0 and variance =1
    private double[] transform(double[] inputArray, boolean inverse, boolean normalize, boolean standardize, boolean log) {

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

        double meanDistance = 1;
        double stdev = 0;

        if (normalize || standardize) {
            meanDistance = DiscreteStatistics.mean(transformedDistances);
        }
        if (standardize) {
            stdev = Math.sqrt(DiscreteStatistics.variance(transformedDistances));
        }

        for(int v=0; v<inputArray.length; v++) {
            if (normalize) {
                transformedDistances[v] =  (transformedDistances[v]/meanDistance);
            } else if (standardize) {
                transformedDistances[v] =  ((transformedDistances[v] - meanDistance)/stdev);
            }
        }
        return transformedDistances;
    }

    private void outputArray(double[] array,  OutputFormat outputFormat, Model model, int predictor, boolean nonreversible) {

        StringBuilder sb1 = new StringBuilder();
        StringBuilder sb2 = new StringBuilder();
        String sep;
        if (outputFormat == OutputFormat.TAB ) {
            sep = "\t";
        } else {
            sep = " ";
        }
        //String newLine = "\n";
        //here we break up the jump parameter output in a lenght*by*length matrix
//        int entryCounter = 1;
//        double length = Math.sqrt(array.length);
//        System.out.println(length);
//        for(double anArray : array) {
//            if (model == model.JUMP) {
//                if ( (entryCounter % length) > 0 ) {
//                    sb1.append(anArray + sep);
//                    sb2.append(1 + sep);
//                    entryCounter ++;
//               } else {
//                    System.out.println("return");
//                    sb1.append(newLine+anArray + sep);
//                    sb2.append(newLine+1 + sep);
//                     entryCounter ++;
//               }
//
//            } else {
//                sb1.append(anArray + sep);
//                sb2.append(1 + sep);
//            }
//       }
         for(double anArray : array) {
                sb1.append(anArray + sep);
                sb2.append(1 + sep);
         }
        
        if (outputFormat == OutputFormat.XML ) {
            if (model == model.GLM) {
                Element parameter = new Element("parameter");
                parameter.setAttribute("id","predictor"+predictor);
                if (nonreversible) {
                    parameter.setAttribute("value",(sb1.toString()+sb1.toString()));
                } else {
                    parameter.setAttribute("value",sb1.toString());
                }
                XMLOutputter xmlOutputter = new XMLOutputter(Format.getPrettyFormat().setTextMode(Format.TextMode.PRESERVE));
                try {
                    xmlOutputter.output(parameter,resultsStream);
                } catch (IOException e) {
                    System.err.println("IO Exception encountered: "+e.getMessage());
                    System.exit(-1);
                }
                resultsStream.print("\r");
            }  else if (model == model.PRIOR) {
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
            }  else if (model == model.JUMP) {

                Element parameter = new Element("parameter");
                parameter.setAttribute("id","jump"+predictor);
                parameter.setAttribute("value",sb1.toString());
                XMLOutputter xmlOutputter = new XMLOutputter(Format.getPrettyFormat().setTextMode(Format.TextMode.PRESERVE));
                try {
                    xmlOutputter.output(parameter,resultsStream);
                } catch (IOException e) {
                    System.err.println("IO Exception encountered: "+e.getMessage());
                    System.exit(-1);
                }
                resultsStream.print("\r");
            }
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

    private double[] objectToPrimitiveArray(Double[] array) {
        double[] returnArray = new double[array.length];
        int counter = 0;
        for(Double anArray : array) {
           returnArray[counter] = anArray.doubleValue();
           counter ++;
        }
        return returnArray;
    }

    private double[] getLocationReward(int indicator, int length) {
        double[] returnArray = new double[length];
        for (int i = 0; i < length; i++) {
            if (i == indicator) {
                returnArray[i] = 1.0;
            } else {
                returnArray[i] = 0.0;
            }
        }
        return returnArray;
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
        Model model = Model.PRIOR;

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.StringOption(COORDINATES, "coordinate file", "specifies a tab-delimited file with coordinates for the locations"),
                        new Arguments.StringOption(DENSITIES, "density file", "specifies a tab-delimited file with densities for the locations"),
                        new Arguments.StringOption(GENERICS, "generics file", "specifies a tab-delimited file-list to use as measures for the locations"),
                        new Arguments.StringOption(FORMAT, TimeSlicer.enumNamesToStringArray(OutputFormat.values()),false,
                                "prior output format [default = XML]"),
                        new Arguments.StringOption(MODEL, TimeSlicer.enumNamesToStringArray(Model.values()),false,
                                "model output [default = rate priors]"),
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

        String modelComponent = arguments.getStringOption(MODEL);
        if (modelComponent != null) {
            model = Model.valueOf(modelComponent.toUpperCase());
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
        DiscreteRatePriorGenerator rates = new DiscreteRatePriorGenerator(locations, latitudes, longitudes, densities);
        rates.output(outputFileName, outputFormat, model);

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

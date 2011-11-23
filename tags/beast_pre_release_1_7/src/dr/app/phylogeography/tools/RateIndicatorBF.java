/*
 * RateIndicatorBF.java
 *
 * Copyright (C) 2002-2010 BEAST Development Team
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
import dr.util.HeapSort;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;

import org.jdom.Element;
import org.jdom.output.XMLOutputter;
import org.jdom.output.Format;

/**
 * @author Philippe Lemey
 * @author Andrew Rambaut
 * @author Marc A. Suchard
 */
public class RateIndicatorBF {

    public static final String LOCATIONSFILE = "locationsfile";
    public static final String HELP = "help";
    public static final String BURNIN = "burnin";
    public static final String KML = "kml";
    public static final String PMEAN = "pmean";
    public static final String POFFSET = "poffset";
    public static final String ISTRING = "istring";
    public static final String RSTRING = "rstring";
    public static final String PSTRING = "pstring";
    public static final String GSTRING = "gstring";
    public static final String LOWCOLOR = "lowcolor";
    public static final String UPCOLOR = "upcolor";
    public static final String WIDTH = "width";
    public static final String KMLFILE = "kmlfile";
    public static final String BFCUTOFF = "bfcutoff";
    public static final String ICUTOFF = "icutoff";
    public static final String LOCATIONSTATES = "locationstates";
    public static final String BWC = "bwc";
    public static final String BWM = "bwm";
    public static final String ALTITUDE = "altitude";

    public static final String[] falseTrue = new String[] {"false","true"};

    private static PrintStream progressStream = System.out;
    private static final String commandName = "rateIndicatorBF";
    public static void printUsage(Arguments arguments) {

        arguments.printUsage(commandName, "<input-file-name> [<output-file-name>]");
        progressStream.println();
        progressStream.println("  Example: " + commandName + " indicator.log rates.out");
        progressStream.println();
    }

    public RateIndicatorBF (String inputFileName, int burnin, String rateIndicatorString, int numberOfStates, String[][] locations, boolean bayesFactor, double cutoff, double meanPoissonPrior, int offsetPoissonPrior, String actualRateString, String relativeRateString, String geoSiteModelString) {

        //count the number of states in the RateIndicatorLog file
        generationCount = getGenerationCount(inputFileName);

        //find the first rateIndicator in the RateIndicatorLog file
        int firstRateIndicator = getFirstEntryOf(inputFileName, rateIndicatorString);
        //progressStream.println("first rateIndicator is at column "+firstRateIndicator);

        //count the rateIndicators in the RateIndicatorLog file
        int numberOfRateIndicators = getNumberOfEntries(inputFileName, firstRateIndicator, rateIndicatorString);

        if (numberOfStates > 0) {
            statesCounter = numberOfStates;
            progressStream.println("number of states provided = "+statesCounter);
        } else {
            statesCounter = locations.length;
            progressStream.println("number of states in coordinates file = "+statesCounter);
        }

        this.bayesFactor = bayesFactor;
        this.cutoff = cutoff;
        this.meanPoissonPrior = meanPoissonPrior;
        this.offsetPoissonPrior = offsetPoissonPrior;
        this.actualRateString = actualRateString;
        this.relativeRateString = relativeRateString;
        this.geoSiteModelString = geoSiteModelString;


        if (numberOfRateIndicators != ((statesCounter*(statesCounter-1.0))/2.0)) {
            if (numberOfRateIndicators == (statesCounter*(statesCounter-1.0))) {
                progressStream.println("K*(K-1) rateIndicators, with K = "+ statesCounter+". So, nonreversible matrix!");
                nonreversible = true;
            } else {
                System.err.println("the number of rateIndicators ("+numberOfRateIndicators+") does not match (K*(K-1)/2) states ("+((statesCounter*(statesCounter-1.0))/2.0)+"; with K = "+statesCounter+").");
            }
        }

        // so now we know the dimension of the rateIndicator array
        if ((generationCount - burnin)< 10) {
            System.err.println("With burn-in = "+burnin+", there are only "+(generationCount - burnin)+" state(s) in indicator log file??");
        }
        double[][] rateIndicators = new double[((generationCount - 1)-burnin)][numberOfRateIndicators];
        fillRateIndicatorArray(inputFileName,rateIndicators,burnin,firstRateIndicator,numberOfRateIndicators);
        //print2DArray(rateIndicators, "test.txt");
        expectedRateIndicators = meanCol(rateIndicators);

        //compile locations
        locationNames = new String [numberOfRateIndicators][2];
        longitudes = new double [numberOfRateIndicators][2];
        latitudes = new double [numberOfRateIndicators][2];

        this.locations = locations;
        compileLocations(locations,locationNames,latitudes,longitudes);
        supportedRateIndicators = getSupportedRateIndicators();

        //below is for rateSummary

        int numberOfActualRates = numberOfRateIndicators;

        //boolean to see it the rateLog contains productStatitistics, if not, we make 'em ourselves
        boolean containsActualRates = hasEntryOf(inputFileName, actualRateString);
        if (!containsActualRates){
            // if there are no productStatistics, we will look for the relativ rates instead
            actualRateString = relativeRateString;
        }

        int firstActualRate = getFirstEntryOf(inputFileName, actualRateString);

//        int geoSiteModelMuPosition = getFirstEntryOf(inputFileName, geoSiteModelString);
        //System.out.println("geoSiteModelMu is at column "+geoSiteModelMuPosition);



    }

    private double[] getSupportedRateIndicators(){

        int indicatorCounter = 0;
        for (int p = 0; p < expectedRateIndicators.length; p++){

            double indicator = expectedRateIndicators[p];
            if (!bayesFactor) {
                if (indicator > cutoff) {
                    indicatorCounter ++;
                }
            }  else {
                if (indicator > getBayesFactorCutOff(cutoff, meanPoissonPrior, offsetPoissonPrior, statesCounter, nonreversible)) {
                    indicatorCounter ++;
                }
            }
        }

        double[] supportedIndicators = new double[indicatorCounter];
        supportedBFs = new double[indicatorCounter];
        supportedLocations = new String[indicatorCounter][2];
        supportedLatitudes = new double[indicatorCounter][2];
        supportedLongitudes = new double[indicatorCounter][2];
        int[] indices = new int[expectedRateIndicators.length];
        HeapSort.sort(expectedRateIndicators, indices);
        int fillCount = 0;

        for (int o = 0; o < expectedRateIndicators.length; o++){
            //we order rate indicators in decreasing order
            double indicator = expectedRateIndicators[indices[expectedRateIndicators.length - o - 1]];
            double BF = getBayesFactor(indicator, meanPoissonPrior, statesCounter, offsetPoissonPrior, nonreversible, 0);
            if (BF == Double.POSITIVE_INFINITY) {
                BF = getBayesFactor(indicator, meanPoissonPrior, statesCounter, offsetPoissonPrior, nonreversible, generationCount);
            }

            double threshold = (bayesFactor ? getBayesFactorCutOff(cutoff, meanPoissonPrior, offsetPoissonPrior, statesCounter, nonreversible) : cutoff);
            if (indicator > threshold) {
                supportedIndicators[fillCount] = indicator;
                supportedBFs[fillCount] = BF;
                supportedLocations[fillCount][0] =  locationNames[indices[(expectedRateIndicators.length - o - 1)]][0];
                supportedLocations[fillCount][1] =  locationNames[indices[(expectedRateIndicators.length - o - 1)]][1];
                supportedLatitudes[fillCount][0] =  latitudes[indices[(expectedRateIndicators.length - o - 1)]][0];
                supportedLatitudes[fillCount][1] =  latitudes[indices[(expectedRateIndicators.length - o - 1)]][1];
                supportedLongitudes[fillCount][0] =  longitudes[indices[(expectedRateIndicators.length - o - 1)]][0];
                supportedLongitudes[fillCount][1] =  longitudes[indices[(expectedRateIndicators.length - o - 1)]][1];
                fillCount ++;
            }
        }
        return supportedIndicators;
    }

    private void compileLocations(String[][] locations, String[][] locationNames, double[][] latitudes, double[][] longitudes){
        //begin of new code
        int elementCounter = 0;
        int secondCounter = 0;
        for (int i = 0; i < (statesCounter - 1); i++) {

            secondCounter ++;

            for (int j = secondCounter; j < statesCounter; j++) {

                if (locations != null) {
                    locationNames[elementCounter][0] = locations[i][0];
                    longitudes[elementCounter][0] = Double.parseDouble(locations[i][2]);
                    latitudes[elementCounter][0] = Double.parseDouble(locations[i][1]);
                    locationNames[elementCounter][1] = locations[j][0];
                    longitudes[elementCounter][1] = Double.parseDouble(locations[j][2]);
                    latitudes[elementCounter][1] = Double.parseDouble(locations[j][1]);
                } else {
                    locationNames[elementCounter][0] = "location"+(i+1);
                    longitudes[elementCounter][0] = Double.NaN;
                    latitudes[elementCounter][0] = Double.NaN;
                    locationNames[elementCounter][1] = "location"+(j+1);
                    longitudes[elementCounter][1] = Double.NaN;
                    latitudes[elementCounter][1] = Double.NaN;
                }

                elementCounter ++;
            }
        }
        // for nonreversible models, we keep on filling the arrays
        if (nonreversible) {
            for (int k = 0; k < elementCounter; k++) {

                locationNames[elementCounter + k][0] = locationNames[k][1];
                longitudes[elementCounter + k][0] = longitudes[k][1];
                latitudes[elementCounter + k][0] = latitudes[k][1];
                locationNames[elementCounter + k][1] = locationNames[k][0];
                longitudes[elementCounter + k][1] = longitudes[k][0];
                latitudes[elementCounter + k][1] = latitudes[k][0];
            }
        }

    }

    public void outputKML(String KMLoutputFile,String lowerLinkColor,String upperLinkColor, double branchWidthConstant, double branchWidthMultiplier, double altitudeFactor){

        double divider = 100;

        PrintStream resultsStream  = System.out;
        if (KMLoutputFile != null) {
            try {
                resultsStream = new PrintStream(new File(KMLoutputFile));
            } catch (IOException e) {
                System.err.println("Error opening file: "+KMLoutputFile);
                System.exit(-1);
            }
        }

        Element rootElement = new Element("kml");
        Element documentElement =  new Element("Document");
        Element folderElement = new Element("Folder");

        Element documentNameElement = new Element("name");
        documentNameElement.addContent(KMLoutputFile);
        documentElement.addContent(documentNameElement);

        List<Element> schema = new ArrayList<Element>();
        Element locationSchema = new Element("Schema");
        locationSchema.setAttribute("id", "Locations_Schema");
        locationSchema.addContent(new Element("SimpleField")
                .setAttribute("name", "Name")
                .setAttribute("type", "string"));
        locationSchema.addContent(new Element("SimpleField")
                .setAttribute("name", "In")
                .setAttribute("type", "number"));
        locationSchema.addContent(new Element("SimpleField")
                .setAttribute("name", "Out")
                .setAttribute("type", "number"));
        locationSchema.addContent(new Element("SimpleField")
                .setAttribute("name", "Total")
                .setAttribute("type", "number"));
        Element ratesSchema = new Element("Schema");
        ratesSchema.setAttribute("id", "Rates_Schema");
        ratesSchema.addContent(new Element("SimpleField")
                .setAttribute("name", "Name")
                .setAttribute("type", "string"));
        ratesSchema.addContent(new Element("SimpleField")
                .setAttribute("name", "BF")
                .setAttribute("type", "double"));
        ratesSchema.addContent(new Element("SimpleField")
                .setAttribute("name", "Indicator")
                .setAttribute("type", "double"));
        schema.add(ratesSchema);

        documentElement.addContent(schema);

        Element folderNameElement = new Element("name");
        String cutoffString;
        if (bayesFactor) {
            cutoffString = "bayes factor";
        } else {
            cutoffString = "indicator";
        }
        folderNameElement.addContent("discrete rates with "+cutoffString+" larger than "+cutoff);
        folderElement.addContent(folderNameElement);

        double[] minMax = new double[2];
        if (supportedRateIndicators.length > 0) {
            minMax[0] = supportedRateIndicators[supportedRateIndicators.length - 1];
            minMax[1] = supportedRateIndicators[0];
        } else {
            minMax[0] = minMax[1] = 0;
            System.err.println("No rate indicators above the specified cut-off!");
        }

//        for (int p = 0; p < supportedRateIndicators.length; p++){
//            addRateAndStyle(supportedRateIndicators[p], supportedLongitudes[p][0], supportedLatitudes[p][0], supportedLongitudes[p][1], supportedLatitudes[p][1], p+1, branchWidthConstant, branchWidthMultiplier, minMax, altitudeFactor, divider, lowerLinkColor, upperLinkColor, folderElement, documentElement);
//        }

        for (int p = 0; p < supportedRateIndicators.length; p++){
            addRateWithData(supportedRateIndicators[p], supportedBFs[p], supportedLocations[p][0], supportedLocations[p][1], supportedLongitudes[p][0], supportedLatitudes[p][0], supportedLongitudes[p][1], supportedLatitudes[p][1], p+1, folderElement);
        }

        //add locations
        Element folder2Element = new Element("Folder");
        Element folderName2Element = new Element("name");
        folderName2Element.addContent("Locations");
        folder2Element.addContent(folderName2Element);
        addLocations(folder2Element);
        documentElement.addContent(folder2Element);

        documentElement.addContent(folderElement);
        rootElement.addContent(documentElement);

        XMLOutputter xmlOutputter = new XMLOutputter(Format.getPrettyFormat().setTextMode(Format.TextMode.PRESERVE));
        try {
            xmlOutputter.output(rootElement,resultsStream);
        } catch (IOException e) {
            System.err.println("IO Exception encountered: "+e.getMessage());
            System.exit(-1);
        }

    }

    private void addLocations(Element folderElement){


        for (int i = 0; i < locations.length ; i++) {
            Element placemarkElement  = new Element("Placemark");
            Element placemarkNameElement = new Element("name");
            placemarkNameElement.addContent(locations[i][0]);
            placemarkElement.addContent(placemarkNameElement);

            int inCount = 0;
            int outCount = 0;

            for (int j = 0; j < supportedRateIndicators.length; j++) {
                if (supportedLocations[j][0].equals(locations[i][0])) {
                    inCount ++;
                }
                if (supportedLocations[j][1].equals(locations[i][0])) {
                    outCount ++;
                }
            }
            int totalCount = inCount + outCount;

            Element data = new Element("ExtendedData");
            Element schemaData = new Element("SchemaData");
            schemaData.setAttribute("schemaUrl", "#Location_Schema");
            schemaData.addContent(new Element("SimpleData").setAttribute("name", "Name").addContent(locations[i][0]));
            schemaData.addContent(new Element("SimpleData").setAttribute("name", "In").addContent(Integer.toString(inCount)));
            schemaData.addContent(new Element("SimpleData").setAttribute("name", "Out").addContent(Integer.toString(outCount)));
            schemaData.addContent(new Element("SimpleData").setAttribute("name", "Total").addContent(Integer.toString(totalCount)));
            data.addContent(schemaData);
            placemarkElement.addContent(data);

            Element pointElement  = new Element("Point");
            Element altitude  = new Element("altitudeMode");
            altitude.addContent("relativeToGround");
            pointElement.addContent(altitude);
            Element coordinates  = new Element("coordinates");
            coordinates.addContent(locations[i][2]+","+locations[i][1]+",0");
            pointElement.addContent(coordinates);
            placemarkElement.addContent(pointElement);

            folderElement.addContent(placemarkElement);
        }
    }

    private void  addRateAndStyle(double rateIndicator, double startLongitude, double startLatitude, double endLongitude, double endLatitude, int number, double branchWidthConstant, double branchWidthMultiplier, double[] minAndMaxRateIndicator, double altitudeFactor, double divider, String lowerLinkColor, String upperLinkColor, Element folderElement, Element documentElement) {

        String opacity = "FF";

        double distance = (3958*Math.PI*Math.sqrt((endLatitude-startLatitude)*(endLatitude-startLatitude)+Math.cos(endLatitude/57.29578)*Math.cos(startLatitude/57.29578)*(endLongitude-startLongitude)*(endLongitude-startLongitude))/180);

        double maxAltitude = distance*altitudeFactor;
        double latitudeDifference = endLatitude - startLatitude;
        double longitudeDifference = endLongitude - startLongitude;


        boolean longitudeBreak = false; //if we go through the 180
        if (endLongitude*startLongitude < 0) {

            double trialDistance = 0;

            if (endLongitude < 0) {
                trialDistance += (endLongitude + 180);
                trialDistance += (180 - startLongitude);
                //System.out.println(parentLongitude+"\t"+longitude+"\t"+trialDistance+"\t"+longitudeDifference);
            } else {
                trialDistance += (startLongitude + 180);
                trialDistance += (180 - endLongitude);
                //System.out.println(parentLongitude+"\t"+longitude+"\t"+trialDistance+"\t"+longitudeDifference);
            }

            if (trialDistance < Math.abs(longitudeDifference)) {
                longitudeDifference = trialDistance;
                longitudeBreak = true;
                //System.out.println("BREAK!"+longitudeDifference);
            }
        }

        Element styleElement = new Element("Style");
        styleElement.setAttribute("id","rate"+number+"_style");
        Element lineStyle = new Element("LineStyle");
        Element width = new Element("width");
        width.addContent(Double.toString(branchWidthConstant+branchWidthMultiplier*rateIndicator));
        Element color = new Element("color");
        String colorString = TimeSlicer.getKMLColor(rateIndicator,minAndMaxRateIndicator,lowerLinkColor,upperLinkColor);
        color.addContent(opacity+colorString);
        lineStyle.addContent(width);
        lineStyle.addContent(color);
        styleElement.addContent(lineStyle);
        documentElement.addContent(styleElement);

        double currentLongitude1 = 0;  //we need this if we go through the 180
        double currentLongitude2 = 0;  //we need this if we go through the 180

        for (int a = 0; a < divider; a ++) {

            Element placemarkElement  = new Element("Placemark");

            Element placemarkNameElement = new Element("name");
            String name = "rate"+number+"_part"+(a+1);
            placemarkNameElement.addContent(name);
            placemarkElement.addContent(placemarkNameElement);

            Element placemarkStyleElement = new Element("styleUrl");
            placemarkStyleElement.addContent("#rate"+number+"_style");
            placemarkElement.addContent(placemarkStyleElement);

            Element lineStringElement = new Element("LineString");
            Element altitudeMode = new Element("altitudeMode");
            altitudeMode.addContent("relativeToGround");
            lineStringElement.addContent(altitudeMode);
            Element tessellate = new Element("tessellate");
            tessellate.addContent("1");
            lineStringElement.addContent(tessellate);
            Element coordinatesElement = new Element("coordinates");
            StringBuffer coordinatesStartBuffer = new StringBuffer();
            StringBuffer coordinatesEndBuffer = new StringBuffer();

            if (longitudeBreak) {

                if (startLongitude > 0) {
                    currentLongitude1 = startLongitude+a*(longitudeDifference/divider);

                    if (currentLongitude1 < 180) {
                        coordinatesStartBuffer.append(currentLongitude1+",");
                        //System.out.println("1 currentLongitude1 < 180\t"+currentLongitude1+"\t"+longitude);

                    } else {
                        coordinatesStartBuffer.append((-180-(180-currentLongitude1))+",");
                        //System.out.println("2 currentLongitude1 > 180\t"+currentLongitude1+"\t"+(-180-(180-currentLongitude1))+"\t"+longitude);
                    }
                } else {
                    currentLongitude1 = startLongitude-a*(longitudeDifference/divider);

                    if (currentLongitude1 > (-180)) {
                        coordinatesStartBuffer.append(currentLongitude1+",");
                        //System.out.println("currentLongitude1 > -180\t"+currentLongitude1+"\t"+longitude);
                    } else {
                        coordinatesStartBuffer.append((180+(currentLongitude1+180))+",");
                        //System.out.println("currentLongitude1 > -180\t"+(180+(currentLongitude1+180))+"\t"+longitude);
                    }
                }

            } else {
                coordinatesStartBuffer.append((startLongitude+a*(longitudeDifference/divider))+",");
            }

            coordinatesStartBuffer.append((startLatitude+a*(latitudeDifference/divider))+","+(maxAltitude*Math.sin(Math.acos(1 - a*(1.0/(divider/2.0))))));

            if (longitudeBreak) {

                if (startLongitude > 0) {
                    currentLongitude2 = startLongitude+(a+1)*(longitudeDifference/divider);

                    if (currentLongitude2 < 180) {
                        coordinatesEndBuffer.append((currentLongitude2)+",");
                    } else {
                        coordinatesEndBuffer.append((-180-(180-currentLongitude2))+",");
                    }
                } else {
                    currentLongitude2 = startLongitude-(a+1)*(longitudeDifference/divider);

                    if (currentLongitude2 > (-180)) {
                        coordinatesEndBuffer.append(currentLongitude2+",");
                    } else {
                        coordinatesEndBuffer.append((180+(currentLongitude2+180))+",");
                    }
                }

            } else {
                coordinatesEndBuffer.append((startLongitude+(a+1)*(longitudeDifference/divider))+",");
            }

            coordinatesEndBuffer.append((startLatitude+(a+1)*(latitudeDifference/divider))+","+(maxAltitude*Math.sin(Math.acos(1 - (a+1)*(1.0/(divider/2.0))))));

            coordinatesElement.addContent(coordinatesStartBuffer.toString());
            coordinatesElement.addContent(" ");
            coordinatesElement.addContent(coordinatesEndBuffer.toString());
            lineStringElement.addContent(coordinatesElement);
            placemarkElement.addContent(lineStringElement);
            folderElement.addContent(placemarkElement);
        }
    }

    private void  addRateWithData(double rateIndicator, double BF, String fromLocation, String toLocation, double startLongitude, double startLatitude, double endLongitude, double endLatitude, int number, Element folderElement) {



        Element placemarkElement  = new Element("Placemark");

        Element placemarkNameElement = new Element("name");
        String name = "rate"+number;
        placemarkNameElement.addContent(name);
        placemarkElement.addContent(placemarkNameElement);

        Element data = new Element("ExtendedData");
        Element schemaData = new Element("SchemaData");
        schemaData.setAttribute("schemaUrl", "#Rate_Schema");
        schemaData.addContent(new Element("SimpleData").setAttribute("name", "Name").addContent(name));
        schemaData.addContent(new Element("SimpleData").setAttribute("name", "From").addContent(fromLocation));
        schemaData.addContent(new Element("SimpleData").setAttribute("name", "To").addContent(toLocation));
        schemaData.addContent(new Element("SimpleData").setAttribute("name", "BF").addContent(Double.toString(BF)));
        schemaData.addContent(new Element("SimpleData").setAttribute("name", "Indicator").addContent(Double.toString(rateIndicator)));
        data.addContent(schemaData);
        placemarkElement.addContent(data);

        Element lineStringElement = new Element("LineString");
        Element altitudeMode = new Element("altitudeMode");
        altitudeMode.addContent("clampToGround");
        lineStringElement.addContent(altitudeMode);
        Element tessellate = new Element("tessellate");
        tessellate.addContent("1");
        lineStringElement.addContent(tessellate);
        Element coordinatesElement = new Element("coordinates");
        coordinatesElement.addContent(startLongitude + "," + startLatitude + " " + endLongitude + "," + endLatitude);
        lineStringElement.addContent(coordinatesElement);
        placemarkElement.addContent(lineStringElement);
        folderElement.addContent(placemarkElement);
    }

    public void outputTextFile(String outFileName) {

//        double[] minMax = new double[2];
//        minMax[0] = supportedRateIndicators[0];
//        minMax[1] = supportedRateIndicators[supportedRateIndicators.length - 1];

        try {
            PrintWriter outFile;
            if (outFileName != null) {
                outFile = new PrintWriter(new FileWriter(outFileName), true);
            } else {
                outFile = new PrintWriter(System.out);
            }
            //sort expected rateIndicator
            if (bayesFactor) {
                outFile.println("Indicator cutoff (for BF = "+cutoff+") = "+getBayesFactorCutOff(cutoff, meanPoissonPrior, offsetPoissonPrior, statesCounter, nonreversible));
            } else {
                outFile.println("Indicator cutoff = "+cutoff);
            }
            outFile.println("mean Poisson Prior = "+meanPoissonPrior);
            outFile.println("Poisson Prior offset = "+offsetPoissonPrior);

            for (int o = 0; o < supportedRateIndicators.length; o++){

                outFile.print(
                        "I="+supportedRateIndicators[o]+"\tBF");

                double BF = getBayesFactor(supportedRateIndicators[o], meanPoissonPrior, statesCounter, offsetPoissonPrior, nonreversible, 0);
                if (BF == Double.POSITIVE_INFINITY) {
                    outFile.print(">"+getBayesFactor(supportedRateIndicators[o], meanPoissonPrior, statesCounter, offsetPoissonPrior, nonreversible, generationCount));
                }  else {
                    outFile.print("="+BF);
                }

                StringBuilder sb = new StringBuilder();
                sb.append(" : between ")
                        .append(supportedLocations[o][0]);
                if (!Double.isNaN(supportedLongitudes[o][0])) {
                    sb.append(" (long: ")
                            .append(supportedLongitudes[o][0])
                            .append("; lat: ")
                            .append(supportedLatitudes[o][0])
                            .append(")");
                }
                sb.append(" and ")
                        .append(supportedLocations[o][1]);
                if (!Double.isNaN(supportedLongitudes[o][1])) {
                    sb.append(" (long: ")
                            .append(supportedLongitudes[o][1])
                            .append("; lat: ")
                            .append(supportedLatitudes[o][1])
                            .append(")");
                }
                outFile.print(sb.toString());
                outFile.println();
            }
            outFile.close();
        } catch(IOException io) {
            System.err.print("Error writing to file: " + outFileName);
        }

    }

    private double[] expectedRateIndicators;
    private boolean nonreversible = false;
    private String[][] locations; // contains start and end locations for transitions
    private String[][] locationNames; // contains start and end locations for transitions
    private double[][] longitudes;  // contains start and end longitudes for transitions
    private double[][] latitudes;  // contains start and end latitudes for transitions
    private String[][] supportedLocations; // contains start and end locations for transitions
    private double[][] supportedLongitudes;  // contains start and end longitudes for transitions
    private double[][] supportedLatitudes;  // contains start and end latitudes for transitions
    private int statesCounter;
    private double[] supportedBFs;
    private double[] supportedRateIndicators;
    private boolean bayesFactor;
    private double cutoff;
    private double meanPoissonPrior;
    private int offsetPoissonPrior;
    private String actualRateString;
    private String relativeRateString;
    private String geoSiteModelString;
    private int generationCount;

    private static double getBayesFactor(double meanIndicator, double meanPoissonPrior, int numberOfLocations, int offset, boolean nonreversible, double generations) {

        double bayesFactor = 0;
        double priorProbabilityDenominator = 0;
        int numberOfRatesMultiplier = 2;
        priorProbabilityDenominator = meanPoissonPrior + offset;
        if (nonreversible) numberOfRatesMultiplier = 1;

        double priorProbability = priorProbabilityDenominator/((numberOfLocations*(numberOfLocations-1))/numberOfRatesMultiplier);

        double priorOdds = priorProbability/(1.0 - priorProbability);
        double posteriorProbability = meanIndicator;
        double posteriorOdds;
        if (generations > 0) {
            posteriorOdds =  (posteriorProbability - (1/generations))/(1 - (posteriorProbability - (1/generations)));
        }  else {
            posteriorOdds =  posteriorProbability/(1 - posteriorProbability);
        }
        bayesFactor = posteriorOdds/priorOdds;

        return bayesFactor;

    }

    private static double getBayesFactorCutOff(double bayesFactor, double meanPoissonPrior, int offset, int numberOfLocations, boolean nonreversible) {

        double bayesFactorCutoff = 0;
        double posteriorOdds = 0;
        int numberOfRatesMultiplier = 2;
        if (nonreversible) numberOfRatesMultiplier = 1;

        double priorProbability = (meanPoissonPrior + offset)/((numberOfLocations*(numberOfLocations-1))/numberOfRatesMultiplier);
        double priorOdds = priorProbability/(1.0 - priorProbability);
        posteriorOdds = priorOdds*bayesFactor;
        bayesFactorCutoff = posteriorOdds/(1.0 + posteriorOdds);

        return bayesFactorCutoff;

    }

    private static double[] meanCol(double[][] x)    {
        double[] returnArray = new double[x[0].length];

        for (int i = 0; i < x[0].length; i++) {

            double m = 0.0;
            int len = 0;

            for (int j = 0; j < x.length; j++) {

                m += x[j][i];
                len += 1;

            }

            returnArray[i] = m / (double) len;

        }
        return returnArray;
    }

    private static void fillRateIndicatorArray(String RateIndicatorLog, double[][] rateIndicators, int burnin, int firstRateIndicator, int numberOfRateIndicators){

        try {
            BufferedReader indicatorReader = new BufferedReader(new FileReader(RateIndicatorLog));
            String rateIndicatorCurrent = indicatorReader.readLine();
            while (rateIndicatorCurrent.startsWith("#")) {
                rateIndicatorCurrent = indicatorReader.readLine();
            }

            // skip the headers in the rateIndicator file
            while (rateIndicatorCurrent.startsWith("state")) {
                rateIndicatorCurrent = indicatorReader.readLine();
            }

            double rateIndicator;
            int linesRead = 0;
            //skip burnin
            while (linesRead < burnin){
                rateIndicatorCurrent = indicatorReader.readLine();
                linesRead ++;
            }

            int rowCounter = 0;

            while (rateIndicatorCurrent!= null && !indicatorReader.equals("")) {

                int columnCounter = 0;

                int startCounter = 1;
                int rateIndicatorCounter = 0;

                StringTokenizer tokens = new StringTokenizer(rateIndicatorCurrent);
                rateIndicator = Double.parseDouble(tokens.nextToken());

                //skip until we encounter rateIndicators
                while (startCounter < firstRateIndicator) {
                    rateIndicator = Double.parseDouble(tokens.nextToken());
                    startCounter ++;

                }

                // read all rateIndicators
                while (rateIndicatorCounter < numberOfRateIndicators) {
                    rateIndicators[rowCounter][columnCounter] = rateIndicator;
                    columnCounter ++;
                    rateIndicatorCounter ++;
                    if (rateIndicatorCounter < numberOfRateIndicators) {
                        rateIndicator = Double.parseDouble(tokens.nextToken());
                    }

                }

                rowCounter ++;
                rateIndicatorCurrent = indicatorReader.readLine();

            }
        } catch (IOException e) {
            System.err.println("Error reading " + RateIndicatorLog);
            System.exit(1);
        }

    }

    private static int getNumberOfEntries(String file, int firstRateIndicator, String rateIndicatorString)    {

        int numberOfRateIndicators = 0;

        try {
            BufferedReader reader = new BufferedReader(new FileReader(file));
            String current3 = reader.readLine();

            //skip comment lines
            while (current3.startsWith("#")) {
                current3 = reader.readLine();
            }

            String rateIndicator = null;
            int startCounter = 1;
            StringTokenizer tokens = new StringTokenizer(current3);
            rateIndicator = tokens.nextToken();
            //find first rateIndicator..
            while (startCounter < firstRateIndicator) {
                rateIndicator = tokens.nextToken();
                startCounter ++;
            }
            // and continue counting from thereon
            while(rateIndicator.contains(rateIndicatorString)) {
                rateIndicator = tokens.nextToken();
                numberOfRateIndicators ++;
            }

        } catch (IOException e) {
            System.err.println("Error reading " + file);
            System.exit(1);
        }

        return numberOfRateIndicators;

    }

    private static boolean hasEntryOf(String file, String entryString)    {

        boolean hasEntry = false;

        try {
            BufferedReader reader2 = new BufferedReader(new FileReader(file));
            String current2 = reader2.readLine();

            //skip comment lines
            while (current2.startsWith("#")) {
                current2 = reader2.readLine();
            }

            String parameter1 = null;
            StringTokenizer tokens1 = new StringTokenizer(current2);

            while(tokens1.hasMoreTokens()) {
                parameter1 = tokens1.nextToken();
                if (parameter1.contains(entryString)) {
                    hasEntry = true;
                }
            }

        } catch (IOException e) {
            System.err.println("Error reading " + file);
            System.exit(1);
        }

        return hasEntry;

    }

    private static int getFirstEntryOf(String file, String entryString)    {

        int firstRateIndicator = 1;

        try {
            BufferedReader reader2 = new BufferedReader(new FileReader(file));
            String current2 = reader2.readLine();

            //skip comment lines
            while (current2.startsWith("#")) {
                current2 = reader2.readLine();
            }

            String parameter1 = null;
            StringTokenizer tokens1 = new StringTokenizer(current2);
            parameter1 = tokens1.nextToken();
            while(!parameter1.contains(entryString)) {
                parameter1 = tokens1.nextToken();
                firstRateIndicator ++;
            }

        } catch (IOException e) {
            System.err.println("Error reading " + file);
            System.exit(1);
        }

        return firstRateIndicator;

    }

    private static int getGenerationCount(String file)    {

        int states = 0;

        try {
            BufferedReader reader1 = new BufferedReader(new FileReader(file));
            String current1 = reader1.readLine();
            while (current1 != null && !reader1.equals("")) {
                while (current1.startsWith("#")) {
                    current1 = reader1.readLine();
                }
                current1 = reader1.readLine();
                states++;
            }

        } catch (IOException e) {
            System.err.println("Error reading " + file);
            System.exit(1);
        }

        return states;

    }

    private static int[] countLinesAndTokens(String coordinatesFileString){
        int lineCounter = 0;
        int tokenCounter = 0;
        int[] container = new int[2];
        try{
            BufferedReader reader1 = new BufferedReader(new FileReader(coordinatesFileString));
            String current1 = reader1.readLine();
            while (current1 != null && !reader1.equals("")) {
                lineCounter++;
                if (lineCounter == 1) {
                    StringTokenizer tokens = new StringTokenizer(current1);
                    while (tokens.hasMoreTokens()) {
                        tokenCounter++;
                        tokens.nextToken();
                    }
                }
                current1 = reader1.readLine();
            }

        } catch (IOException e) {
            System.err.println("Error reading " + coordinatesFileString);
            System.exit(1);
        }
        container[0] = lineCounter;
        container[1] = tokenCounter;
        return container;
    }

    private static void readLocationsCoordinates(String coordinatesFileString, String[][] locationsAndCoordinates){
        try {
            BufferedReader reader2 = new BufferedReader(new FileReader(coordinatesFileString));
            String current2 = reader2.readLine();
            int counter2 = 0;
            while (current2 != null && !reader2.equals("")) {
                StringTokenizer tokens2 = new StringTokenizer(current2);
                for (int i = 0; i < locationsAndCoordinates[0].length; i++) {
                    locationsAndCoordinates[counter2][i] = tokens2.nextToken();
                    progressStream.print(locationsAndCoordinates[counter2][i]+"\t");
                }
                progressStream.print("\r");
                counter2 ++;
                current2 = reader2.readLine();
            }
        } catch (IOException e) {
            e.printStackTrace();
            return;
        }
    }

    public static void main(String[] args) throws IOException {

        String inputFileName = null;
        String outputFileName = null;

        String locationsFileName = null;
        String[][] locations = null;
        boolean kml = false;
        String lowerLinkColor = "FFFFFF"; //red: 0000FF green: 00FF00 magenta: FF00FF white: FFFFFF yellow: 00FFFF cyan: FFFF00
        String upperLinkColor = "FF00FF";
        String KMLoutputFile = "KMLrates.kml";
        double branchWidthConstant = 2.5;
        double branchWidthMultiplier = 7.0;
        double altitudeFactor = 500;
        //Double width = 3.0;

        int burnin = -1;
        double meanPoissonPrior = 0.693;
        int offsetPoissonPrior = 0;
        int numberOfStates = 0;

        double cutoff = 3.0;
        boolean bayesFactor = true; // if false, we will use an indicator cut off value

        boolean rateSummary = false;
        String rateIndicatorString	= "indicators";
        String actualRateString = "productStatistic";
        String relativeRateString = "rates";
        //this is for rate (dist/time) summaries
        String geoSiteModelString = "geoSiteModel";

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.IntegerOption(BURNIN, "the number of states to be considered as 'burn-in' [default = 0]"),
                        new Arguments.StringOption(LOCATIONSFILE,"coordinates file","a file with latitudes and longitudes for each location (required for a kml output)"),
                        //boolean for KML
                        new Arguments.StringOption(KML, falseTrue, false,
                                "generate a KML file including well-supported rates [default = false]"),
                        new Arguments.IntegerOption(LOCATIONSTATES,"the number of locations states used in the analyses [requires a coordinates file if not provided]"),
                        new Arguments.IntegerOption(POFFSET,"the offset of the (truncated) Poisson prior [default=locations-1]"),
                        new Arguments.RealOption(PMEAN,"the mean of the (truncated) Poisson prior  [default=0.693 (log2)]"),
                        new Arguments.RealOption(BFCUTOFF,"the Bayes Factor values above which we consider rates to be well supported  [default=3.0]"),
                        new Arguments.RealOption(ICUTOFF,"the indicator values above which we consider rates to be well supported  [default uses a Bayes factor cut off of 3.0]"),
                        new Arguments.StringOption(ISTRING, "indicator_string", "prefix string used for outputting the rate indicators in the log file [default = indicators]"),
                        new Arguments.StringOption(RSTRING, "relativeRate_string", "prefix string used for outputting the relative rates in the log file [default = rates]"),
                        new Arguments.StringOption(PSTRING, "rate*indicator_string", "prefix string used for outputting the product statistic for rates*indicators [default = productStatistic]"),
                        new Arguments.StringOption(GSTRING, "geoSiteModel.mu_string", "string used for outputting geoSiteModel.mu in the log file [default = geoSiteModel]"),
                        new Arguments.StringOption(KMLFILE,"KML output file","KML output file name [default=KMLrates.kml]"),
                        new Arguments.StringOption(LOWCOLOR, "lower link strength color", "specifies an lower link color for the links [default=FF00FF]"),
                        new Arguments.StringOption(UPCOLOR, "upper link strength color", "specifies an upper link color for the links [default=FFFF00]"),
                        new Arguments.RealOption(BWC,"specifies the connection (rate) width constant [default=2.5]"),
                        new Arguments.RealOption(BWM,"specifies the connection (rate)  width multiplier [default=7.0]"),
                        new Arguments.RealOption(ALTITUDE,"specifies the altitudefactor for the connections (rate) [default=500]"),
                        //new Arguments.RealOption(WIDTH,"width for KML rates [default=3.0]"),
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

        // Make sense of arguments

        if (arguments.hasOption(BURNIN)) {
            burnin = arguments.getIntegerOption(BURNIN);
        }
        progressStream.println("Ignoring "+burnin+" states as burn-in");

        if (arguments.hasOption(LOCATIONSTATES)) {
            numberOfStates = arguments.getIntegerOption(LOCATIONSTATES);
        }

        locationsFileName = arguments.getStringOption(LOCATIONSFILE);
        if (locationsFileName != null) {
            int counts[] = countLinesAndTokens(locationsFileName);
            //System.out.println(counts[0]+"\t"+counts[1]);
            //read in locations
            if (numberOfStates > 0) {
                if (numberOfStates != counts[0]) {
                    System.err.println("number of states provided ("+numberOfStates+") does not match lines in coordinates file ("+counts[0]+") ??");
                }
            }
            locations = new String[counts[0]][counts[1]];
            readLocationsCoordinates(locationsFileName,locations);

            if (numberOfStates == 0) {
                numberOfStates = counts[0];
            }

        } else {
            if (numberOfStates == 0) {
                System.err.println("no states provided, nor coordinates file ??");
            }
        }

        String kmlBooleanString = arguments.getStringOption(KML);
        if (kmlBooleanString != null && kmlBooleanString.compareToIgnoreCase("true") == 0) {
            kml = true;
            if (locationsFileName == null) {
                System.err.println("you want a KML file without a coordinates file??");
            }
        }

        String kmlFileName = arguments.getStringOption(KMLFILE);
        if (kmlFileName != null) {
            KMLoutputFile = kmlFileName;
        }

        if (arguments.hasOption(PMEAN)) {
            meanPoissonPrior = arguments.getRealOption(PMEAN);
            progressStream.println("Poisson prior with mean "+meanPoissonPrior);
        }  else {
            progressStream.println("Poisson prior with mean "+meanPoissonPrior+" (default)");
        }

        if (arguments.hasOption(POFFSET)) {
            offsetPoissonPrior = arguments.getIntegerOption(POFFSET);
            progressStream.println("Poisson offset = "+offsetPoissonPrior);
        }   else {
            offsetPoissonPrior = numberOfStates - 1;
            progressStream.println("Poisson offset = "+offsetPoissonPrior+" (locations - 1)");
        }

        if (arguments.hasOption(BFCUTOFF)) {
            cutoff = arguments.getRealOption(BFCUTOFF);
        }
        if (arguments.hasOption(ICUTOFF)) {
            cutoff = arguments.getRealOption(ICUTOFF);
            bayesFactor = false;
        }
        if (bayesFactor) {
            progressStream.println("Bayes factor cutoff = "+cutoff);
        }   else {
            progressStream.println("indicator factor cutoff = "+cutoff);
        }


        if (arguments.hasOption(BWC)) {
            branchWidthConstant = arguments.getRealOption(BWC);
        }

        if (arguments.hasOption(BWM)) {
            branchWidthMultiplier = arguments.getRealOption(BWM);
        }

        if (arguments.hasOption(ALTITUDE)) {
            altitudeFactor = arguments.getRealOption(ALTITUDE);
        }

        String indicatorString = arguments.getStringOption(ISTRING);
        if (indicatorString != null) {
            rateIndicatorString = indicatorString;
        }

        String rateString = arguments.getStringOption(PSTRING);
        if (rateString != null) {
            actualRateString = rateString;
        }

        String relRateString = arguments.getStringOption(RSTRING);
        if (relRateString != null) {
            relativeRateString = relRateString;
        }

        String geoString = arguments.getStringOption(GSTRING);
        if (geoString != null) {
            geoSiteModelString = geoString;
        }

        String color1String = arguments.getStringOption(LOWCOLOR);
        if (color1String != null) {
            lowerLinkColor = color1String;
            if (locationsFileName == null) {
                System.err.print("color string but no coordinates file for KML output??");
            }
        }

        String color2String = arguments.getStringOption(UPCOLOR);
        if (color2String != null) {
            upperLinkColor = color2String;
            if (locationsFileName == null) {
                System.err.print("color string but no coordinates file for KML output??");
            }
        }


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

        RateIndicatorBF rateIndicatorBF = new RateIndicatorBF(inputFileName,burnin,rateIndicatorString,numberOfStates,locations,bayesFactor,cutoff,meanPoissonPrior,offsetPoissonPrior,actualRateString,relativeRateString,geoSiteModelString);
        rateIndicatorBF.outputTextFile(outputFileName);
        if (kml) {
            rateIndicatorBF.outputKML(KMLoutputFile,lowerLinkColor,upperLinkColor, branchWidthConstant, branchWidthMultiplier, altitudeFactor);
        }
        System.exit(0);


    }
}

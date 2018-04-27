/*
 * TimeSlicer.java
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

package dr.app.tools;

import dr.app.beast.BeastVersion;
import dr.app.util.Arguments;
import dr.evolution.io.Importer;
import dr.evolution.io.NewickImporter;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.TreeImporter;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeUtils;
import dr.geo.KMLCoordinates;
import dr.geo.KernelDensityEstimator2D;
import dr.geo.Polygon2D;
import dr.geo.contouring.*;
import dr.geo.math.SphericalPolarCoordinates;
import dr.inference.trace.TraceDistribution;
import dr.inference.trace.TraceType;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.util.HeapSort;
import dr.util.Version;
import org.jdom.Element;
import org.jdom.output.Format;
import org.jdom.output.XMLOutputter;

import java.awt.geom.Point2D;
import java.io.*;
import java.text.SimpleDateFormat;
import java.util.*;

/**
 * @author Marc A. Suchard
 * @author Philippe Lemey
 *         <p/>
 *         example for location slices through time:
 *         -burnin500 -trait location -sliceFileHeights sliceTimes -summary true -mrsd 2007.63 -format KML -progress true -impute true -noise true -HPD 95 WNV_relaxed_gamma_geo.trees output.kml
 *         <p/>
 *         example for Markov reward o backbone summary (note we cannot imput here):
 *         -burnin 50 -trait AC1_R,AC2_R,AC3_R,AC4_R,AC5_R,AC6_R,AC7_R,AC8_R,AC9_R,AC10_R,AC11_R,AC12_R,AC13_R,AC14_R -sliceTimes 2002,2002.25,2002.5,2002.75,2003,2003.25,2003.5,2003.75,2004,2004.25,2004.5,2004.75,2005,2005.25,2005.5,2005.75,2006 -summary true -mrsd 2007 -branchNorm true -format TAB -progress true -backbonetaxa december2006taxa.txt -branchset backbone airCommunitiesMM_rewards.sub.trees TSRewardOutput.txt
 *         example of 1D antigenic summary
 *         -burnin0 -trait antigenic -sliceFileTimes timeInterval -summary true -mrsd 2007 -format TAB -progress true -impute true -noise false -HPD 95 H3N2_antigenic_sub.trees antigenic.txt
 *
 *         other example
 *         -burnin 100 -trait location -sliceHeights 0,5,15 -summary true -mrsd 2004.7 -rateAttribute location.rate -format KML -progress true -impute true -noise true -HPD 80 -sdr sliceDispersalRates.txt RacRABV_Gamma.trees test.kml
 */

public class TimeSlicer {

    public static final String sep = "\t";
    public static final String PRECISION_STRING = "precision";
    public static final String RATE_ATTRIBUTE = "rateAttribute";

    public static final String SLICE_ELEMENT = "slice";
    public static final String REGIONS_ELEMENT = "hpdRegion";
    public static final String NODE_ELEMENT = "node";
    public static final String ROOT_ELEMENT = "root";
    public static final String TRAIT = "trait";
    public static final String NAME = "name";
    public static final String DENSITY_VALUE = "density";
    public static final String SLICE_VALUE = "time";

    public static final String STYLE = "Style";
    public static final String ID = "id";
    public static final String WIDTH = "0.5";
    public static final String startHPDColor = "00F1D6"; //blue=B36600 yellow="00F1D6"
    public static final String endHPDColor = "00FF00"; //red=0000FF green="00FF00"
    public static final String opacity = "6f";

    public static final String BURNIN = "burnin";
    public static final String SKIP = "skip";
    public static final String SLICE_TIMES = "sliceTimes";
    public static final String SLICE_HEIGHTS = "sliceHeights";
    public static final String SLICE_COUNT = "sliceCount";
    public static final String START_TIME = "startTime";
    public static final String SLICE_FILE_HEIGHTS = "sliceFileHeights";
    public static final String SLICE_FILE_TIMES = "sliceFileTimes";
    public static final String SLICE_MODE = "sliceMode";
    public static final String ROOT = "root";
    public static final String TIPS = "tips";
    public static final String CONTOURS = "contours";
    public static final String POINTS = "points";
    public static final String MRSD = "mrsd";
    public static final String HELP = "help";
    public static final String NOISE = "noise";
    public static final String IMPUTE = "impute";
    public static final String SUMMARY = "summary";
    public static final String FORMAT = "format";
    public static final String CONTOUR_MODE = "contourMode";
    public static final String NORMALIZATION = "normalization";
    public static final String HPD = "hpd";
    public static final String SDR = "sdr";
    public static final String SNR = "snr";
    public static final String PROGRESS = "progress";
    public static final double treeLengthPercentage = 0.00;
    public static final String BRANCH_NORMALIZE = "branchnorm";
    public static final String BRANCHSET = "branchset";
    public static final String BACKBONETAXA = "backbonetaxa";
    public static final String CLADETAXA = "cladetaxa";
    public static final String LATMAX = "latmax";
    public static final String LATMIN = "latmin";
    public static final String LONGMAX = "longmax";
    public static final String LONGMIN = "longmin";
    public static final String ICON = "http://maps.google.com/mapfiles/kml/pal4/icon49.png";
    public static final String GRIDSIZE = "gridsize";
    public static final double[] BANDWIDTHS = new double[]{1.0,1.0};
    public static final boolean BANDWIDTHLIMIT = true;
    public static final boolean GREATCIRCLEDISTANCE = true;
    public static final String SUBSTITUTION = "N";
    public static final String DESCENDENTS = "descendents";

    public static final String[] falseTrue = {"false", "true"};

    private final static Calendar calendar = GregorianCalendar.getInstance();
    private final static SimpleDateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd");

    public TimeSlicer(String treeFileName, int burnin, int skipEvery, String[] traits, double[] sliceHeights, boolean impute,
                      boolean trueNoise, double mrsd, ContourMode contourMode, SliceMode sliceMode,
                      final boolean summarizeRoot, final boolean summarizeTips, Normalization normalize, boolean getSDR, boolean getSNR,
                      String progress, boolean branchNormalization, BranchSet branchset, Set taxaSet, int grid,
                      double latMin, double latMax, double longMin, double longMax, Set descendentTaxaSet, String rateString) {

        this.traits = traits;
        traitCount = traits.length;

        sliceCount = 1;
        doSlices = false;
        mostRecentSamplingDate = mrsd;
        this.contourMode = contourMode;
        this.sliceMode = sliceMode;
        sdr = getSDR;
        snr = getSNR;
        if (sdr && snr){
            System.err.println("both SDR and SNR are requested; only one can be summarized.. prioritizing SDR");
        }

        this.latMin = latMin;
        this.latMax = latMax;
        this.longMin = longMin;
        this.longMax = longMax;
//        System.out.println("latMin: "+latMin);
//        System.out.println("latMax: "+latMax);
//        System.out.println("longMin: "+longMin);
//        System.out.println("longMax: "+longMax);
        this.descendentTaxaSet = descendentTaxaSet;

        rateAttributeString = rateString;

        gridSize = grid;

        if (progress != null) {
            if (progress.equalsIgnoreCase("true")) {
                sliceProgressReport = true;
            } else if (progress.equalsIgnoreCase("check")) {
                sliceProgressReport = true;
                checkSliceContours = true;
            }
        }

        if (sliceHeights != null) {
            sliceCount = sliceHeights.length;
            doSlices = true;

            this.sliceHeights = sliceHeights;

            if ((mostRecentSamplingDate - sliceHeights[sliceHeights.length - 1]) < 0) {
                ancient = true;
            }
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
        if (summarizeRoot) {
            rootValues = new ArrayList<List<Trait>>(traitCount);
            for (int k = 0; k < traitCount; k++) {
                List<Trait> thisTrait = new ArrayList<Trait>();
                rootValues.add(thisTrait);
            }
        }
        if (summarizeTips) {
            tipValues = new ArrayList<List<List<Trait>>>();
            tipNames = new ArrayList<String>();
        }

        try {
            readAndAnalyzeTrees(treeFileName, burnin, skipEvery, traits, sliceHeights, impute, trueNoise, normalize, branchNormalization, branchset, taxaSet);
        } catch (IOException e) {
            System.err.println("Error reading file: " + treeFileName);
            System.exit(-1);
        } catch (Importer.ImportException e) {
            System.err.println("Error parsing trees in file: " + treeFileName);
            System.exit(-1);
        }

//        if (values.get(0).get(0).size() == 0) {
//            System.err.print("Trait(s), ");
//            for (int a = 0; a < traits.length; a++){
//                System.err.print(traits[a]);
//            }
//            System.err.println(",values missing from trees.");
//
//            System.exit(-1);
//        }

        progressStream.println(treesRead + " trees read.");
        progressStream.println(treesAnalyzed + " trees analyzed.");

    }


    private Element rootElement;
    private Element documentElement;
    private Element contourFolderElement;
    private Element pointsFolderElement;
    private Element nodeFolderElement;
    private StringBuffer tabOutput = new StringBuffer();

    public void output(String outFileName, boolean summaryOnly, final boolean summarizeRoot, final boolean summarizeTips, boolean contours, boolean points, OutputFormat outputFormat, double[] hpdValues, String sdrFile, String snrFile) {

        resultsStream = System.out;

        if (outFileName != null) {
            try {
                resultsStream = new PrintStream(new File(outFileName));
            } catch (IOException e) {
                System.err.println("Error opening file: " + outFileName);
                System.exit(-1);
            }
        }

        if (!summaryOnly) {
            outputHeader(traits);

            if (sliceHeights == null || sliceHeights.length == 0) {
                outputSlice(0, Double.NaN);
            } else {
                for (int i = 0; i < sliceHeights.length; i++) {
                    outputSlice(i, sliceHeights[i]);
                }

            }
        } else { // Output summaries

            if (outputFormat == OutputFormat.XML) {
                rootElement = new Element("xml");

            } else if (outputFormat == OutputFormat.KML) {

                Element hpdSchema = new Element("Schema");
                hpdSchema.setAttribute("id", "HPD_Schema");
                hpdSchema.addContent(new Element("SimpleField")
                        .setAttribute("name", "Name")
                        .setAttribute("type", "string")
                        .addContent(new Element("displayName").addContent("Name")));
                hpdSchema.addContent(new Element("SimpleField")
                        .setAttribute("name", "Description")
                        .setAttribute("type", "string")
                        .addContent(new Element("displayName").addContent("Description")));
                hpdSchema.addContent(new Element("SimpleField")
                        .setAttribute("name", "Time")
                        .setAttribute("type", "double")
                        .addContent(new Element("displayName").addContent("Time")));
                hpdSchema.addContent(new Element("SimpleField")
                        .setAttribute("name", "Height")
                        .setAttribute("type", "double")
                        .addContent(new Element("displayName").addContent("Height")));
                hpdSchema.addContent(new Element("SimpleField")
                        .setAttribute("name", "HPD")
                        .setAttribute("type", "double")
                        .addContent(new Element("displayName").addContent("HPD")));

                Element nodeSchema = new Element("Schema");
                nodeSchema.setAttribute("id", "Point_Schema");
                nodeSchema.addContent(new Element("SimpleField")
                        .setAttribute("name", "Name")
                        .setAttribute("type", "string")
                        .addContent(new Element("displayName").addContent("Name")));
                nodeSchema.addContent(new Element("SimpleField")
                        .setAttribute("name", "Description")
                        .setAttribute("type", "string")
                        .addContent(new Element("displayName").addContent("Description")));
                nodeSchema.addContent(new Element("SimpleField")
                        .setAttribute("name", "Time")
                        .setAttribute("type", "double")
                        .addContent(new Element("displayName").addContent("Time")));
                nodeSchema.addContent(new Element("SimpleField")
                        .setAttribute("name", "Height")
                        .setAttribute("type", "double")
                        .addContent(new Element("displayName").addContent("Height")));

                if (contours) {
                    contourFolderElement = new Element("Folder");
                    Element contourFolderNameElement = new Element("name");
                    contourFolderNameElement.addContent("surface HPD regions");
                    //rootElement.setAttribute("xmlns","http://earth.google.com/kml/2.2");
                    contourFolderElement.addContent(contourFolderNameElement);
                }

                if (points) {
                    pointsFolderElement = new Element("Folder");
                    Element pointsFolderNameElement = new Element("name");
                    pointsFolderNameElement.addContent("points");
                    //rootElement.setAttribute("xmlns","http://earth.google.com/kml/2.2");
                    pointsFolderElement.addContent(pointsFolderNameElement);
                }

                if (sliceMode == TimeSlicer.SliceMode.NODES) {
                    nodeFolderElement = new Element("Folder");
                    Element nodeFolderNameElement = new Element("name");
                    nodeFolderNameElement.addContent("nodes");
                    //rootElement.setAttribute("xmlns","http://earth.google.com/kml/2.2");
                    nodeFolderElement.addContent(nodeFolderNameElement);
                }
                Element documentNameElement = new Element("name");
                String documentName = outFileName;
                if (documentName == null)
                    documentName = "default";
                if (documentName.endsWith(".kml"))
                    documentName = documentName.replace(".kml", "");
                documentNameElement.addContent(documentName);

                documentElement = new Element("Document");
                documentElement.addContent(documentNameElement);
                documentElement.addContent(hpdSchema);
                documentElement.addContent(nodeSchema);
                if (contourFolderElement != null) {
                    documentElement.addContent(contourFolderElement);
                }
                if (pointsFolderElement != null) {
                    documentElement.addContent(pointsFolderElement);
                }
                if (nodeFolderElement != null) {
                    documentElement.addContent(nodeFolderElement);
                }

                rootElement = new Element("kml");
                rootElement.addContent(documentElement);
            }

            if (sliceHeights == null) {
                for (double hpdValue : hpdValues) {
                    summarizeSlice(0, Double.NaN, contours, points, outputFormat, hpdValue);
                }
            } else {
                if (outputFormat == OutputFormat.TAB) {
                    if (mostRecentSamplingDate > 0) {
                        tabOutput.append("trait\t" + "sliceTime\t" + "mean\t" + "stdev\t" + "HPDlow\t" + "HPDup");
                    } else {
                        tabOutput.append("trait\t" + "sliceHeight\t" + "mean\t" + "stdev\t" + "HPDlow\t" + "HPDup");
                    }
                }
                for (int i = 0; i < sliceHeights.length; i++) {
                    for (double hpdValue : hpdValues) {
                        summarizeSlice(i, sliceHeights[i], contours, points, outputFormat, hpdValue);
                    }
                }
            }

            if (summarizeRoot) {
                for (double hpdValue : hpdValues) {
                    summarizeRoot(contours, points, outputFormat, hpdValue);
                }
            }

            if (summarizeTips) {
                for (double hpdValue : hpdValues) {
                    summarizeTips(contours, points, outputFormat, hpdValue);
                }
            }

            if (outputFormat == OutputFormat.TAB) {
                resultsStream.println(tabOutput);
            } else {
                XMLOutputter xmlOutputter = new XMLOutputter(Format.getPrettyFormat().setTextMode(Format.TextMode.PRESERVE));
                try {
                    xmlOutputter.output(rootElement, resultsStream);
                } catch (IOException e) {
                    System.err.println("IO Exception encountered: " + e.getMessage());
                    System.exit(-1);
                }
            }

        }
// writes out dispersal rate summaries for the whole tree, there is a beast xml statistic that can do this now
//        if (containsLocation) {
//            try{
//                PrintWriter dispersalRateFile = new PrintWriter(new FileWriter("dispersalRates.log"), true);
//                dispersalRateFile.print("state"+"\t"+"dispersalRate(native units)"+"\t"+"dispersalRate(great circle distance, km)\r");
//                for(int x = 0; x < dispersalrates.size(); x++ ) {
//                    dispersalRateFile.print(x+"\t"+dispersalrates.get(x)+"\r");
//                }
//
//            } catch (IOException e) {
//                System.err.println("IO Exception encountered: "+e.getMessage());
//                System.exit(-1);
//            }
//        }


        if (sdr || snr) {
            double[][] sliceTreeWeightedAverageRates = new double[sliceTreeDistanceArrays.size()][sliceCount];
            double[][] sliceTreeDistancesToSummarize = new double[sliceTreeDistanceArrays.size()][sliceCount];
//          can be useful for E[N]-E[S]  stuff
            double[][] sliceTreeIntervalDistances = new double[sliceTreeDistanceArrays.size()][sliceCount];

            double[][] sliceTreeMaxPathRates = new double[sliceTreeDistanceArrays.size()][sliceCount];
            double[][] sliceTreeMaxDistanceFromRootRates = new double[sliceTreeDistanceArrays.size()][sliceCount];
             double[][] sliceTreeDistances = new double[sliceCount][sliceTreeDistanceArrays.size()];
            double[][] sliceTreeTimes = new double[sliceCount][sliceTreeDistanceArrays.size()];
            //double[][] sliceTreeMaxDistances = new double[sliceCount][sliceTreeMaxDistanceArrays.size()];
            double[][] sliceTreeMaxDistancesFromRoot = new double[sliceTreeMaxDistanceFromRootArrays.size()][sliceCount];
            double[][] sliceTreeMaxPathDistances = new double[sliceTreeMaxPathDistanceArrays.size()][sliceCount];
            double[][] sliceTreeTimesFromRoot = new double[sliceCount][sliceTreeTimeFromRootArrays.size()];
            double[][] sliceTreeDiffusionCoefficients = new double[sliceTreeDiffusionCoefficientArrays.size()][sliceCount];
            double[][] sliceTreeDiffusionCoefficientVariances = new double[sliceTreeDiffusionCoefficientVarianceArrays.size()][sliceCount];
            //double[][] sliceTreeWeightedAverageDiffusionCoefficients = new double[sliceTreeDistanceArrays.size()][sliceCount];
            for (int q = 0; q < sliceTreeDistanceArrays.size(); q++) {
                double[] distanceArray = (double[]) sliceTreeDistanceArrays.get(q);
                double[] timeArray = (double[]) sliceTreeTimeArrays.get(q);

                double[] maxDistanceFromRootArray = null;
                double[] maxPathDistanceArray = null;
                double[] timeFromRootArray = null;
                double[] diffusionCoefficientArray = null;
                double[] diffusionCoefficientVarianceArray = null;

                if (sdr) {
                    maxPathDistanceArray = (double[]) sliceTreeMaxPathDistanceArrays.get(q);
                    maxDistanceFromRootArray = (double[]) sliceTreeMaxDistanceFromRootArrays.get(q);
                    timeFromRootArray = (double[]) sliceTreeTimeFromRootArrays.get(q);
                    diffusionCoefficientArray = (double[]) sliceTreeDiffusionCoefficientArrays.get(q);
                    diffusionCoefficientVarianceArray = (double[]) sliceTreeDiffusionCoefficientVarianceArrays.get(q);

                }

                for (int r = 0; r < distanceArray.length; r++) {
                    sliceTreeDistances[r][q] = distanceArray[r];
                    sliceTreeTimes[r][q] = timeArray[r];
                    if (sdr){
                        sliceTreeMaxPathDistances[q][r] = maxPathDistanceArray[r];
                        sliceTreeMaxDistancesFromRoot[q][r] = maxDistanceFromRootArray[r];
                        sliceTreeTimesFromRoot[r][q] = timeFromRootArray[r];
                        sliceTreeDiffusionCoefficients[q][r] = diffusionCoefficientArray[r];
                        sliceTreeDiffusionCoefficientVariances[q][r] = diffusionCoefficientVarianceArray[r];
                    }
                }
            }

            //print2DArray(sliceTreeDistances,"sliceTreeDistances.txt");
            //print2DArray(sliceTreeTimes,"sliceTreeTimes.txt");
            print2DTransposedArray(sliceTreeDiffusionCoefficients,"sliceTreeDiffusionCoefficients.txt");
            print2DTransposedArray(sliceTreeDiffusionCoefficientVariances,"sliceTreeDiffusionCoefficientVariances.txt");
            //print2DArray(sliceTreeTimesFromRoot,"sliceTreeTimesFromRoot.txt");

            if (sliceCount > 1) {
                for (int s = 0; s < sliceTreeDistanceArrays.size(); s++) {
                    double[] distanceArray = (double[]) sliceTreeDistanceArrays.get(s);
                    double[] timeArray = (double[]) sliceTreeTimeArrays.get(s);
                    double[] maxDistanceFromRootArray = null;
                    double[] maxPathDistanceArray = null;
                    double[] timeFromRoot = null;
                    if (sdr){
                        maxPathDistanceArray = (double[]) sliceTreeMaxPathDistanceArrays.get(s);
                        maxDistanceFromRootArray = (double[]) sliceTreeMaxDistanceFromRootArrays.get(s);
                        timeFromRoot = (double[]) sliceTreeTimeFromRootArrays.get(s);
                    }

                    for (int t = 0; t < (sliceCount - 1); t++) {
                        if (sdr){
                            sliceTreeMaxDistanceFromRootRates[s][t] = (maxDistanceFromRootArray[t]) / (timeFromRoot[t]);
                            sliceTreeMaxPathRates[s][t] = (maxPathDistanceArray[t]) / (timeFromRoot[t]);
                        }
                        if ((timeArray[t] - timeArray[t + 1]) > ((Double) treeLengths.get(s) * treeLengthPercentage)) {
                            sliceTreeWeightedAverageRates[s][t] = (distanceArray[t] - distanceArray[t + 1]) / (timeArray[t] - timeArray[t + 1]);
                            sliceTreeDistancesToSummarize[s][t] = (distanceArray[t]);
                            sliceTreeIntervalDistances[s][t] = (distanceArray[t] - distanceArray[t + 1]);
                            //sliceTreeWeightedAverageRates[s][t] = (sliceTreeDistances[t][s] - sliceTreeDistances[t+1][s])/(sliceTreeTimes[t][s] - sliceTreeTimes[t+1][s]);
                            //sliceTreeWeightedAverageDiffusionCoefficients[s][t] = Math.pow((distanceArray[t] - distanceArray[t+1]),2.0)/(4.0*(timeArray[t] - timeArray[t+1]));
                        } else {
                            if (timeArray[t] > 0) {
                                if (t == 0) {
                                    throw new RuntimeException("Philippe to fix: the time slices are expected in ascending height order");
                                }
                                sliceTreeWeightedAverageRates[s][t] = sliceTreeWeightedAverageRates[s][t - 1];
                                sliceTreeDistancesToSummarize[s][t] = sliceTreeDistancesToSummarize[s][t - 1];
                                sliceTreeIntervalDistances[s][t] = sliceTreeIntervalDistances[s][t - 1];
                                //sliceTreeWeightedAverageDiffusionCoefficients[s][t] = sliceTreeWeightedAverageDiffusionCoefficients[s][t-1];
                            } else {
                                //set it to NaN, we ignore NaNs when getting summary stats
                                sliceTreeWeightedAverageRates[s][t] = Double.NaN;
                                sliceTreeDistancesToSummarize[s][t] = Double.NaN;
                                sliceTreeIntervalDistances[s][t] = Double.NaN;
                                //sliceTreeWeightedAverageDiffusionCoefficients[s][t] = Double.NaN;
                            }
                        }
                    }

                    if ((timeArray[sliceCount - 1]) > ((Double) treeLengths.get(s) * treeLengthPercentage)) {
                        sliceTreeWeightedAverageRates[s][sliceCount - 1] = (distanceArray[sliceCount - 1]) / (timeArray[sliceCount - 1]);
                        sliceTreeDistancesToSummarize[s][sliceCount - 1] = (distanceArray[sliceCount - 1]);
                        sliceTreeIntervalDistances[s][sliceCount - 1] = (distanceArray[sliceCount - 1]);
                        //sliceTreeWeightedAverageDiffusionCoefficients[s][sliceCount-1] = Math.pow(distanceArray[sliceCount-1],2.0)/(4.0*timeArray[sliceCount-1]);
                    } else {
                        if ((timeArray[sliceCount - 1]) > 0) {
                            sliceTreeWeightedAverageRates[s][sliceCount - 1] = sliceTreeWeightedAverageRates[s][sliceCount - 2];
                            sliceTreeDistancesToSummarize[s][sliceCount - 1] = sliceTreeDistancesToSummarize[s][sliceCount - 2];
                            sliceTreeIntervalDistances[s][sliceCount - 1] = sliceTreeIntervalDistances[s][sliceCount - 2];
                                  //sliceTreeWeightedAverageDiffusionCoefficients[s][sliceCount-1] = sliceTreeWeightedAverageDiffusionCoefficients[s][sliceCount-2];
                        } else {
                            sliceTreeWeightedAverageRates[s][sliceCount - 1] = (distanceArray[sliceCount - 1]) / (timeArray[sliceCount - 1]);
                            sliceTreeDistancesToSummarize[s][sliceCount - 1] = (distanceArray[sliceCount - 1]);
                            sliceTreeIntervalDistances[s][sliceCount - 1] = distanceArray[sliceCount - 1];
                            //sliceTreeWeightedAverageDiffusionCoefficients[s][sliceCount-1] = Math.pow(distanceArray[sliceCount-1],2.0)/(4.0*timeArray[sliceCount-1]);
                        }
                    }
                    if(sdr){
                        sliceTreeMaxPathRates[s][sliceCount - 1] = maxPathDistanceArray[sliceCount - 1] / timeFromRoot[sliceCount - 1];
                        sliceTreeMaxDistanceFromRootRates[s][sliceCount - 1] = maxDistanceFromRootArray[sliceCount - 1] / timeFromRoot[sliceCount - 1];
                     }
                }
            } else {
                for (int s = 0; s < sliceTreeDistanceArrays.size(); s++) {
                    double[] distanceArray = (double[]) sliceTreeDistanceArrays.get(s);
                    double[] timeArray = (double[]) sliceTreeTimeArrays.get(s);
                    sliceTreeWeightedAverageRates[s][0] = distanceArray[0] / timeArray[0];
                    sliceTreeDistancesToSummarize[s][0] = distanceArray[0];
                    sliceTreeIntervalDistances[s][0] = distanceArray[0];
                    //sliceTreeWeightedAverageDiffusionCoefficients[s][0] = Math.pow(distanceArray[0],2.0)/(4.0*timeArray[0]);
                    if (sdr){
                        double[] maxPathDistanceArray = (double[]) sliceTreeMaxPathDistanceArrays.get(s);
                        double[] maxDistanceFromRootArray = (double[]) sliceTreeMaxDistanceFromRootArrays.get(s);
                        double[] timeFromRoot = (double[]) sliceTreeTimeFromRootArrays.get(s);
                        sliceTreeMaxPathRates[s][0] = maxPathDistanceArray[0] / timeFromRoot[0];
                        sliceTreeMaxDistanceFromRootRates[s][0] = maxDistanceFromRootArray[0] / timeFromRoot[0];
                     }
                }
            }

            //print2DArray(sliceTreeDistancesToSummarize,"sliceTreeDistancesToSummarize.txt");
            //print2DArray(sliceTreeWeightedAverageRates,"sliceTreeWeightedAverageRates.txt");
            //print2DArray(sliceTreeIntervalDistances,"sliceTreeIntervalDistances.txt");


            try {
                String outfile = sdrFile;
                if (!sdr) {
                    outfile = snrFile;
                }
                PrintWriter sliceDispersalRateFile = new PrintWriter(new FileWriter(outfile), true);
                sliceDispersalRateFile.print("sliceTime" + "\t");
                if (mostRecentSamplingDate > 0) {
                    sliceDispersalRateFile.print("realTime" + "\t");
                }
                sliceDispersalRateFile.print("mean dispersalRate" + "\t" + "hpd low" + "\t" + "hpd up");
                if (sdr){
                    sliceDispersalRateFile.print("\t" +"mean wavefrontRate" + "\t" + "hpd low" + "\t" + "hpd up" + "\t" + "mean wavefrontDistance" + "\t" + "hpd low" + "\t" + "hpd up" + "\t" + "mean cumulative DiffusionCoefficient" + "\t" + "hpd low" + "\t" + "hpd up" + "\r");
                }  else if (snr){
                    sliceDispersalRateFile.print("\r");
                }

                double[] meanWeightedAverageDispersalRates = meanColNoNaN(sliceTreeWeightedAverageRates);
                double[][] hpdWeightedAverageDispersalRates = getArrayHPDintervals(sliceTreeWeightedAverageRates);
//                double[] meanSliceTreeDistancesToSummarize = meanColNoNaN(sliceTreeDistancesToSummarize);
//                double[][] hpdSliceTreeDistancesToSummarize = getArrayHPDintervals(sliceTreeDistancesToSummarize);

//                double[] meanMaxDispersalDistances = null;
//                double[][] hpdMaxDispersalDistances = null;
//                double[] meanMaxDispersalRates = null;
//                double[][] hpdMaxDispersalRates = null;

                double[] meanMaxPathDispersalDistances = null;
                double[][] hpdMaxPathDispersalDistances = null;
                double[] meanMaxPathDispersalRates = null;
                double[][] hpdMaxPathDispersalRates = null;

                double[] meanMaxRootDispersalDistances = null;
                double[][] hpdMaxRootDispersalDistances = null;
                double[] meanMaxRootDispersalRates = null;
                double[][] hpdMaxRootDispersalRates = null;

                double[] meanDiffusionCoefficients = null;
                double[][] hpdDiffusionCoefficients = null;

                if (sdr) {
                    meanMaxPathDispersalDistances = meanColNoNaN(sliceTreeMaxPathDistances);
                    hpdMaxPathDispersalDistances = getArrayHPDintervals(sliceTreeMaxPathDistances);
                    meanMaxPathDispersalRates = meanColNoNaN(sliceTreeMaxPathRates);
                    hpdMaxPathDispersalRates = getArrayHPDintervals(sliceTreeMaxPathRates);

                    meanMaxRootDispersalDistances = meanColNoNaN(sliceTreeMaxDistancesFromRoot);
                    hpdMaxRootDispersalDistances = getArrayHPDintervals(sliceTreeMaxDistancesFromRoot);
                    meanMaxRootDispersalRates = meanColNoNaN(sliceTreeMaxDistanceFromRootRates);
                    hpdMaxRootDispersalRates = getArrayHPDintervals(sliceTreeMaxDistanceFromRootRates);

                    meanDiffusionCoefficients = meanColNoNaN(sliceTreeDiffusionCoefficients);
                    hpdDiffusionCoefficients = getArrayHPDintervals(sliceTreeDiffusionCoefficients);
                }

                //double[] meanWeightedAverageDiffusionCoefficients = meanColNoNaN(sliceTreeWeightedAverageDiffusionCoefficients);
                //double[][] hpdWeightedAverageDiffusionCoefficients = getArrayHPDintervals(sliceTreeWeightedAverageDiffusionCoefficients);
                for (int u = 0; u < sliceCount; u++) {
                    sliceDispersalRateFile.print(sliceHeights[u] + "\t");
                    if (mostRecentSamplingDate > 0) {
                        sliceDispersalRateFile.print((mostRecentSamplingDate - sliceHeights[u]) + "\t");
                    }
                    sliceDispersalRateFile.print(meanWeightedAverageDispersalRates[u] + "\t" + hpdWeightedAverageDispersalRates[u][0] + "\t" + hpdWeightedAverageDispersalRates[u][1]);

                    if (sdr){
                           sliceDispersalRateFile.print("\t" + meanMaxRootDispersalRates[u] + "\t" + hpdMaxRootDispersalRates[u][0] + "\t" + hpdMaxRootDispersalRates[u][1] + "\t" + meanMaxRootDispersalDistances[u] + "\t" + hpdMaxRootDispersalDistances[u][0] + "\t" + hpdMaxRootDispersalDistances[u][1] + "\t" + meanDiffusionCoefficients[u] + "\t" + hpdDiffusionCoefficients[u][0] + "\t" + hpdDiffusionCoefficients[u][1] + "\r");
//                        sliceDispersalRateFile.print("\t" + meanMaxPathDispersalRates[u] + "\t" + hpdMaxPathDispersalRates[u][0] + "\t" + hpdMaxPathDispersalRates[u][1] + "\t" + meanMaxPathDispersalDistances[u] + "\t" + hpdMaxPathDispersalDistances[u][0] + "\t" + hpdMaxPathDispersalDistances[u][1] + "\t" + meanDiffusionCoefficients[u] + "\t" + hpdDiffusionCoefficients[u][0] + "\t" + hpdDiffusionCoefficients[u][1] + "\r");
                    }  else if (snr){
                           sliceDispersalRateFile.print("\r");
                    }

                }
                sliceDispersalRateFile.close();
            } catch (IOException e) {
                System.err.println("IO Exception encountered: " + e.getMessage());
                System.exit(-1);
            }
        }
    }

    enum Normalization {
        LENGTH,
        HEIGHT,
        NONE
    }

    enum OutputFormat {
        TAB,
        KML,
        XML
    }

    enum BranchSet {
        ALL,
        INT,
        EXT,
        BACKBONE,
        CLADE
    }

    enum SliceMode {
        BRANCHES,
        NODES,
    }

    public static <T extends Enum<T>> String[] enumNamesToStringArray(T[] values) {
        int i = 0;
        String[] result = new String[values.length];
        for (T value : values) {
            result[i++] = value.name();
        }
        return result;
    }


    private void addDimInfo(Element element, int j, int dim) {
        if (dim > 1)
            element.setAttribute("dim", Integer.toString(j + 1));
    }

    private void summarizeRoot(boolean contours, boolean points, OutputFormat outputFormat, double hpdValue){

        if (sliceProgressReport) {
            progressStream.print("summarizing root" + "\t");
            progressStream.print("hpd " + (hpdValue * 100) + "\t");
        }

        Element contourElement = null;
        Element pointsElement = null;
        if (contours) {
            contourElement = new Element("Folder");
            Element name =  new Element("name");
            name.addContent("root_hpd" + (hpdValue * 100));
            contourElement.addContent(name);
        }

        if (points) {
            pointsElement = new Element("Folder");
            Element name =  new Element("name");
            name.addContent("root_points");
            pointsElement.addContent(name);
        }

        for (int traitIndex = 0; traitIndex < rootValues.size(); traitIndex++) {

            List<Trait> thisTrait = rootValues.get(traitIndex);
            if (thisTrait.size() == 0) {
                return;
            }
            boolean isNumber = thisTrait.get(0).isNumber();
            boolean isMultivariate = thisTrait.get(0).isMultivariate();
            int dim = thisTrait.get(0).getDim();
            boolean isBivariate = isMultivariate && dim == 2;

            if (isNumber) {
                if (isBivariate) {

                    if (outputFormat == OutputFormat.KML) {

                        int count = thisTrait.size();
                        double[][] y = new double[dim][count];
                        double[] h = new double[count];
                        for (int i = 0; i < count; i++) {
                            Trait trait = thisTrait.get(i);
                            double[] value = trait.getValue();
                            h[i] = trait.getHeight();
                            for (int j = 0; j < dim; j++) {
                                y[j][i] = value[j];
                            }
                        }

                        if(sliceMode == SliceMode.NODES) {
                            Element nodeSliceElement = generateNodeSliceElement(Double.NaN, y, -1);
                            nodeFolderElement.addContent(nodeSliceElement);
                        }

                        if (contourElement != null) {
                            String name = "root_hpd" + (hpdValue * 100);
                            generateContours(name, contourElement, null, y, -1, Double.NaN, Double.NaN, hpdValue);
                        }

                        if (pointsElement != null) {
                            String name = "root_points";
                            generatePoints(name, pointsElement, y, 0, 0, h) ;
                        }

                    }
                }

            }
        }

        if (contourFolderElement != null) {
            contourFolderElement.addContent(contourElement);
        }

        if (pointsFolderElement != null) {
            pointsFolderElement.addContent(pointsElement);
        }

        if (sliceProgressReport) {
            progressStream.print("\r");
        }

    }

    private void summarizeTips(boolean contours, boolean points, OutputFormat outputFormat, double hpdValue){

        if (sliceProgressReport) {
            progressStream.print("summarizing tips" + "\t");
            progressStream.print("hpd " + (hpdValue * 100) + "\t");
        }

        Element contourElement = null;
        Element pointsElement = null;
        if (contours) {
            contourElement = new Element("Folder");
            Element name =  new Element("name");
            name.addContent("tips_hpd" + (hpdValue * 100));
            contourElement.addContent(name);
        }

        if (points) {
            pointsElement = new Element("Folder");
            Element name =  new Element("name");
            name.addContent("tips_points");
            pointsElement.addContent(name);
        }

        for (int traitIndex = 0; traitIndex < tipValues.size(); traitIndex++) {

            List<List<Trait>> thisTrait = tipValues.get(traitIndex);
            for (int tipIndex = 0; tipIndex < thisTrait.size(); tipIndex++) {
                List<Trait> thisTip = thisTrait.get(tipIndex);
                if (thisTip.size() == 0) {
                    return;
                }
                boolean isNumber = thisTip.get(0).isNumber();
                boolean isMultivariate = thisTip.get(0).isMultivariate();
                int dim = thisTip.get(0).getDim();
                boolean isBivariate = isMultivariate && dim == 2;

                if (isNumber) {
                    if (isBivariate) {

                        if (outputFormat == OutputFormat.KML) {

                            int count = thisTrait.size();
                            double[][] y = new double[dim][count];
                            double[] h = new double[count];

                            for (int i = 0; i < count; i++) {
                                Trait trait = thisTip.get(i);
                                h[i] = trait.getHeight();
                                double[] value = trait.getValue();
                                for (int j = 0; j < dim; j++) {
                                    y[j][i] = value[j];
                                }
                            }


                            if (contourElement != null) {
                                String name = tipNames.get(tipIndex) + "_hpd";
                                generateContours(name, contourElement, null, y, -1, Double.NaN, Double.NaN, hpdValue);
                            }

                            if (pointsElement != null) {
                                String name = tipNames.get(tipIndex) + "_points";
                                generatePoints(name, pointsElement, y, 0, 0, h) ;
                            }

                        }
                    }

                }

            }

        }

        if (contourFolderElement != null) {
            contourFolderElement.addContent(contourElement);
        }

        if (pointsFolderElement != null) {
            pointsFolderElement.addContent(pointsElement);
        }

        if (sliceProgressReport) {
            progressStream.print("\r");
        }

    }

    private void summarizeSlice(int slice, double sliceValue, boolean contours, boolean points, OutputFormat outputFormat, double hpdValue) {

        //if (outputFormat == OutputFormat.TAB)
        //    throw new RuntimeException("Only XML/KML output is implemented");

        Element contourElement = null;
        Element pointsElement = null;

        if (outputFormat == OutputFormat.XML) {
            if (contours) {
                contourElement = new Element(SLICE_ELEMENT);
                contourElement.setAttribute(SLICE_VALUE, Double.toString(sliceValue));
            }
        } else if (outputFormat == OutputFormat.KML) {
            if (contours) {
                contourElement = new Element("Folder");
                Element name =  new Element("name");
                name.addContent("slice" + Double.toString(sliceValue) + "_hpd" + (hpdValue * 100));
                contourElement.addContent(name);
            }

            if (points) {
                pointsElement = new Element("Folder");
                Element name =  new Element("name");
                name.addContent("slice" + Double.toString(sliceValue) + "_points");
                pointsElement.addContent(name);
            }
        }

        List<List<Trait>> thisSlice = values.get(slice);
        int traitCount = thisSlice.size();

        for (int traitIndex = 0; traitIndex < traitCount; traitIndex++) {

//            if (outputFormat == OutputFormat.KML) {
//                summarizeSliceTrait(folder, slice, thisSlice.get(traitIndex), traitIndex, sliceValue,
//                        outputFormat,
//                        hpdValue);
//
//            } else {
            summarizeSliceTrait(contourElement, pointsElement, slice, thisSlice.get(traitIndex), traitIndex, sliceValue,
                    outputFormat,
                    hpdValue);

//            }
        }

        if (outputFormat == OutputFormat.KML) {
            if (contourFolderElement != null) {
                contourFolderElement.addContent(contourElement);
            }

            if (pointsFolderElement != null) {
                pointsFolderElement.addContent(pointsElement);
            }
        } else if (outputFormat == OutputFormat.XML && contourElement != null) {
            if (contourElement != null) {
                rootElement.addContent(contourElement);
            }

            if (pointsElement != null) {
                rootElement.addContent(pointsElement);
            }
        }
    }

    private void summarizeSliceTrait(Element contourElement, Element pointsElement, int slice, List<Trait> thisTrait, int traitIndex, double sliceValue,
                                     OutputFormat outputFormat,
                                     double hpdValue) {

        if (thisTrait.size() == 0) {
            return;
        }

        boolean isNumber = thisTrait.get(0).isNumber();
        boolean isMultivariate = thisTrait.get(0).isMultivariate();
        int dim = thisTrait.get(0).getDim();
        boolean isBivariate = isMultivariate && dim == 2;
        if (sliceProgressReport) {
            progressStream.print("slice " + sliceValue + "\t");
            progressStream.print("hpd " + (hpdValue * 100) + "\t");
            if (mostRecentSamplingDate > 0) {
                progressStream.print("time=" + (mostRecentSamplingDate - sliceValue) + "\t");
            }
            progressStream.print("trait=" + traits[traitIndex] + "\t");
        }
        if (isNumber) {

            Element traitElement = null;
            if (outputFormat == OutputFormat.XML || outputFormat == OutputFormat.TAB) {
                traitElement = new Element(TRAIT);
                traitElement.setAttribute(NAME, traits[traitIndex]);
            }

            if (outputFormat == OutputFormat.KML) {
                if (useStyles) {
                    Element styleElement = new Element(STYLE);
                    constructPolygonStyleElement(styleElement, sliceValue);
                    documentElement.addContent(styleElement);
                }
            }

            int count = thisTrait.size();
//            System.out.println("count = "+count+", dim = "+dim);
            double[][] y = new double[dim][count];
            for (int i = 0; i < count; i++) {
                Trait trait = thisTrait.get(i);
                double[] value = trait.getValue();
                for (int j = 0; j < dim; j++) {
                    y[j][i] = value[j];
                }
            }
//            System.out.println(y.length+"\t"+y[0].length);

            if (outputFormat == OutputFormat.XML || outputFormat == OutputFormat.TAB) {
                // Compute marginal means and standard deviations
                for (int j = 0; j < dim; j++) {
                    List<Double> x = new ArrayList();
                    for (int k = 0; k < y[j].length; k++) {
                        x.add(y[j][k]);
                    }
                    TraceDistribution trace = new TraceDistribution(x, TraceType.REAL);
                    Element statsElement = new Element("stats");
                    addDimInfo(statsElement, j, dim);
                    StringBuffer sb = new StringBuffer();
                    sb.append(KMLCoordinates.NEWLINE);
                    tabOutput.append(KMLCoordinates.NEWLINE);
                    tabOutput.append(traits[traitIndex] + "\t");
                    if (mostRecentSamplingDate > 0) {
                        tabOutput.append((mostRecentSamplingDate - sliceValue) + "\t");
                    } else {
                        tabOutput.append(sliceValue + "\t");
                    }
                    sb.append(String.format(KMLCoordinates.FORMAT,
                            trace.getMean())).append(KMLCoordinates.SEPARATOR);
                    tabOutput.append(String.format(KMLCoordinates.FORMAT,
                            trace.getMean())).append("\t");
                    sb.append(String.format(KMLCoordinates.FORMAT,
                            trace.getStdError())).append(KMLCoordinates.SEPARATOR);
                    tabOutput.append(String.format(KMLCoordinates.FORMAT,
                            trace.getStdError())).append("\t");
                    sb.append(String.format(KMLCoordinates.FORMAT,
                            trace.getLowerHPD())).append(KMLCoordinates.SEPARATOR);
                    tabOutput.append(String.format(KMLCoordinates.FORMAT,
                            trace.getLowerHPD())).append("\t");
                    sb.append(String.format(KMLCoordinates.FORMAT,
                            trace.getUpperHPD())).append(KMLCoordinates.NEWLINE);
                    tabOutput.append(String.format(KMLCoordinates.FORMAT,
                            trace.getUpperHPD())).append("\t");
                    statsElement.addContent(sb.toString());
                    traitElement.addContent(statsElement);
                }
            }

            if (isBivariate) {

                Element nodeSliceElement = null;

                if(sliceMode == SliceMode.NODES) {
                    nodeSliceElement = generateNodeSliceElement(sliceValue, y, slice);
                    nodeFolderElement.addContent(nodeSliceElement);

                    if (useStyles) {
                        Element styleElement = new Element(STYLE);
                        constructNodeStyleElement(styleElement, sliceValue);
                        documentElement.addContent(styleElement);
                    }
                }

                double date = mostRecentSamplingDate - sliceValue;

                if (pointsElement != null) {
                    String name = "" + date + "_points";
                    generatePointsElement(name, pointsElement, y, date, sliceValue);
                }

                if (contourElement != null) {
                    String name = "" + date + "_hpd" + hpdValue;
                    generateContours(name, contourElement, traitElement, y, slice, date, sliceValue, hpdValue);
                }

            }
            if (outputFormat == OutputFormat.XML)
                contourElement.addContent(traitElement);
        } // else skip
        if (sliceProgressReport) {
            progressStream.print("\r");
        }
    }

    private void generatePointsElement(String name, Element pointsFolderElement, double[][] y, double date, double height) {
        Element pointsElement = new Element("Folder");
        Element nameElement = new Element("name");
        nameElement.addContent(name);
        pointsElement.addContent(nameElement);

        for (int a = 0; a < y[0].length; a++)  {
            Element placemarkElement = new Element("Placemark");

            if (sliceCount > 1) {
                Element timeSpan = new Element("TimeSpan");
                Element begin = new Element("begin");

                if (!ancient) {
                    calendar.set(Calendar.YEAR, (int) Math.floor(date));
                    calendar.set(Calendar.DAY_OF_YEAR, (int) (365 * (date - Math.floor(date))));

                    begin.addContent(dateFormat.format(calendar.getTime()));
                } else {
                    begin.addContent(Integer.toString((int)Math.round(date)));
                }
                timeSpan.addContent(begin);
                placemarkElement.addContent(timeSpan);
            }

            placemarkElement.addContent(generatePointData(name, height));

            Element pointElement = new Element("Point");
            Element coordinates = new Element("coordinates");
            coordinates.addContent(y[1][a]+","+y[0][a]+",0");
            pointElement.addContent(coordinates);
            placemarkElement.addContent(pointElement);

            pointsElement.addContent(placemarkElement);
        }

        pointsFolderElement.addContent(pointsElement);
    }

    private void generatePoints(String name, Element pointsFolderElement, double[][] y, double date, double height, double[] heights) {
        for (int a = 0; a < y[0].length; a++)  {
            Element placemarkElement = new Element("Placemark");

            if (heights == null) {
                placemarkElement.addContent(generatePointData(name, height));
            } else {
                placemarkElement.addContent(generatePointData(name, heights[a]));
            }

            Element pointElement = new Element("Point");
            Element coordinates = new Element("coordinates");
            coordinates.addContent(y[1][a]+","+y[0][a]+",0");
            pointElement.addContent(coordinates);
            placemarkElement.addContent(pointElement);

            pointsFolderElement.addContent(placemarkElement);
        }
    }

    private void generateContours(String name, Element sliceElement, Element traitElement, double[][] y, int slice, double date, double height, double hpdValue) {
        //to test how much points are within the polygons
        double numberOfPointsInPolygons = 0;
        double totalArea = 0;

        ContourMaker contourMaker;
        if (contourMode == ContourMode.JAVA)
//            contourMaker = new KernelDensityEstimator2D(y[0], y[1], gridSize);
            contourMaker = new KernelDensityEstimator2D(y[0], y[1], BANDWIDTHLIMIT);
        else if (contourMode == ContourMode.R)
            contourMaker = new ContourWithR(y[0], y[1], gridSize);
        else if (contourMode == ContourMode.SNYDER)
//            contourMaker = new ContourWithSynder(y[0], y[1], gridSize);
            contourMaker = new ContourWithSynder(y[0], y[1], BANDWIDTHLIMIT);
        else
            throw new RuntimeException("Unimplemented ContourModel!");

        ContourPath[] paths = contourMaker.getContourPaths(hpdValue);
        int pathCounter = 1;
        for (ContourPath path : paths) {

            KMLCoordinates coords = new KMLCoordinates(path.getAllX(), path.getAllY());
            if (traitElement != null) {
                Element regionElement = new Element(REGIONS_ELEMENT);
                regionElement.setAttribute(DENSITY_VALUE, Double.toString(hpdValue));
                regionElement.addContent(coords.toXML());
                traitElement.addContent(regionElement);
            }

            // only if the trait is location we will write KML
            if (sliceElement != null) {
                //because KML polygons require long,lat,alt we need to switch lat and long first
                coords.switchXY();

                String name1 = name + "_path_"+pathCounter;

                Element placemarkElement = generatePlacemarkElementWithPolygon(name1, date, hpdValue, height, coords, slice, pathCounter);

                //testing how many points are within the polygon
                if (checkSliceContours) {
                    Element testElement = new Element("test");
                    testElement.addContent(coords.toXML());
                    Polygon2D testPolygon = new Polygon2D(testElement);
                    totalArea += testPolygon.calculateArea();

                    double[][] dest = new double[y.length][y[0].length];
                    for (int i = 0; i < y.length; i++) {
                        for (int j = 0; j < y[0].length; j++) {
                            dest[i][j] = y[i][j];
                        }
                    }
                    numberOfPointsInPolygons += getNumberOfPointsInPolygon(dest, testPolygon);
                }
                sliceElement.addContent(placemarkElement);
            }
            pathCounter ++;
        }

        if (checkSliceContours) {
            progressStream.print("numberOfContours=" + paths.length + "\tfreqOfPointsInContour=" + numberOfPointsInPolygons / y[0].length + "\ttotalArea = " + totalArea);
        }
        if (paths.length == 0 && !sliceProgressReport) {
            progressStream.println("Warning: slice at height " + height + ", contains no contours.");
        }
    }

    public static int getNumberOfPointsInPolygon(double[][] pointsArray, Polygon2D testPolygon) {
        int numberOfPointsInPolygon = 0;
        for (int x = 0; x < pointsArray[0].length; x++) {
            if (testPolygon.containsPoint2D(new Point2D.Double(pointsArray[0][x], pointsArray[1][x]))) {
                numberOfPointsInPolygon++;
            }
        }
        return numberOfPointsInPolygon;
    }

    private void constructPolygonStyleElement(Element styleElement, double sliceValue) {
        double date;
        if (Double.isNaN(sliceValue)){
            date = Double.NaN;
            styleElement.setAttribute(ID, ROOT_ELEMENT + "_hpd" + "_style");
        } else {
            date = mostRecentSamplingDate - sliceValue;
            styleElement.setAttribute(ID, REGIONS_ELEMENT + date + "_style");
        }
        Element lineStyle = new Element("LineStyle");
        Element width = new Element("width");
        width.addContent(WIDTH);
        lineStyle.addContent(width);
        Element polyStyle = new Element("PolyStyle");
        Element color = new Element("color");
        double[] minMax = new double[2];
        minMax[0] = sliceHeights[0];
        minMax[1] = sliceHeights[(sliceHeights.length - 1)];
        String colorString;
        if (Double.isNaN(sliceValue)){
            colorString = startHPDColor;
        } else {
            colorString = getKMLColor(sliceValue, minMax, endHPDColor, startHPDColor);
        }
        color.addContent(opacity + colorString);
        Element outline = new Element("outline");
        outline.addContent("0");
        polyStyle.addContent(color);
        polyStyle.addContent(outline);
        styleElement.addContent(lineStyle);
        styleElement.addContent(polyStyle);
    }

    private void constructNodeStyleElement(Element styleElement, double sliceValue) {
        double date;
        if (Double.isNaN(sliceValue)){
            date = Double.NaN;
            styleElement.setAttribute(ID, ROOT_ELEMENT + date + "_style");
        } else {
            date = mostRecentSamplingDate - sliceValue;
            styleElement.setAttribute(ID, NODE_ELEMENT + date + "_style");
        }
        Element iconStyle = new Element("IconStyle");
        Element scale = new Element("scale");
        scale.addContent("1");
        iconStyle.addContent(scale);

        Element icon = new Element("Icon");
        Element href = new Element("href");
        href.addContent(ICON);
        icon.addContent(href);
        iconStyle.addContent(icon);

        Element color = new Element("color");
        double[] minMax = new double[2];
        minMax[0] = sliceHeights[0];
        minMax[1] = sliceHeights[(sliceHeights.length - 1)];
        String colorString;
        if (Double.isNaN(sliceValue)){
            colorString = startHPDColor;
        } else {
            colorString = getKMLColor(sliceValue, minMax, endHPDColor, startHPDColor);
        }
        color.addContent(opacity + colorString);
        iconStyle.addContent(color);

        Element colorMode = new Element("colorMode");
        colorMode.addContent("normal");
        iconStyle.addContent(colorMode);

        styleElement.addContent(iconStyle);
    }

    private Element generatePlacemarkElementWithPolygon(String name, double date, double hpdValue, double sliceValue, KMLCoordinates coords, int sliceInteger, int pathCounter) {
        Element placemarkElement = new Element("Placemark");

        Element placemarkNameElement = new Element("name");
        placemarkNameElement.addContent(name);
        placemarkElement.addContent(placemarkNameElement);


        Element visibility = new Element("visibility");
        if (sliceInteger == -1){
            visibility.addContent("0");
        } else {
            visibility.addContent("1");
        }
        placemarkElement.addContent(visibility);

        if (!Double.isNaN(date)) {
            Element timeSpan = new Element("TimeSpan");
            Element begin = new Element("begin");

            if (!ancient) {
                calendar.set(Calendar.YEAR, (int) Math.floor(date));
                calendar.set(Calendar.DAY_OF_YEAR, (int) (365 * (date - Math.floor(date))));

                begin.addContent(dateFormat.format(calendar.getTime()));
            } else {
                begin.addContent(Integer.toString((int)Math.round(date)));
            }
            timeSpan.addContent(begin);
//            if (sliceInteger > 1) {
//                Element end = new Element("end");
//                end.addContent(Double.toString(mostRecentSamplingDate- sliceHeights[(sliceInteger-1)]));
//                timeSpan.addContent(end);
//            }
            placemarkElement.addContent(timeSpan);
        }

        if (useStyles) {
            Element style = new Element("styleUrl");
            if (sliceInteger == -1){
                style.addContent("#" + ROOT_ELEMENT + "_hpd" + "_style");
            }  else {
                style.addContent("#" + REGIONS_ELEMENT + date + "_style");
            }

            placemarkElement.addContent(style);
        }

        placemarkElement.addContent(generateContourData(name, date, sliceValue, hpdValue));

        Element polygonElement = new Element("Polygon");
        Element altitudeMode = new Element("altitudeMode");
        altitudeMode.addContent("clampToGround");
        polygonElement.addContent(altitudeMode);
        Element tessellate = new Element("tessellate");
        tessellate.addContent("1");
        polygonElement.addContent(tessellate);
        Element outerBoundaryIs = new Element("outerBoundaryIs");
        Element LinearRing = new Element("LinearRing");
        LinearRing.addContent(coords.toXML());
        outerBoundaryIs.addContent(LinearRing);
        polygonElement.addContent(outerBoundaryIs);
        placemarkElement.addContent(polygonElement);

        return placemarkElement;
    }

    private Element generateNodeSliceElement(double sliceValue, double[][] nodes, int sliceInteger) {
        double date;
        Element nodeSliceElement = new Element("Folder");
        Element nameNodeSliceElement = new Element("name");
        String name;
        if (sliceInteger == -1){
            date = Double.NaN;
            name = ROOT_ELEMENT;
        } else {
            date = mostRecentSamplingDate - sliceValue;
            name = "nodeSlice_"+sliceInteger+"_"+date;
        }
        nameNodeSliceElement.addContent(name);
        nodeSliceElement.addContent(nameNodeSliceElement);

        Element initialVisibility = new Element("visibility");
        initialVisibility.addContent("0");
        if (sliceInteger == -1){
            nodeSliceElement.addContent(initialVisibility);
        }

        for (int a = 0; a < nodes[0].length; a++)  {
            Element placemarkElement = new Element("Placemark");

//            Element visibility = new Element("visibility");
//            visibility.addContent("1");
//            placemarkElement.addContent(visibility);

            if (sliceCount > 1) {
                Element timeSpan = new Element("TimeSpan");
                Element begin = new Element("begin");

                if (!ancient) {
                    calendar.set(Calendar.YEAR, (int) Math.floor(date));
                    calendar.set(Calendar.DAY_OF_YEAR, (int) (365 * (date - Math.floor(date))));

                    begin.addContent(dateFormat.format(calendar.getTime()));
                } else {
                    begin.addContent(Integer.toString((int)Math.round(date)));
                }
                timeSpan.addContent(begin);
                placemarkElement.addContent(timeSpan);
            }

            placemarkElement.addContent(generatePointData(name, sliceValue));

            if (useStyles) {
                Element style = new Element("styleUrl");
                if (sliceInteger == -1){
                    style.addContent("#" + ROOT_ELEMENT + date + "_style");
                } else {
                    style.addContent("#" + NODE_ELEMENT + date + "_style");
                }

                placemarkElement.addContent(style);
            }

            Element pointElement = new Element("Point");
            //        Element altitudeMode = new Element("altitudeMode");
            //        altitudeMode.addContent("clampToGround");
            //        polygonElement.addContent(altitudeMode);
            //        Element tessellate = new Element("tessellate");
            //        tessellate.addContent("1");
            //        polygonElement.addContent(tessellate);
            Element coordinates = new Element("coordinates");
            coordinates.addContent(nodes[1][a]+","+nodes[0][a]+",0");
            pointElement.addContent(coordinates);
            placemarkElement.addContent(pointElement);

            nodeSliceElement.addContent(placemarkElement);

        }

        return nodeSliceElement;
    }

    private Element generateContourData(String name, double date, double height, double hpd) {
        Element data = new Element("ExtendedData");
        Element schemaData = new Element("SchemaData");
        schemaData.setAttribute("schemaUrl", "#HPD_Schema");
        schemaData.addContent(new Element("SimpleData").setAttribute("name", "Name").addContent(name));
//        schemaData.addContent(new Element("SimpleData").setAttribute("name", "Description"));
        schemaData.addContent(new Element("SimpleData").setAttribute("name", "Time").addContent(Double.toString(date)));
        schemaData.addContent(new Element("SimpleData").setAttribute("name", "Height").addContent(Double.toString(height)));
        if (hpd > 0) {
            schemaData.addContent(new Element("SimpleData").setAttribute("name", "HPD").addContent(Double.toString(hpd)));
        }
        data.addContent(schemaData);
        return data;
    }

    private Element generatePointData(String name, double height) {
        Element data = new Element("ExtendedData");
        Element schemaData = new Element("SchemaData");
        schemaData.setAttribute("schemaUrl", "#Point_Schema");
        double date = mostRecentSamplingDate - height;
        schemaData.addContent(new Element("SimpleData").setAttribute("name", "Name").addContent(name));
        schemaData.addContent(new Element("SimpleData").setAttribute("name", "Time").addContent(Double.toString(date)));
        schemaData.addContent(new Element("SimpleData").setAttribute("name", "Height").addContent(Double.toString(height)));
        data.addContent(schemaData);
        return data;
    }
    private void outputHeader(String[] traits) {
        StringBuffer sb = new StringBuffer("slice");
        for (int i = 0; i < traits.length; i++) {
            // Load first value to check dimensionality
            Trait trait = values.get(0).get(i).get(0);
            if (trait.isMultivariate()) {
                int dim = trait.getDim();
                for (int j = 1; j <= dim; j++)
                    sb.append(sep).append(traits[i]).append(j);
            } else
                sb.append(sep).append(traits[i]);
        }
        sb.append("\n");
        resultsStream.print(sb);
    }

//    private List<Tree> importTrees(String treeFileName, int burnin) throws IOException, Importer.ImportException {
//
//        int totalTrees = 10000;
//
//        progressStream.println("Reading trees (bar assumes 10,000 trees)...");
//        progressStream.println("0              25             50             75            100");
//        progressStream.println("|--------------|--------------|--------------|--------------|");
//
//        int stepSize = totalTrees / 60;
//        if (stepSize < 1) stepSize = 1;
//
//        List<Tree> treeList = new ArrayList<Tree>();
//        BufferedReader reader1 = new BufferedReader(new FileReader(treeFileName));
//
//        String line1 = reader1.readLine();
//        TreeImporter importer1;
//        if (line1.toUpperCase().startsWith("#NEXUS")) {
//            importer1 = new NexusImporter(new FileReader(treeFileName));
//        } else {
//            importer1 = new NewickImporter(new FileReader(treeFileName));
//        }
//        totalTrees = 0;
//        while (importer1.hasTree()) {
//            Tree treeTime = importer1.importNextTree();
//
//            if (totalTrees > burnin)
//                treeList.add(treeTime);
//
//            if (totalTrees > 0 && totalTrees % stepSize == 0) {
//                progressStream.print("*");
//                progressStream.flush();
//            }
//            totalTrees++;
//        }
//        return treeList;
//    }

    private void readAndAnalyzeTrees(String treeFileName, int burnin, int skipEvery,
                                     String[] traits, double[] slices,
                                     boolean impute, boolean trueNoise, Normalization normalize,
                                     boolean divideByBranchLength, BranchSet branchset, Set taxaSet)
            throws IOException, Importer.ImportException {

        int totalTrees = 10000;
        int totalStars = 0;

        progressStream.println("Reading and analyzing trees (bar assumes 10,000 trees)...");
        progressStream.println("0              25             50             75            100");
        progressStream.println("|--------------|--------------|--------------|--------------|");

        int stepSize = totalTrees / 60;
        if (stepSize < 1) stepSize = 1;

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
            if (totalTrees % skipEvery == 0) {
                treesRead++;
                if (totalTrees >= burnin) {
                    analyzeTree(treeTime, traits, slices, impute, trueNoise, normalize, divideByBranchLength, branchset, taxaSet);
                }
            }
            if (totalTrees > 0 && totalTrees % stepSize == 0) {
                progressStream.print("*");
                totalStars++;
                if (totalStars % 61 == 0)
                    progressStream.print("\n");
                progressStream.flush();
            }
            totalTrees++;
        }
        progressStream.print("\n");
    }

    class Trait {

        Trait(Object obj) {
            this.obj = obj;
            if (obj instanceof Object[]) {
                isMultivariate = true;
                array = (Object[]) obj;
            }
            this.height = 0.0;
        }

        Trait(Object obj, double height) {
            this.obj = obj;
            if (obj instanceof Object[]) {
                isMultivariate = true;
                array = (Object[]) obj;
            }
            this.height = height;
        }

        public boolean isMultivariate() {
            return isMultivariate;
        }

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
            if (!isMultivariate) {
                result[0] = (Double) obj;
            } else {
                for (int i = 0; i < dim; i++)
                    result[i] = (Double) array[i];
            }
            return result;
        }

        public double getHeight() {
            return height;
        }

        public void multiplyBy(double factor) {
            if (!isMultivariate) {
                obj = ((Double) obj * factor);
            } else {
                for (int i = 0; i < array.length; i++) {
                    array[i] = ((Double) array[i] * factor);
                }
            }
        }

        private Object obj;
        private Object[] array;
        private boolean isMultivariate = false;
        private double height;

        public String toString() {
            if (!isMultivariate)
                return obj.toString();
            StringBuffer sb = new StringBuffer(array[0].toString());
            for (int i = 1; i < array.length; i++)
                sb.append(sep).append(array[i]);
            return sb.toString();
        }
    }

    private List<List<List<Trait>>> values;
    private List<List<Trait>> rootValues;
    private List<List<List<Trait>>> tipValues;
    private List<String> tipNames;

    private void outputSlice(int slice, double sliceValue) {

        List<List<Trait>> thisSlice = values.get(slice);
        int traitCount = thisSlice.size();
        int valueCount = thisSlice.get(0).size();

        StringBuffer sb = new StringBuffer();

        for (int v = 0; v < valueCount; v++) {
            if (Double.isNaN(sliceValue))
                sb.append("All");
            else
                sb.append(sliceValue);
            for (int t = 0; t < traitCount; t++) {
                sb.append(sep);
                sb.append(thisSlice.get(t).get(v));
            }
            sb.append("\n");
        }

        resultsStream.print(sb);
    }

    private static boolean onBackbone(Tree tree, NodeRef node, Set targetSet) {

        if (tree.isExternal(node)) return false;

        Set leafSet = TreeUtils.getDescendantLeaves(tree, node);
        int size = leafSet.size();

        leafSet.retainAll(targetSet);

        if (leafSet.size() > 0) {

            // if all leaves below are in target then check just above.
            if (leafSet.size() == size) {

                Set superLeafSet = TreeUtils.getDescendantLeaves(tree, tree.getParent(node));
                superLeafSet.removeAll(targetSet);

                // the branch is on ancestral path if the super tree has some non-targets in it
                return (superLeafSet.size() > 0);

            } else return true;

        } else return false;
    }

    private static boolean inClade(Tree tree, NodeRef node, Set targetSet, boolean includeStem) {

        Set leafSet = TreeUtils.getDescendantLeaves(tree, node);
        int size = leafSet.size();

        leafSet.retainAll(targetSet);

        if (leafSet.size() > 0) {

            // check if node is not ancestral to mrca of targetSet
            if (leafSet.size() == size) {

                if (!includeStem){
                    return false;
                } else {
                    Set newLeafSet = TreeUtils.getDescendantLeaves(tree, node);
                    newLeafSet.removeAll(targetSet);
                    if (newLeafSet.size() == 0){
                        return true;
                    } else {
                        return false;
                    }
                }

            } else return true;

        } else return false;
    }

    private void analyzeTree(Tree treeTime, String[] traits, double[] slices, boolean impute,
                             boolean trueNoise, Normalization normalize, boolean divideByBranchlength,
                             BranchSet branchset, Set taxaSet) {

        double[][] precision = null;

        if (impute) {

            Object o = treeTime.getAttribute(PRECISION_STRING);
            double treeNormalization = 1; // None
            if (normalize == Normalization.LENGTH) {
                treeNormalization = TreeUtils.getTreeLength(treeTime, treeTime.getRoot());
            } else if (normalize == Normalization.HEIGHT) {
                treeNormalization = treeTime.getNodeHeight(treeTime.getRoot());
            }

            if (o != null) {
                Object[] array = (Object[]) o;
                int dim = (int) Math.sqrt(1 + 8 * array.length) / 2;
                precision = new double[dim][dim];
                int c = 0;
                for (int i = 0; i < dim; i++) {
                    for (int j = i; j < dim; j++) {
                        precision[j][i] = precision[i][j] = ((Double) array[c++]) * treeNormalization;
                    }
                }
            }
            //System.out.println(treeNormalization+"\t"+precision[0][0]+"\t"+precision[0][1]+"\t"+precision[1][0]+"\t"+precision[1][1]);
        }


        if (tipValues != null && tipValues.size() == 0) {
            // this is the first tree so initialize the tip value lists
            for (int i = 0; i < treeTime.getExternalNodeCount(); i++) {
                List<List<Trait>> thisTip = new ArrayList<List<Trait>>(traitCount);
                tipValues.add(thisTip);
                for (int j = 0; j < traitCount; j++) {
                    List<Trait> thisTipTrait = new ArrayList<Trait>();
                    thisTip.add(thisTipTrait);
                }

                tipNames.add(treeTime.getNodeTaxon(treeTime.getExternalNode(i)).getId());
            }

        }

//  employed to get dispersal rates across the whole tree
//        double treeNativeDistance = 0;
//        double treeKilometerGreatCircleDistance = 0;

        double[] treeSliceTime = new double[sliceCount];
        double[] treeSliceDistance = new double[sliceCount];
        double[] treeSliceMaxDistance = new double[sliceCount];
        double[] treeTimeFromRoot = new double[sliceCount];
        double[] maxDistanceFromRoot = new double[sliceCount];
        //this one will used for weighted average (WA)
        //double[] treeSliceDiffusionCoefficientWA = new double[sliceCount];
        //this one is used for simple average
        double[] treeSliceDiffusionCoefficientA = new double[sliceCount];
        // this is for the variance
        double[] treeSliceDiffusionCoefficientV = new double[sliceCount];
        double[][] treeSliceDiffusionCoefficients = new double[sliceCount][treeTime.getNodeCount() - 1];
        double[] treeSliceBranchCount = new double[sliceCount];

        treeLengths.add(TreeUtils.getTreeLength(treeTime, treeTime.getRoot()));

        for (int x = 0; x < treeTime.getNodeCount(); x++) {

            NodeRef node = treeTime.getNode(x);

            if (!(treeTime.isRoot(node))) {

                double nodeHeight = treeTime.getNodeHeight(node);
                double parentHeight = treeTime.getNodeHeight(treeTime.getParent(node));

                double oneOverBranchLength = 1.0;
                if (divideByBranchlength) {
                    oneOverBranchLength = 1.0 / treeTime.getBranchLength(node);
                }

                boolean proceed = false;
                // we allow for both branch constraints as well as constraints based on locations, so we define a temp boolean for the first type of constraints to combine with the second one
                boolean tempBranchProceed = false;
                if (branchset == BranchSet.ALL) {
                    tempBranchProceed = true;
                } else if (branchset == BranchSet.EXT && treeTime.isExternal(node)) {
                    tempBranchProceed = true;
                } else if (branchset == BranchSet.INT && !treeTime.isExternal(node)) {
                    tempBranchProceed = true;
                } else if (branchset == BranchSet.BACKBONE && onBackbone(treeTime, node, taxaSet)) {
                    tempBranchProceed = true;
                }  else if (branchset == BranchSet.CLADE && inClade(treeTime, node, taxaSet, true)) {    //TODO: make the includeStem an argument
                    tempBranchProceed = true;
                }

                if (tempBranchProceed) {
                    boolean coordinatesOK = false;
                    if ((latMin > -Double.MAX_VALUE) || (latMax < Double.MAX_VALUE) || (longMin > -Double.MAX_VALUE) || (longMax < Double.MAX_VALUE)){
                        Trait locTrait = new Trait(treeTime.getNodeAttribute(node, traits[0]));
                        double[] loc = locTrait.getValue();
//                        System.out.println("loc0 and loc1 = "+loc[0]+","+loc[1]);
                        if ((latMin < loc[0]) && (latMax > loc[0]) && (longMin < loc[1]) && (longMax > loc[1])){
                            coordinatesOK = true;
//                            Checking the constraints by looking at descendents
//                            Set leafSet = Tree.Utils.getDescendantLeaves(treeTime, node);
//                            Iterator iter = leafSet.iterator();
//                            System.out.println("descendents: ");
//                            while (iter.hasNext()) {
//                                System.out.print(iter.next()+" ");
//                             }
//                            System.out.println(";");
                        }
                    }  else {
                        coordinatesOK = true;
                    }
                    boolean descendentsOK = false;
                    if (descendentTaxaSet!=null){

                        NodeRef setNode = TreeUtils.getCommonAncestorNode(treeTime, descendentTaxaSet);

                        if (setNode==null){
                            System.err.println("no common ancestor node for taxa you have defined:");
                            Iterator iter = descendentTaxaSet.iterator();
                            while (iter.hasNext()) {
                                System.out.print(iter.next()+" ");
                            }
                            System.out.println(";");
                            System.exit(-1);
                        }

//                        Set leafSet = Tree.Utils.getDescendantLeaves(treeTime, node);
//                        Iterator iter2 = leafSet.iterator();
//                        System.out.println("leafs: ");
//                        while (iter2.hasNext()) {
//                            System.out.print(iter2.next()+" ");
//                        }
//                        System.out.println(";");


                        if (node.equals(setNode)){
                            descendentsOK = true;
//                            System.out.println("descendenst are .... OK");
                        }
                    } else {
                        descendentsOK = true;
                    }

                    if (coordinatesOK && descendentsOK){
                        proceed = true;
                    }

                }

//                if(proceed){
//                    System.out.println("proceeding");
//                }

//  employed to get dispersal rates across the whole tree
//                if (containsLocation){
//
//                    Trait nodeLocationTrait = new Trait (treeTime.getNodeAttribute(node, LOCATIONTRAIT));
//                    Trait parentNodeLocationTrait = new Trait (treeTime.getNodeAttribute(treeTime.getParent(node), LOCATIONTRAIT));
//                    treeNativeDistance += getNativeDistance(nodeLocationTrait.getValue(),parentNodeLocationTrait.getValue());
//                    treeKilometerGreatCircleDistance += getKilometerGreatCircleDistance(nodeLocationTrait.getValue(),parentNodeLocationTrait.getValue());
//
//                }

                for (int i = 0; i < sliceCount; i++) {
                    //System.out.println(slices[i]);
                    if (sdr || snr) {
                        if (!doSlices ||
                                (slices[i] < nodeHeight)
                                ) {
                            if (proceed) {

                                treeSliceTime[i] += (parentHeight - nodeHeight);
                                double diffusionCoefficient = 0;

                                if (sdr) {
                                    Trait nodeLocationTrait = new Trait(treeTime.getNodeAttribute(node, traits[0]));
                                    Trait parentNodeLocationTrait = new Trait(treeTime.getNodeAttribute(treeTime.getParent(node), traits[0]));
                                    if (GREATCIRCLEDISTANCE){
                                        treeSliceDistance[i] += getGeographicalDistance(nodeLocationTrait.getValue(), parentNodeLocationTrait.getValue());
                                        diffusionCoefficient =  (Math.pow((getGeographicalDistance(nodeLocationTrait.getValue(), parentNodeLocationTrait.getValue())), 2.0) / (4.0 * (parentHeight - nodeHeight)));
                                    }else{
                                        treeSliceDistance[i] += getNativeDistance(nodeLocationTrait.getValue(), parentNodeLocationTrait.getValue());
                                        diffusionCoefficient =  (Math.pow((getNativeDistance(nodeLocationTrait.getValue(), parentNodeLocationTrait.getValue())), 2.0) / (4.0 * (parentHeight - nodeHeight)));
                                    }

                                } else if (snr) {
                                    treeSliceDistance[i] += (Double)treeTime.getNodeAttribute(node, SUBSTITUTION);
                                }

                                //TreeModel model = new TreeModel(treeTime, true);
                                //treeSliceDistance[i] += getKilometerGreatCircleDistance(model.getMultivariateNodeTrait(node, LOCATIONTRAIT),model.getMultivariateNodeTrait(model.getParent(node), LOCATIONTRAIT));

                                //treeSliceDiffusionCoefficientWA[i] += (Math.pow((getGeographicalDistance(nodeLocationTrait.getValue(),parentNodeLocationTrait.getValue())),2.0)/(4.0*(parentHeight-nodeHeight)))*(parentHeight-nodeHeight);
                                treeSliceDiffusionCoefficientA[i] += diffusionCoefficient;
                                treeSliceDiffusionCoefficients[i][x] = diffusionCoefficient;
                                treeSliceBranchCount[i]++;
                            }
                        }
                    }

                    double height = Double.MAX_VALUE;
                    if (i < sliceCount - 1) {
                        height = slices[i + 1];
                    }

                    if (!doSlices ||
                            ((sliceMode == SliceMode.BRANCHES) && (slices[i] >= nodeHeight && slices[i] < parentHeight)) ||
                            ((sliceMode == SliceMode.NODES) && (slices[i] < nodeHeight && height >= nodeHeight))
                            ) {

                        if (proceed) {

                            List<List<Trait>> thisSlice = values.get(i);
                            for (int j = 0; j < traitCount; j++) {

                                List<Trait> thisTraitSlice = thisSlice.get(j);
                                Object tmpTrait = treeTime.getNodeAttribute(node, traits[j]);
                                if (tmpTrait == null) {
                                    System.err.println("Trait '" + traits[j] + "' not found on branch.");
                                    System.exit(-1);
                                }
                                Trait trait = new Trait(tmpTrait);
                                //System.out.println("trees "+treesAnalyzed+"\tslice "+slices[i]+"\t"+trait.toString());
                                if (divideByBranchlength) {
                                    trait.multiplyBy(oneOverBranchLength);
                                }
                                if (impute && (sliceMode == SliceMode.BRANCHES)) {
                                    double rate = 1.0;
                                    if (!rateAttributeString.equals("none")){
                                        Double rateAttribute = (Double) treeTime.getNodeAttribute(node, rateAttributeString);
                                        if (rateAttribute != null) {
                                            rate = rateAttribute;
                                            if (outputRateWarning) {
                                                progressStream.println("Warning: using "+rateAttributeString+" as rate attribute during imputation!");
                                                outputRateWarning = false;
                                            }
                                        }
                                    }
                                    if (trueNoise && precision == null) {
                                        progressStream.println("Error: not precision available for imputation with correct noise!");
                                        System.exit(-1);
                                    }
//                                    if (slices[i] > nodeHeight) {
                                    trait = imputeValue(trait, new Trait(treeTime.getNodeAttribute(treeTime.getParent(node), traits[j])),
                                            slices[i], nodeHeight, parentHeight, precision, rate, trueNoise);
//                                    System.out.println(slices[i]+"\t"+nodeHeight+"\t"+parentHeight+"\t"+precision[0][0]+"\t"+precision[0][1]+"\t"+precision[1][0]+"\t"+precision[1][1]+"\t"+rate+"\t"+trait);
////
//  }
                                    // QUESTION to PL: MAS does not see how slices[i] is ever less than nodeHeight
//                                } else if (impute && (sliceMode == SliceMode.NODES)) {
//                                    progressStream.println("no imputation for slice mode = nodes");
                                }
                                thisTraitSlice.add(trait);
                                //System.out.println("trees "+treesAnalyzed+"\tslice "+slices[i]+"\t"+trait.toString());

                                treeSliceTime[i] += (parentHeight - slices[i]);
                                double diffusionCoefficient = 0;
                               //if trait is location
                                 if (sdr) {
                                    Trait parentTrait = new Trait(treeTime.getNodeAttribute(treeTime.getParent(node), traits[j]));
                                    if (GREATCIRCLEDISTANCE){
                                        treeSliceDistance[i] += getGeographicalDistance(trait.getValue(), parentTrait.getValue());
                                        diffusionCoefficient =  (Math.pow((getGeographicalDistance(trait.getValue(), parentTrait.getValue())), 2.0) / (4.0 * (parentHeight - slices[i])));
                                    }   else{
                                        treeSliceDistance[i] += getNativeDistance(trait.getValue(), parentTrait.getValue());
                                        diffusionCoefficient =  (Math.pow((getNativeDistance(trait.getValue(), parentTrait.getValue())), 2.0) / (4.0 * (parentHeight - slices[i])));
//                                        System.out.println(getNativeDistance(trait.getValue(), parentTrait.getValue())+"\t"+(parentHeight - slices[i])+"\t"+getNativeDistance(trait.getValue(), parentTrait.getValue())/(parentHeight - slices[i]));
                                    }
                                 } else if (snr) {
                                     treeSliceDistance[i] += (Double)treeTime.getNodeAttribute(node, SUBSTITUTION)*((parentHeight - slices[i])/(parentHeight - nodeHeight));

                                 }


                                 treeSliceDiffusionCoefficientA[i] += diffusionCoefficient;
                                 treeSliceDiffusionCoefficients[i][x] = diffusionCoefficient;
                                 treeSliceBranchCount[i]++;

                                 if(sdr) {
                                    double tempDistanceFromRoot = getDistanceFromRoot(treeTime, traits[j], trait.getValue());
                                    if (maxDistanceFromRoot[i] < tempDistanceFromRoot) {
                                        maxDistanceFromRoot[i] = tempDistanceFromRoot;
                                        treeSliceMaxDistance[i] = getPathDistance(treeTime, node, traits[j], trait.getValue());
                                        //putting this below here ensures that treeTimeFromRoot is never < 0
                                        treeTimeFromRoot[i] = treeTime.getNodeHeight(treeTime.getRoot()) - slices[i];
                                    }


                                }
                            }
                        }
                    }
                }

                if (tipValues != null && treeTime.isExternal(node)) {
                    List<List<Trait>> thisTip = tipValues.get(x);

                    for (int j = 0; j < traitCount; j++) {
                        Object tmpTrait = treeTime.getNodeAttribute(node, traits[j]);
                        if (tmpTrait == null) {
                            System.err.println("Trait '" + traits[j] + "' not found for tip.");
                            System.exit(-1);
                        }
                        thisTip.get(j).add(new Trait(tmpTrait, treeTime.getNodeHeight(node)));

                    }
                }

            } else {
                if (sliceMode == SliceMode.NODES) {
                    double nodeHeight = treeTime.getNodeHeight(node);
                    for (int i = 0; i < sliceCount; i++) {
                        double height = Double.MAX_VALUE;
                        if (i < sliceCount - 1) {
                            height = slices[i + 1];
                        }
                        if ((slices[i] < nodeHeight && height >= nodeHeight)){
                            List<List<Trait>> thisSlice = values.get(i);
                            for (int j = 0; j < traitCount; j++) {
                                List<Trait> thisTraitSlice = thisSlice.get(j);
                                Object tmpTrait = treeTime.getNodeAttribute(node, traits[j]);
                                if (tmpTrait == null) {
                                    System.err.println("Trait '" + traits[j] + "' not found on node.");
                                    System.exit(-1);
                                }
                                Trait trait = new Trait(tmpTrait);
                                thisTraitSlice.add(trait);

                            }

                        }

                    }
                }

                if (rootValues != null) {
                    for (int j = 0; j < traitCount; j++) {
                        List<Trait> thisRootTrait = rootValues.get(j);
                        Object tmpTrait = treeTime.getNodeAttribute(node, traits[j]);
                        if (tmpTrait == null) {
                            System.err.println("Trait '" + traits[j] + "' not found on root node.");
                            System.exit(-1);
                        }
                        Trait trait = new Trait(tmpTrait, treeTime.getNodeHeight(node));
                        thisRootTrait.add(trait);

                    }
                }
            }
        }

        //System.out.println(Tree.Utils.getTreeLength(treeTime, treeTime.getRoot())+"\t"+test);

        if (sdr || snr) {
            sliceTreeDistanceArrays.add(treeSliceDistance);
            sliceTreeTimeArrays.add(treeSliceTime);
            if (sdr){
                sliceTreeMaxPathDistanceArrays.add(treeSliceMaxDistance);
                sliceTreeMaxDistanceFromRootArrays.add(maxDistanceFromRoot);
                sliceTreeTimeFromRootArrays.add(treeTimeFromRoot);
                for (int i = 0; i < treeSliceDiffusionCoefficientA.length; i++) {
                    //treeSliceDiffusionCoefficientWA[i] = treeSliceDiffusionCoefficientWA[i]/treeSliceTime[i];
                    treeSliceDiffusionCoefficientA[i] = treeSliceDiffusionCoefficientA[i] / treeSliceBranchCount[i];
                    for (int j = 0; j < treeSliceDiffusionCoefficients[0].length; j++) {
                        treeSliceDiffusionCoefficientV[i] += Math.pow((treeSliceDiffusionCoefficients[i][j] - treeSliceDiffusionCoefficientA[i]),2);
                    }
                    treeSliceDiffusionCoefficientV[i] = treeSliceDiffusionCoefficientV[i] / treeSliceBranchCount[i];
                    //System.out.println(treeSliceTime[i]+"\t"+treeLengths.get(i));
                }
                sliceTreeDiffusionCoefficientArrays.add(treeSliceDiffusionCoefficientA);
                sliceTreeDiffusionCoefficientVarianceArrays.add(treeSliceDiffusionCoefficientV);
            }
        }


//  employed to get dispersal rates across the whole tree
//        if (containsLocation) {
//           double treelength = Tree.Utils.getTreeLength(treeTime, treeTime.getRoot());
//            double dispersalNativeRate = treeNativeDistance/treelength;
//            double dispersalKilometerRate = treeKilometerGreatCircleDistance/treelength;
        //System.out.println(dispersalNativeRate+"\t"+dispersalKilometerRate);
//            dispersalrates.add(dispersalNativeRate+"\t"+dispersalKilometerRate);
//        }

        treesAnalyzed++;

    }

    private static double getNativeDistance(double[] location1, double[] location2) {
        return Math.sqrt(Math.pow((location2[0]-location1[0]),2.0)+Math.pow((location2[1]-location1[1]),2.0));
    }


    private static double getGeographicalDistance(double[] location1, double[] location2) {
        if (location1.length == 1) {
            // assume we only have latitude so put them on the prime meridian
            return getKilometerGreatCircleDistance(new double[] { location1[0], 0.0}, new double[] { location2[0], 0.0});
        } else if (location1.length == 2) {
            return getKilometerGreatCircleDistance(location1, location2);
        }
        throw new RuntimeException("Distances can only be calculated for longitude and latitude (or just latitude)");
    }

    private static double getKilometerGreatCircleDistance(double[] location1, double[] location2) {
        SphericalPolarCoordinates coord1 = new SphericalPolarCoordinates(location1[0], location1[1]);
        SphericalPolarCoordinates coord2 = new SphericalPolarCoordinates(location2[0], location2[1]);
        return (coord1.distance(coord2));
    }

    private double getPathDistance(Tree tree, NodeRef node, String locationTrait, double[] sliceTrait) {

        double pathDistance = 0;
        NodeRef parentNode = tree.getParent(node);
        Trait parentTrait = new Trait(tree.getNodeAttribute(parentNode, locationTrait));
        if (GREATCIRCLEDISTANCE){
            pathDistance += getGeographicalDistance(sliceTrait, parentTrait.getValue());
        }  else {
            pathDistance += getNativeDistance(sliceTrait, parentTrait.getValue());
        }

        while (parentNode != tree.getRoot()) {
            node = tree.getParent(node);
            parentNode = tree.getParent(parentNode);
            Trait nodeTrait = new Trait(tree.getNodeAttribute(node, locationTrait));
            parentTrait = new Trait(tree.getNodeAttribute(parentNode, locationTrait));
            if (GREATCIRCLEDISTANCE){
                pathDistance += getGeographicalDistance(nodeTrait.getValue(), parentTrait.getValue());
            }  else {
                pathDistance += getNativeDistance(nodeTrait.getValue(), parentTrait.getValue());
            }
        }
        return pathDistance;
    }

    private double getDistanceFromRoot(Tree tree, String locationTrait, double[] sliceTrait) {
        NodeRef rootNode = tree.getRoot();
        Trait rootTrait = new Trait(tree.getNodeAttribute(rootNode, locationTrait));
        if (GREATCIRCLEDISTANCE){
            return getGeographicalDistance(sliceTrait, rootTrait.getValue());
        } else {
            return getNativeDistance(sliceTrait, rootTrait.getValue());
        }

    }

    private int traitCount;
    private int sliceCount;
    private String[] traits;
    private double[] sliceHeights;
    private boolean sliceProgressReport;
    private boolean checkSliceContours;
    private boolean doSlices;
    private int treesRead = 0;
    private int treesAnalyzed = 0;
    private double mostRecentSamplingDate;
    private ContourMode contourMode;
    private SliceMode sliceMode;
    private boolean ancient = false;
    private boolean useStyles = true;
    private int gridSize;
    private double latMin;
    private double latMax;
    private double longMin;
    private double longMax;
    private Set descendentTaxaSet;
    private String rateAttributeString;


//  employed to get dispersal rates across the whole tree
//    private static boolean containsLocation = false;
//    private static ArrayList dispersalrates = new ArrayList();

//    private void run(List<Tree> trees, String[] traits, double[] slices, boolean impute, boolean trueNoise) {
//
//        for (Tree treeTime : trees) {
//            analyzeTree(treeTime, traits, slices, impute, trueNoise);
//        }
//    }

    private ArrayList sliceTreeDistanceArrays = new ArrayList();
    private ArrayList sliceTreeTimeArrays = new ArrayList();
    private ArrayList sliceTreeMaxPathDistanceArrays = new ArrayList();
    private ArrayList sliceTreeMaxDistanceFromRootArrays = new ArrayList();
    private ArrayList sliceTreeTimeFromRootArrays = new ArrayList();
    private ArrayList sliceTreeDiffusionCoefficientArrays = new ArrayList();
    private ArrayList sliceTreeDiffusionCoefficientVarianceArrays = new ArrayList();
    private boolean sdr;
    private boolean snr;
    private ArrayList treeLengths = new ArrayList();

    private boolean outputRateWarning = true;


    private Trait imputeValue(Trait nodeTrait, Trait parentTrait, double time, double nodeHeight, double parentHeight, double[][] precision, double rate, boolean trueNoise) {
        if (!nodeTrait.isNumber()) {
            System.err.println("Can only impute numbers!");
            System.exit(-1);
        }

        int dim = nodeTrait.getDim();
        double[] nodeValue = nodeTrait.getValue();
        double[] parentValue = parentTrait.getValue();

        final double scaledTimeChild = (time - nodeHeight) * rate;
        final double scaledTimeParent = (parentHeight - time) * rate;
        final double scaledWeightTotal = 1.0 / scaledTimeChild + 1.0 / scaledTimeParent;

        if (scaledTimeChild == 0)
            return nodeTrait;

        if (scaledTimeParent == 0)
            return parentTrait;

        // Find mean value, weighted average
        double[] mean = new double[dim];
        double[][] scaledPrecision = new double[dim][dim];

        for (int i = 0; i < dim; i++) {
            mean[i] = (nodeValue[i] / scaledTimeChild + parentValue[i] / scaledTimeParent) / scaledWeightTotal;
            if (trueNoise) {
                for (int j = i; j < dim; j++)
                    scaledPrecision[j][i] = scaledPrecision[i][j] = precision[i][j] * scaledWeightTotal;
            }
        }

//        System.out.print(time+"\t"+nodeHeight+"\t"+parentHeight+"\t"+scaledTimeChild+"\t"+scaledTimeParent+"\t"+scaledWeightTotal+"\t"+mean[0]+"\t"+mean[1]+"\t"+scaledPrecision[0][0]+"\t"+scaledPrecision[0][1]+"\t"+scaledPrecision[1][0]+"\t"+scaledPrecision[1][1]);

        if (trueNoise) {
            mean = MultivariateNormalDistribution.nextMultivariateNormalPrecision(mean, scaledPrecision);
        }
//        System.out.println("\t"+mean[0]+"\t"+mean[1]+"\r");

        Object[] result = new Object[dim];
        for (int i = 0; i < dim; i++)
            result[i] = mean[i];
        return new Trait(result);
    }

    // for debugging purposes
//    private Trait imputeValue(Trait nodeTrait, Trait parentTrait, double time, double nodeHeight, double parentHeight, double[][] precision, double rate, boolean trueNoise) {
//        Object[] result = new Object[2];
//        double[] nodeValue = nodeTrait.getValue();
//        double[] parentValue = parentTrait.getValue();
//        result[0] = (nodeValue[0] + parentValue[0])/2;
//        result[1] = (nodeValue[1] + parentValue[1])/2;
//        return new Trait(result);
//    }

    // Messages to stderr, output to stdout
    private static PrintStream progressStream = System.err;
    private PrintStream resultsStream;

    private final static Version version = new BeastVersion();

    private static final String commandName = "timeslicer";


    public static void printUsage(Arguments arguments) {

        arguments.printUsage(commandName, "<input-file-name> [<output-file-name>]");
        progressStream.println();
        progressStream.println("  Example: " + commandName + " test.trees out.kml");
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
        centreLine("Katholieke Universiteit Leuven", 60);
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

    public static double[] parseVariableLengthDoubleArray(String inString) throws Arguments.ArgumentException {

        List<Double> returnList = new ArrayList<Double>();
        StringTokenizer st = new StringTokenizer(inString, ",");
        while (st.hasMoreTokens()) {
            try {
                returnList.add(Double.parseDouble(st.nextToken()));
            } catch (NumberFormatException e) {
                throw new Arguments.ArgumentException();
            }

        }

        if (returnList.size() > 0) {
            double[] doubleArray = new double[returnList.size()];
            for (int i = 0; i < doubleArray.length; i++)
                doubleArray[i] = returnList.get(i);

            return doubleArray;
        }
        return null;
    }

    private static String[] parseVariableLengthStringArray(String inString) {

        List<String> returnList = new ArrayList<String>();
        StringTokenizer st = new StringTokenizer(inString, ",");
        while (st.hasMoreTokens()) {
            returnList.add(st.nextToken());
        }

        if (returnList.size() > 0) {
            String[] stringArray = new String[returnList.size()];
            stringArray = returnList.toArray(stringArray);
            return stringArray;
        }
        return null;
    }

    private static Set parseVariableLengthStringSet(String inString) {

        Set targetSet = new HashSet();

        StringTokenizer st = new StringTokenizer(inString, ",");
//        System.out.println(inString);
        while (st.hasMoreTokens()) {
            targetSet.add(st.nextToken());
        }

        return targetSet;
    }

    private static double[] parseFileWithArray(String file) {
        List<Double> returnList = new ArrayList<Double>();
        try {
            BufferedReader readerTimes = new BufferedReader(new FileReader(file));
            String line = readerTimes.readLine();
            while (line != null && !line.equals("")) {
                returnList.add(Double.valueOf(line));
                line = readerTimes.readLine();
            }
        } catch (IOException e) {
            System.err.println("Error reading " + file);
            System.exit(1);
        }

        if (returnList.size() > 0) {
            double[] doubleArray = new double[returnList.size()];
            for (int i = 0; i < doubleArray.length; i++)
                doubleArray[i] = returnList.get(i);

            //System.out.println(doubleArray.length);
            return doubleArray;
        }
        return null;
    }

    public static String getKMLColor(double value, double[] minMaxMedian, String startColor, String endColor) {

        startColor = startColor.toLowerCase();
        String startBlue = startColor.substring(0, 2);
        String startGreen = startColor.substring(2, 4);
        String startRed = startColor.substring(4, 6);

        endColor = endColor.toLowerCase();
        String endBlue = endColor.substring(0, 2);
        String endGreen = endColor.substring(2, 4);
        String endRed = endColor.substring(4, 6);

        double proportion = (value - minMaxMedian[0]) / (minMaxMedian[1] - minMaxMedian[0]);

        // generate an array with hexadecimal code for each RGB entry number
        String[] colorTable = new String[256];

        int colorTableCounter = 0;

        for (int a = 0; a < 10; a++) {

            for (int b = 0; b < 10; b++) {

                colorTable[colorTableCounter] = a + "" + b;
                colorTableCounter++;
            }

            for (int c = (int) ('a'); c < 6 + (int) ('a'); c++) {
                colorTable[colorTableCounter] = a + "" + (char) c;
                colorTableCounter++;
            }

        }
        for (int d = (int) ('a'); d < 6 + (int) ('a'); d++) {

            for (int e = 0; e < 10; e++) {

                colorTable[colorTableCounter] = (char) d + "" + e;
                colorTableCounter++;
            }

            for (int f = (int) ('a'); f < 6 + (int) ('a'); f++) {
                colorTable[colorTableCounter] = (char) d + "" + (char) f;
                colorTableCounter++;
            }

        }


        int startBlueInt = 0;
        int startGreenInt = 0;
        int startRedInt = 0;

        int endBlueInt = 0;
        int endGreenInt = 0;
        int endRedInt = 0;

        for (int i = 0; i < colorTable.length; i++) {

            if (colorTable[i].equals(startBlue)) {
                startBlueInt = i;
            }
            if (colorTable[i].equals(startGreen)) {
                startGreenInt = i;
            }
            if (colorTable[i].equals(startRed)) {
                startRedInt = i;
            }
            if (colorTable[i].equals(endBlue)) {
                endBlueInt = i;
            }
            if (colorTable[i].equals(endGreen)) {
                endGreenInt = i;
            }
            if (colorTable[i].equals(endRed)) {
                endRedInt = i;
            }

        }

        int blueInt = startBlueInt + (int) Math.round((endBlueInt - startBlueInt) * proportion);
        int greenInt = startGreenInt + (int) Math.round((endGreenInt - startGreenInt) * proportion);
        int redInt = startRedInt + (int) Math.round((endRedInt - startRedInt) * proportion);

        String blue = null;
        String green = null;
        String red = null;

        for (int j = 0; j < colorTable.length; j++) {

            if (j == blueInt) {
                blue = colorTable[j];
            }
            if (j == greenInt) {
                green = colorTable[j];
            }
            if (j == redInt) {
                red = colorTable[j];
            }

        }

        return blue + green + red;
    }

    private static double[] getHPDInterval(double proportion, double[] array, int[] indices) {

        double returnArray[] = new double[2];
        double minRange = Double.MAX_VALUE;
        int hpdIndex = 0;

        int diff = (int) Math.round(proportion * (double) array.length);
        for (int i = 0; i <= (array.length - diff); i++) {
            double minValue = array[indices[i]];
            double maxValue = array[indices[i + diff - 1]];
            double range = Math.abs(maxValue - minValue);
            if (range < minRange) {
                minRange = range;
                hpdIndex = i;
            }
        }
        returnArray[0] = array[indices[hpdIndex]];
        returnArray[1] = array[indices[hpdIndex + diff - 1]];
        return returnArray;
    }

    private static void print2DArray(double[][] array, String name) {
        try {
            PrintWriter outFile = new PrintWriter(new FileWriter(name), true);

            for (double[] anArray : array) {
                for (int j = 0; j < array[0].length; j++) {
                    outFile.print(anArray[j] + "\t");
                }
                outFile.println("");
            }
            outFile.close();

        } catch (IOException io) {
            System.err.print("Error writing to file: " + name);
        }
    }

    private static void print2DTransposedArray(double[][] array, String name) {
        try {
            PrintWriter outFile = new PrintWriter(new FileWriter(name), true);

            for (int i = 0; i < array[0].length; i++) {
                for (double[] anArray : array) {
                    outFile.print(anArray[i] + "\t");
                }
                outFile.println("");
            }
            outFile.close();

        } catch (IOException io) {
            System.err.print("Error writing to file: " + name);
        }
    }

    private static double[] meanColNoNaN(double[][] x) {
        double[] returnArray = new double[x[0].length];

        for (int i = 0; i < x[0].length; i++) {

            double m = 0.0;
            int lenNoZero = 0;


            for (int j = 0; j < x.length; j++) {

                if (!(((Double) x[j][i]).isNaN())) {
                    m += x[j][i];
                    lenNoZero += 1;
                }
            }
            returnArray[i] = m / (double) lenNoZero;

        }
        return returnArray;
    }

    private static double[][] getArrayHPDintervals(double[][] array) {

        double[][] returnArray = new double[array[0].length][2];

        for (int col = 0; col < array[0].length; col++) {

            int counter = 0;

            for (int row = 0; row < array.length; row++) {

                if (!(((Double) array[row][col]).isNaN())) {
                    counter += 1;
                }
            }

            if (counter > 0) {
                double[] columnNoNaNArray = new double[counter];

                int index = 0;
                for (int row = 0; row < array.length; row++) {

                    if (!(((Double) array[row][col]).isNaN())) {
                        columnNoNaNArray[index] = array[row][col];
                        index += 1;
                    }
                }
                int[] indices = new int[counter];
                HeapSort.sort(columnNoNaNArray, indices);
                double hpdBinInterval[] = getHPDInterval(0.95, columnNoNaNArray, indices);

                returnArray[col][0] = hpdBinInterval[0];
                returnArray[col][1] = hpdBinInterval[1];
            } else {
                returnArray[col][0] = Double.NaN;
                returnArray[col][1] = Double.NaN;
            }

        }

        return returnArray;
    }

    private static Set getTargetSet(String x) {
        Set targetSet = null;
        targetSet = new HashSet();
        try {
            BufferedReader reader = new BufferedReader(new FileReader(x));
            try {
                String line = reader.readLine().trim();
                while (line != null && !line.equals("")) {
                    targetSet.add(line);
                    line = reader.readLine();
                    if (line != null) line = line.trim();
                }
            }
            catch (IOException io) {
                System.err.println("Error reading " + x);
            }
        }
        catch (FileNotFoundException a) {
            System.err.println("Error finding " + x);
        }
        return targetSet;
    }

    public static void main(String[] args) throws IOException {

        String inputFileName = null;
        String outputFileName = null;
        String[] traitNames = null;
        double[] sliceHeights = null;
        OutputFormat outputFormat = OutputFormat.KML;
        boolean impute = false;
        boolean trueNoise = true;
        boolean summaryOnly = true;
        ContourMode contourMode = ContourMode.SNYDER;
        Normalization normalize = Normalization.LENGTH;
        int burnin = -1;
        int skipEvery = 1;
        double mrsd = 0;
        boolean summarizeRoot = false;
        boolean summarizeTips = false;
        boolean contours = true;
        boolean points = false;
        double[] hpdValues = { 0.80 };
        String outputFileSDR = null;
        String outputFileSNR = null;
        boolean getSDR = false;
        boolean getSNR = false;
        String progress = null;
        boolean branchNormalization = false;
        BranchSet set = BranchSet.ALL;
        Set taxaSet = null;
        SliceMode sliceMode = SliceMode.BRANCHES;
        int grid = 200;
        double latMax = Double.MAX_VALUE;
        double latMin = -Double.MAX_VALUE;
        double longMax = Double.MAX_VALUE;
        double longMin = -Double.MAX_VALUE;
        String rateString = "location.rate";
        Set descendents = null;

//        if (args.length == 0) {
//          // TODO Make flash GUI
//        }

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.IntegerOption(BURNIN, "the number of states to be considered as 'burn-in' [default = 0]"),
                        new Arguments.IntegerOption(SKIP, "skip every i'th tree [default = 0]"),
                        new Arguments.StringOption(TRAIT, "trait_name", "specifies an attribute-list to use to create a density map [default = location.rate]"),
                        new Arguments.StringOption(SLICE_TIMES, "slice_times", "specifies a slice time-list [default=none]"),
                        new Arguments.StringOption(SLICE_HEIGHTS, "slice_heights", "specifies a slice height-list [default=none]"),
                        new Arguments.StringOption(SLICE_FILE_HEIGHTS, "heights_file", "specifies a file with a slice heights-list, is overwritten by command-line specification of slice heights [default=none]"),
                        new Arguments.StringOption(SLICE_FILE_TIMES, "Times_file", "specifies a file with a slice Times-list, is overwritten by command-line specification of slice times [default=none]"),
                        new Arguments.StringOption(RATE_ATTRIBUTE, "rate_attribute", "specifies the trait rate attribute string [default=location.rate]; use 'none' when no rate needs to be used (homogeneous Brownian model)"),
                        new Arguments.IntegerOption(SLICE_COUNT, "the number of time slices to use [default=0]"),
                        new Arguments.StringOption(SLICE_MODE, "Slice_mode", "specifies how to perform the slicing [default=branches]"),
                        new Arguments.StringOption(ROOT, falseTrue, false, "include a summary for the root [default=off]"),
                        new Arguments.StringOption(TIPS, falseTrue, false, "include a summary for the tips [default=off]"),
                        new Arguments.StringOption(CONTOURS, falseTrue, true, "include contours in summary [default=true]"),
                        new Arguments.StringOption(POINTS, falseTrue, false, "include all points for the summary [default=off]"),
                        new Arguments.RealOption(START_TIME, "the time of the earliest slice [default=0]"),
                        new Arguments.RealOption(MRSD, "specifies the most recent sampling data in fractional years to rescale time [default=0]"),
                        new Arguments.Option(HELP, "option to print this message"),
                        new Arguments.StringOption(NOISE, falseTrue, false,
                                "add true noise [default = true])"),
                        new Arguments.StringOption(IMPUTE, falseTrue, false,
                                "impute trait at time-slice [default = false]"),
                        new Arguments.StringOption(SUMMARY, falseTrue, false,
                                "compute summary statistics [default = true]"),
                        new Arguments.StringOption(FORMAT, enumNamesToStringArray(OutputFormat.values()), false,
                                "summary output format [default = KML]"),
                        new Arguments.StringOption(HPD, "hpd", "mass (1 - 99%) to include in HPD regions (or list) [default = 80]"),
                        new Arguments.StringOption(CONTOUR_MODE, enumNamesToStringArray(ContourMode.values()), false,
                                "contouring model [default = snyder]"),
                        new Arguments.StringOption(NORMALIZATION, enumNamesToStringArray(Normalization.values()), false,
                                "tree normalization [default = length"),
                        new Arguments.StringOption(SDR, "sliceDispersalRate", "specifies output file name for dispersal rates for each slice (from previous sliceTime[or root of the trees] up to current sliceTime"),
                        new Arguments.StringOption(SNR, "sliceNonynonymousRate", "specifies output file name for Nonynsonymous rates for each slice (from previous sliceTime[or root of the trees] up to current sliceTime"),
                        new Arguments.StringOption(PROGRESS, "progress report", "reports slice progress and checks the bivariate contour HPD regions by calculating what fraction of points the polygons for a given slice contain  [default = false]"),
                        new Arguments.StringOption(BRANCH_NORMALIZE, falseTrue, false,
                                "devide a branch trait by branch length (can be useful for 'rewards' [default = false]"),
                        new Arguments.StringOption(BRANCHSET, TimeSlicer.enumNamesToStringArray(BranchSet.values()), false,
                                "branch set [default = all]"),
                        new Arguments.StringOption(BACKBONETAXA, "Backbone taxa file", "specifies a file with taxa that define the backbone"),
                        new Arguments.RealOption(LATMAX, "specifies the maximum latitude for a child node for a branch to be included in the summary [default=MAX_VALUE]"),
                        new Arguments.RealOption(LATMIN, "specifies the minimum latitude for a child node for a branch to be included in the summary [default=MIN_VALUE]"),
                        new Arguments.RealOption(LONGMAX, "specifies the maximum longitude for a child node for a branch to be included in the summary [default=MAX_VALUE]"),
                        new Arguments.RealOption(LONGMIN, "specifies the minimum longitude for a child node for a branch to be included in the summary [default=MIN_VALUE]"),
                        new Arguments.IntegerOption(GRIDSIZE, "the grid size for contouring [default=200]"),
                        new Arguments.StringOption(DESCENDENTS, "descendent taxa", "specifies a branch based on the descendent taxa [default=all branches]")

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

        try { // Make sense of arguments

            String sliceHeightsFileString = arguments.getStringOption(SLICE_FILE_HEIGHTS);
            if (sliceHeightsFileString != null) {
                sliceHeights = parseFileWithArray(sliceHeightsFileString);
            }

            if (arguments.hasOption(MRSD)) {
                mrsd = arguments.getRealOption(MRSD);
            }

            if (arguments.hasOption(LATMAX)) {
                latMax = arguments.getRealOption(LATMAX);
            }

            if (arguments.hasOption(LATMIN)) {
                latMin = arguments.getRealOption(LATMIN);
            }

            if (arguments.hasOption(LONGMAX)) {
                longMax = arguments.getRealOption(LONGMAX);
            }

            if (arguments.hasOption(LONGMIN)) {
                longMin = arguments.getRealOption(LONGMIN);
            }

            String taxaString = arguments.getStringOption(DESCENDENTS);
            if (taxaString != null) {
                descendents = parseVariableLengthStringSet(taxaString);
            }

                String sliceTimesFileString = arguments.getStringOption(SLICE_FILE_TIMES);
            if (sliceTimesFileString != null) {
                //System.out.println(sliceTimesFileString);
                double[] sliceTimes = parseFileWithArray(sliceTimesFileString);
                sliceHeights =  new double[sliceTimes.length];
                for (int i = 0; i < sliceTimes.length; i++) {
                    if (mrsd == 0) {
                        sliceHeights[i] = sliceTimes[i];
                    } else {
                        //System.out.println((mrsd - sliceTimes[i]));
                        sliceHeights[i] = mrsd - sliceTimes[i];

                    }
                }
            }

            String sliceModeString = arguments.getStringOption(SLICE_MODE);
            if (sliceModeString != null) {
                try {
                    sliceMode = SliceMode.valueOf(sliceModeString.toUpperCase());
                } catch (IllegalArgumentException iae) {
                    System.err.println("Unrecognized slice mode: " + sliceModeString);
                }
            }

            String traitRateString = arguments.getStringOption(RATE_ATTRIBUTE);
            if (traitRateString != null) {
                rateString = traitRateString;
            }

            if (arguments.hasOption(BURNIN)) {
                burnin = arguments.getIntegerOption(BURNIN);
                System.err.println("Ignoring a burnin of " + burnin + " trees.");
            }

            if (arguments.hasOption(SKIP)) {
                skipEvery = arguments.getIntegerOption(SKIP);
                System.err.println("Skipping every " + skipEvery + " trees.");
            }
            if (skipEvery < 1) {
                skipEvery = 1;
            }

            String hpdString = arguments.getStringOption(HPD);
            if (hpdString != null) {
                hpdValues = parseVariableLengthDoubleArray(hpdString);
                if (hpdValues.length > 0) {
                    for (int i = 0; i < hpdValues.length; i++) {
                        if (hpdValues[i] < 1 || hpdValues[i] > 99) {
                            progressStream.println("HPD Region mass falls outside of 1 - 99% range.");
                            System.exit(-1);
                        }
                        hpdValues[i] = hpdValues[i] / 100.0;
                    }

                } else {
                    hpdValues = new double[] { 80.0 };
                }
            }

            String traitString = arguments.getStringOption(TRAIT);
            if (traitString != null) {
                traitNames = parseVariableLengthStringArray(traitString);
//employed to get dispersal rates across the whole tree
//                for (int y = 0; y < traitNames.length; y++) {
//                    if (traitNames[y].equals(LOCATIONTRAIT)){
//                       containsLocation =  true;
//                    }
//                }
            }

            if (traitNames == null) {
                traitNames = new String[1];
                traitNames[0] = "location.rate";
            }

            String sliceTimeString = arguments.getStringOption(SLICE_TIMES);
            if (sliceTimeString != null) {
                double[] sliceTimes = parseVariableLengthDoubleArray(sliceTimeString);
                sliceHeights = new double[sliceTimes.length];
                for (int i = 0; i < sliceTimes.length; i++) {
                    if (mrsd == 0) {
                        sliceHeights[i] = sliceTimes[i];
                    } else {
                        sliceHeights[i] = mrsd - sliceTimes[i];
                    }
                }
            }

            String sliceHeightString = arguments.getStringOption(SLICE_HEIGHTS);
            if (sliceHeightString != null) {
                if (sliceTimeString != null) {
                    progressStream.println("Either sliceTimes, sliceHeights, timesFile or sliceCount" +
                            "nt.");
                    System.exit(-1);
                }
                sliceHeights = parseVariableLengthDoubleArray(sliceHeightString);
            }

            if (arguments.hasOption(SLICE_COUNT)) {
                int sliceCount = arguments.getIntegerOption(SLICE_COUNT);
                double startTime = arguments.getRealOption(START_TIME);
                double delta;
                if (mrsd != 0) {
                    delta = (mrsd - startTime) / (sliceCount - 1);
                } else {
                    delta = startTime / (sliceCount - 1);
                }
                sliceHeights = new double[sliceCount];
//                double height = mrsd - startTime;
//                for (int i = 0; i < sliceCount; i++) {
//                   sliceHeights[i] = height;
//                   height -= delta;
//                }
                double height = 0;
                for (int i = 0; i < sliceCount; i++) {
                    sliceHeights[i] = height;
                    height += delta;
                }
            }

            String optionString = arguments.getStringOption(ROOT);
            if (optionString != null && optionString.compareToIgnoreCase("true") == 0) {
                summarizeRoot = true;
            }

            optionString = arguments.getStringOption(TIPS);
            if (optionString != null && optionString.compareToIgnoreCase("true") == 0) {
                summarizeTips = true;
            }

            optionString = arguments.getStringOption(CONTOURS);
            if (optionString != null && optionString.compareToIgnoreCase("false") == 0) {
                contours = false;
            }

            optionString = arguments.getStringOption(POINTS);
            if (optionString != null && optionString.compareToIgnoreCase("true") == 0) {
                points = true;
            }

            //sorting sliceHeights
            if (sliceHeights != null) {
                Arrays.sort(sliceHeights);
            }

            String imputeString = arguments.getStringOption(IMPUTE);
            if (imputeString != null && imputeString.compareToIgnoreCase("true") == 0)
                impute = true;

            String branchNormString = arguments.getStringOption(BRANCH_NORMALIZE);
            if (branchNormString != null && branchNormString.compareToIgnoreCase("true") == 0)
                branchNormalization = true;

            String noiseString = arguments.getStringOption(NOISE);
            if (noiseString != null && noiseString.compareToIgnoreCase("false") == 0)
                trueNoise = false;

            String summaryString = arguments.getStringOption(SUMMARY);
            if (summaryString != null && summaryString.compareToIgnoreCase("true") == 0)
                summaryOnly = true;

            String modeString = arguments.getStringOption(CONTOUR_MODE);
            if (modeString != null) {
                contourMode = ContourMode.valueOf(modeString.toUpperCase());
                if (contourMode == ContourMode.R && !ContourWithR.processWithR)
                    contourMode = ContourMode.SNYDER;
            }

            String normalizeString = arguments.getStringOption(NORMALIZATION);
            if (normalizeString != null) {
                normalize = Normalization.valueOf(normalizeString.toUpperCase());
            }

            String summaryFormat = arguments.getStringOption(FORMAT);
            if (summaryFormat != null) {
                outputFormat = OutputFormat.valueOf(summaryFormat.toUpperCase());
            }

            String sdrString = arguments.getStringOption(SDR);
            if (sdrString != null) {
                outputFileSDR = sdrString;
                getSDR = true;
            }

            String snrString = arguments.getStringOption(SNR);
            if (snrString != null) {
                outputFileSNR = snrString;
                getSNR = true;
            }

            String progressString = arguments.getStringOption(PROGRESS);
            if (progressString != null) {
                progress = progressString;
            }

            String branch = arguments.getStringOption(BRANCHSET);
            if (branch != null) {
                set = BranchSet.valueOf(branch.toUpperCase());
                System.out.println("Using the branch set: " + set.name());
            }
            if (set == set.BACKBONE) {
                if (arguments.hasOption(BACKBONETAXA)) {
                    taxaSet = getTargetSet(arguments.getStringOption(BACKBONETAXA));
                } else {
                    System.err.println("you want to summarize the backbone, but have no taxa to define it??");
                }
            }
            if (set == set.CLADE) {
                if (arguments.hasOption(CLADETAXA)) {
                    taxaSet = getTargetSet(arguments.getStringOption(BACKBONETAXA));
                } else {
                    System.err.println("you want to get summaries for a clade, but have no taxa to define it??");
                }

            }

            if (arguments.hasOption(GRIDSIZE)) {
                grid = arguments.getIntegerOption(GRIDSIZE);
            }



        } catch (Arguments.ArgumentException e) {
            progressStream.println(e);
            printUsage(arguments);
            System.exit(-1);
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

        TimeSlicer timeSlicer = new TimeSlicer(inputFileName, burnin, skipEvery, traitNames, sliceHeights, impute,
                trueNoise, mrsd, contourMode, sliceMode,summarizeRoot, summarizeTips, normalize, getSDR, getSNR, progress,
                branchNormalization, set, taxaSet, grid, latMin, latMax, longMin, longMax, descendents, rateString);
        timeSlicer.output(outputFileName, summaryOnly, summarizeRoot, summarizeTips, contours, points, outputFormat, hpdValues, outputFileSDR, outputFileSNR);

        System.exit(0);
    }
}

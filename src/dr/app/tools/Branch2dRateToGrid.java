/*
 * Branch2dRateToGrid.java
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
import dr.geo.color.ChannelColorScheme;
import dr.geo.color.ColorScheme;
import dr.geo.math.SphericalPolarCoordinates;
import dr.math.distributions.MultivariateNormalDistribution;
import dr.util.TIFFWriter;
import dr.util.Version;

import javax.imageio.ImageIO;
import java.io.*;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: phil
 * Date: 08/01/2013
 * Time: 08:27
 * To change this template use File | Settings | File Templates.
 */
public class Branch2dRateToGrid {

    public static final boolean GREATCIRCLEDISTANCE = true;
    public static final String PRECISION_STRING = "precision";
    public static final String BURNIN = "burnin";
    public static final String SKIP = "skip";
    public static final String RATE_ATTRIBUTE = "rateAttribute";
    public static final String TRAIT = "trait";
    public static final String HELP = "help";
    public static final String TRAIT_NOISE = "traitNoise";
    public static final String RATE_NOISE = "rateNoise";
    public static final String FORMAT = "format";
    public static final String NORMALIZATION = "normalization";
    public static final String LATMAX = "latmax";
    public static final String LATMIN = "latmin";
    public static final String LONGMAX = "longmax";
    public static final String LONGMIN = "longmin";
    public static final String[] falseTrue = {"false", "true"};
    public static final String GRIDXCELLS = "gridXcells";
    public static final String GRIDYCELLS = "gridYcells";
    public static final String MAXPATHLENGTH = "maxPathLength";
    public static final String sep = "\t";
    //slicing
    public static final String SLICE_TIMES = "sliceTimes";
    public static final String SLICE_HEIGHTS = "sliceHeights";
    public static final String SLICE_FILE_HEIGHTS = "sliceFileHeights";
    public static final String SLICE_FILE_TIMES = "sliceFileTimes";
    public static final String SLICE_COUNT = "sliceCount";
    public static final String START_TIME = "startTime";
    public static final String STDEVS = "stdevs";
    public static final String CUTOFF = "cutoff";
    public static final String DISCRETE_TRAIT_NAME = "discreteTraitName";
    public static final String DISCRETE_TRAIT_STATES = "discreteTraitStates";
    public static final String HISTORY_ANNOTATION = "history";


    public Branch2dRateToGrid(String treeFileName, int burnin, int skipEvery,
                              String traitString, boolean traitNoise, boolean rateNoise, String rateString, double maxPathLength, Normalization normalize,
                              double latMin, double latMax, double longMin, double longMax, int gridXcells, int gridYcells, double[] sliceHeights, boolean getStdevs,
                              double posteriorCutoff, String discreteTrait, String[] discreteTraitStates, String historyAnnotation) {


        this.latMin = latMin;
        this.latMax = latMax;
        this.longMin = longMin;
        this.longMax = longMax;
        this.gridXcells = gridXcells;
        this.gridYcells = gridYcells;

        cellXWidth = (longMax - longMin) / gridXcells;
        cellYHeight = (latMax - latMin) / gridYcells;

        rateAttributeString = rateString;

        //slicing
//        mostRecentSamplingDate = mrsd;
        sliceCount = sliceHeights.length;

        if (sliceCount > 1) {
            sliceHeights = extractUnique(sliceHeights);
            Arrays.sort(sliceHeights);
            reverse(sliceHeights);
        }

        this.sliceHeights = sliceHeights;
        this.getStdevs = getStdevs;

        if (discreteTrait != null) {
            this.discreteTrait = discreteTrait;
            this.discreteTraitStates = discreteTraitStates;
            this.historyAnnotation = historyAnnotation;
            summarizeByDiscreteTrait = true;
        }

        densities = new double[sliceCount][gridXcells][gridYcells];
        if (summarizeByDiscreteTrait) {
            densitiesByDTrait = new double[sliceCount][discreteTraitStates.length][gridXcells][gridYcells];
        }
        rates = new double[sliceCount][gridXcells][gridYcells];
        if (getStdevs) {
            stdevs = new double[sliceCount][gridXcells][gridYcells];
        }

        this.posteriorCutoff = posteriorCutoff;

        try {
            readAndAnalyzeTrees(treeFileName, burnin, skipEvery, traitString, traitNoise, rateNoise, maxPathLength, normalize, false);
            summarizeMeanRates();
            if (getStdevs) {
//                progressStream.println("\n");
                readAndAnalyzeTrees(treeFileName, burnin, skipEvery, traitString, traitNoise, rateNoise, maxPathLength, normalize, true);
                summarizeStdevs();
            }
        } catch (IOException e) {
            System.err.println("Error reading file: " + treeFileName);
            System.exit(-1);
        } catch (Importer.ImportException e) {
            System.err.println("Error parsing trees in file: " + treeFileName);
            System.exit(-1);
        }

        progressStream.println(treesRead + " trees read.");
        progressStream.println(treesAnalyzed + " trees analyzed.");

    }


    private void readAndAnalyzeTrees(String treeFileName, int burnin, int skipEvery,
                                     String traitString, boolean traitNoise, boolean rateNoise, double maxPathLength, Normalization normalize,
                                     boolean getStdev)
            throws IOException, Importer.ImportException {

        int totalTrees = 10000;
        int totalStars = 0;

        if (!printedBar) {
            progressStream.println("Reading and analyzing trees (bar assumes 10,000 trees)...");
            progressStream.println("0              25             50             75            100");
            progressStream.println("|--------------|--------------|--------------|--------------|");
            printedBar = true;
        }

        if (getStdev) {
            progressStream.println("summarizing standard deviations");
        }

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
                    analyzeTree(treeTime, traitString, traitNoise, rateNoise, maxPathLength, normalize, getStdev);
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

    private void analyzeTree(Tree treeTime, String traitString, boolean traitNoise, boolean rateNoise, double maxPathLength, Normalization normalize, boolean getStdev) {

        double[][] precision = null;

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

        treeLengths.add(TreeUtils.getTreeLength(treeTime, treeTime.getRoot()));

        for (int x = 0; x < treeTime.getNodeCount(); x++) {

            NodeRef node = treeTime.getNode(x);

            if (!(treeTime.isRoot(node))) {

                double nodeHeight = treeTime.getNodeHeight(node);
                double parentHeight = treeTime.getNodeHeight(treeTime.getParent(node));
                Object tmpNodeTrait = treeTime.getNodeAttribute(node, traitString);
                if (tmpNodeTrait == null) {
                    System.err.println("Trait '" + traitString + "' not found on branch.");
                    System.exit(-1);
                }
                Trait nodeTrait = new Trait(tmpNodeTrait);

                Object tmpParentNodeTrait = treeTime.getNodeAttribute(treeTime.getParent(node), traitString);
                if (tmpParentNodeTrait == null) {
                    System.err.println("Trait '" + traitString + "' not found on branch.");
                    System.exit(-1);
                }
                Trait parentNodeTrait = new Trait(tmpParentNodeTrait);

                double rate = 1;
                if (!rateNoise) {
                    if (GREATCIRCLEDISTANCE) {
                        rate = getGeographicalDistance(nodeTrait.getValue(), parentNodeTrait.getValue()) / (parentHeight - nodeHeight);
                    } else {
                        rate = getNativeDistance(nodeTrait.getValue(), parentNodeTrait.getValue()) / (parentHeight - nodeHeight);
                    }
                }


                History history = null;
                if (summarizeByDiscreteTrait) {
                    Object historyObject = treeTime.getNodeAttribute(node, historyAnnotation);
                    String nodeState = ((String) treeTime.getNodeAttribute(node, discreteTrait)).replaceAll("\"", "");
                    String parentNodeState = ((String) treeTime.getNodeAttribute(treeTime.getParent(node), discreteTrait)).replaceAll("\"", "");

                    if (!parentNodeState.equals(nodeState)) {
                        if (historyObject instanceof Object[]) {
                            Object[] histories = (Object[]) historyObject;
                            history = setUpHistory(histories, nodeState, parentNodeState, nodeHeight, parentHeight);
                        } else {
                            System.err.println("History '" + historyAnnotation + "' not found on a branch that has different node and parent node states.");
                            System.exit(-1);
                        }
                    } else {
                        history = setUpHistory(nodeState, nodeHeight, parentHeight);
                    }
                }

                // slicing
                boolean inSliceInterval = false;

                for (int i = 0; i < sliceCount; i++) {
                    double heightAbove = Double.MAX_VALUE;
                    if (i > 0) {
                        heightAbove = sliceHeights[i - 1];
                    }

                    //branch entirely in interval
                    if ((nodeHeight > sliceHeights[i]) && (parentHeight < heightAbove)) {
                        inSliceInterval = true;
                    }
                    if ((parentHeight > heightAbove) && (nodeHeight < heightAbove)) {
                        inSliceInterval = true;
                        //impute parentNode
                        Double rateAttribute = (Double) treeTime.getNodeAttribute(node, rateAttributeString);
                        double diffusionRate = 1.0;
                        if (rateAttribute != null) {
                            diffusionRate = rateAttribute;
                            if (outputRateWarning) {
                                progressStream.println("Warning: using " + rateAttributeString + " as rate attribute during imputation!");
                                outputRateWarning = false;
                            }
                        }
                        if (traitNoise && precision == null) {
                            progressStream.println("Error: not precision available for imputation with correct noise!");
                            System.exit(-1);
                        }
                        parentNodeTrait = imputeValue(nodeTrait, parentNodeTrait, heightAbove, nodeHeight, parentHeight, precision, diffusionRate, traitNoise);
                        parentHeight = heightAbove;

                        // no truncation needed as all we eventually need to do is to find in what state we are at a particular point in time in that history, may need to be reconsidered for other purposes
//                        if (summarizeByDiscreteTrait) {
//                            history.truncateUpper(heightAbove);
//                        }

                    }
                    if ((nodeHeight < sliceHeights[i]) && (parentHeight > sliceHeights[i])) {
                        inSliceInterval = true;
                        //impute node
                        Double rateAttribute = (Double) treeTime.getNodeAttribute(node, rateAttributeString);
                        double diffusionRate = 1.0;
                        if (rateAttribute != null) {
                            diffusionRate = rateAttribute;
                            if (outputRateWarning) {
                                progressStream.println("Warning: using " + rateAttributeString + " as rate attribute during imputation!");
                                outputRateWarning = false;
                            }
                        }
                        if (traitNoise && precision == null) {
                            progressStream.println("Error: no precision available for imputation with correct noise!");
                            System.exit(-1);
                        }
                        nodeTrait = imputeValue(nodeTrait, parentNodeTrait, heightAbove, nodeHeight, parentHeight, precision, diffusionRate, traitNoise);
                        nodeHeight = sliceHeights[i];

                        // no truncation needed as all we eventually need to do is to find in what state we are at a particular point in time in that history, may need to be reconsidered for other purposes
//                        if (summarizeByDiscreteTrait) {
//                            history.truncateLower(sliceHeights[i]);
//                        }
                    }


                    if ((isInGrid(nodeTrait.getValue())) && (isInGrid(parentNodeTrait.getValue())) && inSliceInterval) {
                        //                    System.out.println("found a branch in the grid");

                        double branchLength = parentHeight - nodeHeight;
                        while (branchLength > maxPathLength) {

                            Double rateAttribute = (Double) treeTime.getNodeAttribute(node, rateAttributeString);
                            double diffusionRate = 1.0;
                            if (rateAttribute != null) {
                                diffusionRate = rateAttribute;
                                if (outputRateWarning) {
                                    progressStream.println("Warning: using " + rateAttributeString + " as rate attribute during imputation!");
                                    outputRateWarning = false;
                                }
                            }
                            if (traitNoise && precision == null) {
                                progressStream.println("Error: no precision available for imputation with correct noise!");
                                System.exit(-1);
                            }

                            double imputeTime = parentHeight - maxPathLength;
                            Trait intermediateTrait = imputeValue(nodeTrait, parentNodeTrait, imputeTime, nodeHeight, parentHeight, precision, diffusionRate, traitNoise);

                            if ((isInGrid(intermediateTrait.getValue())) && (isInGrid(parentNodeTrait.getValue()))) {

                                // no truncation needed as all we eventually need to do is to find in what state we are at a particular point in time in that history, may need to be reconsidered for other purposes
//                                if (summarizeByDiscreteTrait) {
//                                    history.truncateUpper(imputeTime);
//                                }

                                int[] nodeGridCell = getCellforPoint(intermediateTrait.getValue());
                                int[] parentNodeGridCell = getCellforPoint(parentNodeTrait.getValue());

                                if (rateNoise) {
                                    if (precision == null) {
                                        progressStream.println("Error: no precision available for imputation with correct noise!");
                                        System.exit(-1);
                                    }
                                    if (GREATCIRCLEDISTANCE) {
                                        rate = getGeographicalDistance(intermediateTrait.getValue(), parentNodeTrait.getValue()) / (parentHeight - imputeTime);
                                    } else {
                                        rate = getNativeDistance(intermediateTrait.getValue(), parentNodeTrait.getValue()) / (parentHeight - imputeTime);
                                    }
                                }

                                if (!getStdev) {
                                    addDensityAndRate(parentNodeGridCell[0], parentNodeGridCell[1], nodeGridCell[0], nodeGridCell[1], i, rate, history, parentHeight, nodeHeight);
                                } else {
                                    addStdev(parentNodeGridCell[0], parentNodeGridCell[1], nodeGridCell[0], nodeGridCell[1], i, rate);
                                }

                                branchLength = branchLength - maxPathLength;
                                parentHeight = parentHeight - maxPathLength;
                                parentNodeTrait = intermediateTrait;
                            }
                        }

                        int[] nodeGridCell = getCellforPoint(nodeTrait.getValue());
                        int[] parentNodeGridCell = getCellforPoint(parentNodeTrait.getValue());

                        if (rateNoise) {
                            if (precision == null) {
                                progressStream.println("Error: no precision available for imputation with correct noise!");
                                System.exit(-1);
                            }
                            if (GREATCIRCLEDISTANCE) {
                                rate = getGeographicalDistance(nodeTrait.getValue(), parentNodeTrait.getValue()) / (parentHeight - nodeHeight);
                            } else {
                                rate = getNativeDistance(nodeTrait.getValue(), parentNodeTrait.getValue()) / (parentHeight - nodeHeight);
                            }
                        }

                        if (!getStdev) {
                            addDensityAndRate(parentNodeGridCell[0], parentNodeGridCell[1], nodeGridCell[0], nodeGridCell[1], i, rate, history, parentHeight, nodeHeight);
                        } else {
                            addStdev(parentNodeGridCell[0], parentNodeGridCell[1], nodeGridCell[0], nodeGridCell[1], i, rate);
                        }

                    }
                }
            }
        }
        treesAnalyzed++;
    }

    private int[] getCellforPoint(double[] location) {
        // For point coordinates (x,y) and a grid with origin (Ox,Oy) and cellsize c, the grid coordinates are found by rounding (x-Ox)/c and (y-Oy)/c down to the nearest integer.
        int[] cell = new int[2];
        cell[0] = (int) Math.floor((location[1] - longMin) / cellXWidth);
        cell[1] = gridYcells - (int) Math.ceil((location[0] - latMin) / cellYHeight);
        return cell;
    }

    public void addDensityAndRate(int x, int y, int x2, int y2, int slice, double rate) {
        int w = x2 - x;
        int h = y2 - y;
        int dx1 = 0, dy1 = 0, dx2 = 0, dy2 = 0;
        if (w < 0) dx1 = -1;
        else if (w > 0) dx1 = 1;
        if (h < 0) dy1 = -1;
        else if (h > 0) dy1 = 1;
        if (w < 0) dx2 = -1;
        else if (w > 0) dx2 = 1;
        int longest = Math.abs(w);
        int shortest = Math.abs(h);
        if (!(longest > shortest)) {
            longest = Math.abs(h);
            shortest = Math.abs(w);
            if (h < 0) dy2 = -1;
            else if (h > 0) dy2 = 1;
            dx2 = 0;
        }
        int numerator = longest >> 1;
        for (int i = 0; i <= longest; i++) {
//            putpixel(x,y,color) ;
            densities[slice][x][y]++;
            rates[slice][x][y] = +rate;
            numerator += shortest;
            if (!(numerator < longest)) {
                numerator -= longest;
                x += dx1;
                y += dy1;
            } else {
                x += dx2;
                y += dy2;
            }
        }
    }

    public void addDensityAndRate(int x, int y, int x2, int y2, int slice, double rate, History history, double parentHeight, double nodeHeight) {
        int w = x2 - x;
        int h = y2 - y;
        double branchLength = parentHeight - nodeHeight;
        int dx1 = 0, dy1 = 0, dx2 = 0, dy2 = 0;
        if (w < 0) dx1 = -1;
        else if (w > 0) dx1 = 1;
        if (h < 0) dy1 = -1;
        else if (h > 0) dy1 = 1;
        if (w < 0) dx2 = -1;
        else if (w > 0) dx2 = 1;
        int longest = Math.abs(w);
        int shortest = Math.abs(h);
        if (!(longest > shortest)) {
            longest = Math.abs(h);
            shortest = Math.abs(w);
            if (h < 0) dy2 = -1;
            else if (h > 0) dy2 = 1;
            dx2 = 0;
        }
        int numerator = longest >> 1;
        for (int i = 0; i <= longest; i++) {

            if (history != null) {
                double height = parentHeight;
                if (longest > 0) {
                    height -= (branchLength / longest) * i;
                }
                double[] test = history.getHeights();
                if (height < test[test.length - 1]) {
//                    System.err.println("encountered height (" + height + ") that is lower than lower of history heights (" + test[test.length - 1] + "). Possible rounding error -- setting the height to upper of history heights");
                    height = test[test.length - 1];
                }

                String state = history.getStateForHeight(height);
                int stateInt = getIntForState(state);
                if (stateInt >= 0) {
                    densitiesByDTrait[slice][stateInt][x][y]++;
                }
            }

            densities[slice][x][y]++;
            rates[slice][x][y] =+ rate;

            numerator += shortest;
            if (!(numerator < longest)) {
                numerator -= longest;
                x += dx1;
                y += dy1;
            } else {
                x += dx2;
                y += dy2;
            }
        }
    }

    public void addStdev(int x, int y, int x2, int y2, int slice, double rate) {
        int w = x2 - x;
        int h = y2 - y;
        int dx1 = 0, dy1 = 0, dx2 = 0, dy2 = 0;
        if (w < 0) dx1 = -1;
        else if (w > 0) dx1 = 1;
        if (h < 0) dy1 = -1;
        else if (h > 0) dy1 = 1;
        if (w < 0) dx2 = -1;
        else if (w > 0) dx2 = 1;
        int longest = Math.abs(w);
        int shortest = Math.abs(h);
        if (!(longest > shortest)) {
            longest = Math.abs(h);
            shortest = Math.abs(w);
            if (h < 0) dy2 = -1;
            else if (h > 0) dy2 = 1;
            dx2 = 0;
        }
        int numerator = longest >> 1;
        for (int i = 0; i <= longest; i++) {
//            System.out.println(i+"\t"+x+"\t"+y);
//            System.out.println(rates[slice][x][y]);
            stdevs[slice][x][y] = +Math.pow((rate - rates[slice][x][y]), 2);
            numerator += shortest;
            if (!(numerator < longest)) {
                numerator -= longest;
                x += dx1;
                y += dy1;
            } else {
                x += dx2;
                y += dy2;
            }
        }
    }

    private static double getNativeDistance(double[] location1, double[] location2) {
        return Math.sqrt(Math.pow((location2[0] - location1[0]), 2.0) + Math.pow((location2[1] - location1[1]), 2.0));
    }

    private static double getGeographicalDistance(double[] location1, double[] location2) {
        if (location1.length == 1) {
            // assume we only have latitude so put them on the prime meridian
            return getKilometerGreatCircleDistance(new double[]{location1[0], 0.0}, new double[]{location2[0], 0.0});
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

    private double[] imputeValue(double[] nodeValue, double[] parentValue, double time, double nodeHeight, double parentHeight, double[] precisionArray, double rate, boolean trueNoise) {

        final double scaledTimeChild = (time - nodeHeight) * rate;
        final double scaledTimeParent = (parentHeight - time) * rate;
        final double scaledWeightTotal = 1.0 / scaledTimeChild + 1.0 / scaledTimeParent;
        final int dim = nodeValue.length;

        double[][] precision = new double[dim][dim];
        int counter = 0;
        for (int a = 0; a < dim; a++) {
            for (int b = 0; b < dim; b++) {
                precision[a][b] = precisionArray[counter];
                counter++;
            }
        }

        if (scaledTimeChild == 0)
            return nodeValue;

        if (scaledTimeParent == 0)
            return parentValue;

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

        double[] result = new double[dim];
        for (int i = 0; i < dim; i++)
            result[i] = mean[i];
        return result;
    }

    private boolean isInGrid(double[] point) {
        boolean in = false;
        if (((point[0] > latMin) && (point[0] < latMax)) &&
                ((point[1] > longMin) && (point[1] < longMax))) {
            in = true;
        }
        return in;
    }

    private int getIntForState(String state) {

        int returnInt = -1;

        for (int x = 0; x < discreteTraitStates.length; x++) {
            if (state.equals(discreteTraitStates[x])) {
                returnInt = x;
            }

        }
//        if (returnInt < 0){
//            System.err.print("state "+state+" is not in the specified discrete states");
//        }
        return returnInt;
    }


    enum Normalization {
        LENGTH,
        HEIGHT,
        NONE
    }


    private String name(String pre, String post) {
        return pre + "." + post;
    }

    public void output(String outFileName) {

        resultsStream = System.out;

        if (outFileName != null) {
            try {
                resultsStream = new PrintStream(new File(outFileName));
            } catch (IOException e) {
                System.err.println("Error opening file: " + outFileName);
                System.exit(-1);
            }
        }

        outputGridInfo();
        resultsStream.print("\n");

        if (posteriorCutoff > 0) {
//            setDensityCutoff();
            setPosteriorCutoff();
        }

        double maxGridDensity = 0;
        if (sliceCount > 1) {
            for (int a = 0; a < sliceCount; a++) {
                double tempMax = getMaxMatrix(densities[a]);
                if (maxGridDensity < tempMax) {
                    maxGridDensity = tempMax;
                }
            }
        }

        double maxGridDensityByDtrait = 0;
        if (summarizeByDiscreteTrait) {
            for (int a = 0; a < sliceCount; a++) {
                for (int b = 0; b < discreteTraitStates.length; b++) {
                    double tempMax = getMaxMatrix(densitiesByDTrait[a][b]);
                    if (maxGridDensityByDtrait < tempMax) {
                        maxGridDensityByDtrait = tempMax;
                    }
                }
            }
        }


        String writerNames[] = ImageIO.getWriterFormatNames();
        for (String name : writerNames) {
            System.err.println("Available format:" + name);
        }
        // set graphicsFormat equal to one of the names get that's printed out above.
        String graphicsFormat = "png";
//        String graphicsFormat = "tiff";
        String suffix = ".worldFile";
        if (graphicsFormat.equals("png")){
            suffix = "pgw";
        }

        for (int i = 0; i < sliceCount; i++) {

            String ratesFile = "gridRates";
            String densityFile = "gridDensity";
            String rateStdevFile = "gridRateStdevs";

            if (sliceCount == 1) {
                //all these tranposation below to accomodate printing X coordinates row by row
                resultsStream.print("grid rates:\n");
                printGrid(transpose(rates[i]));
                resultsStream.print("\n");

                resultsStream.print("grid densities:\n");
                printGrid(transpose(densities[i]));
                resultsStream.print("\n");

                if (getStdevs) {
                    resultsStream.print("grid stdevs:\n");
                    printGrid(transpose(stdevs[i]));
                }
                resultsStream.print("\n");

//                writeAsTIFF(ratesFile + ".tiff", transpose(rates[i]), true);
//                writeWorldFile(ratesFile + ".tfw", cellXWidth, cellYHeight, 0, 0, longMin, latMax);
//                writeAsTIFF(densityFile + ".tiff", transpose(densities[i]), true);
//                writeWorldFile(densityFile + ".tfw", cellXWidth, cellYHeight, 0, 0, longMin, latMax);
//                writeAsTIFF(rateStdevFile + ".tiff", transpose(stdevs[i]), true);
//                writeWorldFile(rateStdevFile + ".tfw", cellXWidth, cellYHeight, 0, 0, longMin, latMax);

//                System.err.println("Start PNGs");

                writeAsAnyFormat(name(ratesFile, graphicsFormat), graphicsFormat, rates[i], true);
                writeAsAnyFormat(name(densityFile, graphicsFormat), graphicsFormat, densities[i], true);
                writeAsAnyFormat(name(rateStdevFile, graphicsFormat), graphicsFormat, stdevs[i], true);

                writeWorldFile(name(ratesFile, suffix), cellXWidth, cellYHeight, 0, 0, longMin, latMax);
                writeWorldFile(name(densityFile, suffix), cellXWidth, cellYHeight, 0, 0, longMin, latMax);
                writeWorldFile(name(rateStdevFile, suffix), cellXWidth, cellYHeight, 0, 0, longMin, latMax);

//                System.err.println("End PNGs");

                if (summarizeByDiscreteTrait) {
                    for (int x = 0; x < discreteTraitStates.length; x++) {
                        resultsStream.print("grid densities for discrete trait " + discreteTraitStates[x] + ":\n");
                        printGrid(transpose(densitiesByDTrait[i][x]));
                        resultsStream.print("\n");

//                        writeAsTIFF(densityFile + ".discreteTrait" + discreteTraitStates[x] + ".tiff", transpose(densitiesByDTrait[i][x]), true, maxGridDensityByDtrait);
                        writeAsAnyFormat(name(densityFile + ".discreteTrait" + discreteTraitStates[x], graphicsFormat),
                                graphicsFormat,
                                densitiesByDTrait[i][x], true, maxGridDensityByDtrait);
                        writeWorldFile(name(densityFile + ".discreteTrait" + discreteTraitStates[x], suffix), cellXWidth, cellYHeight, 0, 0, longMin, latMax);
                    }

                    // Try multiple channels
                    List<double[][]> channels = new ArrayList<double[][]>();
                    for (int x = 0; x < discreteTraitStates.length; ++x) {
                        channels.add(densitiesByDTrait[i][x]);
                    }
                    writeAsAnyFormatMultiChannel(name("channel." + densityFile + ".discreteTraitAll", graphicsFormat),
                            graphicsFormat,
                            channels, true, maxGridDensityByDtrait, ChannelColorScheme.CHANNEL_RED_BLUE);
                    writeWorldFile(name("channel." + densityFile + ".discreteTraitAll", suffix), cellXWidth, cellYHeight, 0, 0, longMin, latMax);

                }

            } else {
                resultsStream.print("grid rates for slice height " + sliceHeights[i] + ":\n");
                printGrid(transpose(rates[i]));
                resultsStream.print("\n");

                resultsStream.print("grid densities for slice height " + sliceHeights[i] + ":\n");
                printGrid(transpose(densities[i]));
                resultsStream.print("\n");

                if (getStdevs){
                    resultsStream.print("grid rate stdevs for slice height " + sliceHeights[i] + ":\n");
                    printGrid(transpose(stdevs[i]));
                    resultsStream.print("\n");
                }

//                writeAsTIFF(ratesFile + ".height" + sliceHeights[i] + ".tiff", transpose(rates[i]), true, maxGridDensity);
//                writeAsTIFF(densityFile + ".height" + sliceHeights[i] + ".tiff", transpose(densities[i]), truemaxGridDensity);
//                writeAsTIFF(rateStdevFile + ".height" + sliceHeights[i] + ".tiff", transpose(stdevs[i]), true, maxGridDensity);
                writeAsAnyFormat(name(ratesFile + ".height" + sliceHeights[i], graphicsFormat), graphicsFormat, rates[i], true, maxGridDensity);
                writeAsAnyFormat(name(densityFile + ".height" + sliceHeights[i], graphicsFormat), graphicsFormat, densities[i], true, maxGridDensity);

                writeWorldFile(name(ratesFile + ".height" + sliceHeights[i], suffix), cellXWidth, cellYHeight, 0, 0, longMin, latMax);
                writeWorldFile(name(densityFile + ".height" + sliceHeights[i], suffix), cellXWidth, cellYHeight, 0, 0, longMin, latMax);

                if (getStdevs){
                    writeAsAnyFormat(name(rateStdevFile + ".height" + sliceHeights[i], graphicsFormat), graphicsFormat, stdevs[i], true, maxGridDensity);
                    writeWorldFile(name(rateStdevFile + ".height" + sliceHeights[i], suffix), cellXWidth, cellYHeight, 0, 0, longMin, latMax);
                }

                if (summarizeByDiscreteTrait) {
                    for (int x = 0; x < discreteTraitStates.length; x++) {
                        resultsStream.print("grid densities for slice height " + sliceHeights[i] + " and for discrete trait " + discreteTraitStates[x] + ":\n");
                        printGrid(transpose(densitiesByDTrait[i][x]));
                        resultsStream.print("\n");

//                        writeAsTIFF(densityFile + ".height" + sliceHeights[i] + ".discreteTrait" + discreteTraitStates[x] + ".tiff", transpose(densitiesByDTrait[i][x]), true, maxGridDensityByDtrait);
                        writeAsAnyFormat(name(densityFile + ".height" + sliceHeights[i] + ".discreteTrait" + discreteTraitStates[x], graphicsFormat), graphicsFormat, densitiesByDTrait[i][x], true, maxGridDensityByDtrait);
                        writeWorldFile(name(densityFile + ".height" + sliceHeights[i] + ".discreteTrait" + discreteTraitStates[x], suffix), cellXWidth, cellYHeight, 0, 0, longMin, latMax);
                    }

                    // Try multiple channels
                    List<double[][]> channels = new ArrayList<double[][]>();
                    for (int x = 0; x < discreteTraitStates.length; ++x) {
                        channels.add(densitiesByDTrait[i][x]);
                    }
                    writeAsAnyFormatMultiChannel(name("channel." + densityFile + ".height" + sliceHeights[i] + ".discreteTraitAll", graphicsFormat),
                            graphicsFormat,
                            channels, true, maxGridDensityByDtrait, ChannelColorScheme.CHANNEL_RED_BLUE);
                    writeWorldFile(name("channel." + densityFile + ".height" + sliceHeights[i] + ".discreteTraitAll", suffix), cellXWidth, cellYHeight, 0, 0, longMin, latMax);

                }
            }
        }
    }

    private void outputGridInfo() {
        StringBuffer sb = new StringBuffer("# grid info\n");
        sb.append("# lat min (Y)").append(sep).append(latMin);
        sb.append("\n");
        sb.append("# lat max (Y)").append(sep).append(latMax);
        sb.append("\n");
        sb.append("# long min (Y)").append(sep).append(longMin);
        sb.append("\n");
        sb.append("# long max (Y)").append(sep).append(longMax);
        sb.append("\n");
        sb.append("# number of X cells (long)").append(sep).append(gridXcells);
        sb.append("\n");
        sb.append("# number of Y cells (lat)").append(sep).append(gridYcells);
        sb.append("\n");
        sb.append("# X cell width (long)").append(sep).append(cellXWidth);
        sb.append("\n");
        sb.append("# Y cell height (lat)").append(sep).append(cellYHeight);
        sb.append("\n");
        if (posteriorCutoff > 0) {
            sb.append("# cut off for reporting grid values").append(sep).append(posteriorCutoff);
            sb.append("\n");
        }
        resultsStream.print(sb);
    }

    private void printGrid(double array[][]) {

        StringBuffer sb = new StringBuffer();
        for (int a = 0; a < array.length; a++) {
            for (int b = 0; b < array[0].length; b++) {
                sb.append(array[a][b]).append(sep);
            }
            sb.append("\n");
        }
        resultsStream.print(sb);
    }

    public void writeAsTIFF(String fileName, double[][] matrix, boolean log) {

//        System.out.print("matrix x = "+matrix.length+"; "+"y = "+matrix[0].length);

        double[][] mat = normalize(matrix, 255, log);
        try {
            DataOutputStream tiffOut = new DataOutputStream(new FileOutputStream(fileName));
            TIFFWriter.writeDoubleArray(tiffOut, mat);
            tiffOut.close();
        } catch (Exception e) {
            throw new RuntimeException(e.getMessage());
        }
    }


    // format = "png", = "tiff" = "gif", etc.

    public void writeAsAnyFormat(String fileName, String format, double[][] matrix, boolean log) {
        double[][] mat = normalize(matrix, 255, log);
        try {
            TIFFWriter.writeDoubleArray(fileName, mat, format, ColorScheme.HEATMAP);
        } catch (Exception e) {
            throw new RuntimeException(e.getMessage());
        }
    }

    public void writeAsTIFF(String fileName, double[][] matrix, boolean log, double maxValue) {

//        System.out.print("matrix x = "+matrix.length+"; "+"y = "+matrix[0].length);

        double[][] mat = normalize(matrix, 255, log, maxValue);
        try {
            DataOutputStream tiffOut = new DataOutputStream(new FileOutputStream(fileName));
            TIFFWriter.writeDoubleArray(tiffOut, mat);
            tiffOut.close();
        } catch (Exception e) {
            throw new RuntimeException(e.getMessage());
        }
    }

    public void writeAsAnyFormat(String fileName, String format, double[][] matrix, boolean log, double maxValue) {

//        System.out.print("matrix x = "+matrix.length+"; "+"y = "+matrix[0].length);

        double[][] mat = normalize(matrix, 255, log, maxValue);
        try {
            TIFFWriter.writeDoubleArray(fileName, mat, format,
//                    TIFFWriter.HEATMAP
//                    TIFFWriter.TRANPARENT_HEATMAP2  // Assumes NaN is transparent
                    ColorScheme.TRANPARENT0_HEATMAP // Assumes 0.0 is transparent
            );
        } catch (Exception e) {
            throw new RuntimeException(e.getMessage());
        }
    }

    public void writeAsAnyFormatMultiChannel(String fileName, String format, List<double[][]> matrices, boolean log, double maxValue,
                                             ChannelColorScheme scheme) {

//        System.out.print("matrix x = "+matrix.length+"; "+"y = "+matrix[0].length);

        List<double[][]> normalizedMatrix = new ArrayList<double[][]>();
        for (double[][] matrix : matrices) {
            normalizedMatrix.add(normalize(matrix, 255, log, maxValue)); // TODO Maybe we want different maxValue for each matrix?
        }

        try {
            TIFFWriter.writeDoubleArrayMultiChannel(fileName, normalizedMatrix, format, scheme);
        } catch (Exception e) {
            throw new RuntimeException(e.getMessage());
        }
    }

    public void writeWorldFile(String fileName, double xDim, double yDim, double rot1, double rot2, double upperLeftLong, double upperLeftLat) {

        try {
            PrintWriter outFile = new PrintWriter(new FileWriter(fileName), true);

            outFile.println(xDim);
            outFile.println(rot1);
            outFile.println(rot2);
            outFile.println(-yDim);
            outFile.println(upperLeftLong+(xDim/2.0));
            outFile.println(upperLeftLat-(yDim/2.0));

            outFile.close();

        } catch(IOException io) {
            System.err.print("Error writing to file: " + fileName);
        }

    }

    // TODO Fix function to handle matrices with NaN entries
    private double sum(double[][] mat) {
        double value = 0.0;
        for (int i = 0; i < mat.length; ++i) {
            for (int j = 0; j < mat[i].length; ++j) {
                value += mat[i][j];
            }
        }
        return value;
    }

    // TODO Fix function to handle matrices with NaN entries
    // TODO Or: remove code duplication in function immediately below
    private double[][] normalize(double inputArray[][], double max, boolean log) {

        double[][] matrix = new double[inputArray.length][inputArray[0].length];

        double maxValue = 0;
        for (int i = 0; i < inputArray[0].length; i++) {
            for (int j = 0; j < inputArray.length; j++) {
                if (log) {
                    //pseudocount
                    if (inputArray[j][i] == 0) {
                        inputArray[j][i] = 1;
                    }
                    inputArray[j][i] = Math.log(inputArray[j][i]);
                }
                if (inputArray[j][i] > maxValue) {
                    maxValue = inputArray[j][i];
                }
            }
        }

        for (int i = 0; i < inputArray[0].length; i++) {
            for (int j = 0; j < inputArray.length; j++) {
                matrix[j][i] = ((double) inputArray[j][i] / maxValue) * max;
            }
        }

        return matrix;

    }

    // TODO Fix function to handle matrices with NaN entries
    private double[][] normalize(double inputArray[][], double max, boolean log, double maxValue) {

        double[][] matrix = new double[inputArray.length][inputArray[0].length];

        if (log) {
            maxValue = Math.log(maxValue);
            for (int i = 0; i < inputArray[0].length; i++) {
                for (int j = 0; j < inputArray.length; j++) {
                    //pseudocount
                    if (inputArray[j][i] == 0) {
                        matrix[j][i] = 0;
                    } else {
                        matrix[j][i] = Math.log(inputArray[j][i]);
                    }
                }
            }
        }


        for (int i = 0; i < inputArray[0].length; i++) {
            for (int j = 0; j < inputArray.length; j++) {
                matrix[j][i] = ((double) matrix[j][i] / maxValue) * max;
            }
        }

        return matrix;

    }

    private double[][] toDoubleArray(int array[][]) {

        double[][] returnArray = new double[array.length][array[0].length];
        for (int a = 0; a < returnArray.length; a++) {
            for (int b = 0; b < array[0].length; b++) {
                returnArray[a][b] = array[a][b];
            }
        }
        return returnArray;
    }

    public double[][] transpose(double[][] matrix) {

        double[][] transposeMatrix = new double[matrix[0].length][matrix.length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                transposeMatrix[j][i] = matrix[i][j];
            }

        }
        return transposeMatrix;
    }

    private void summarizeMeanRates() {
        for (int a = 0; a < rates.length; a++) {
            for (int b = 0; b < rates[0].length; b++) {
                for (int c = 0; c < rates[0][0].length; c++) {
                    rates[a][b][c] = rates[a][b][c] / densities[a][b][c];
                    if (rates[a][b][c] == 0){
                        rates[a][b][c] = Double.NaN;
                    }
                }
            }
        }
    }

    private void summarizeStdevs() {
        for (int a = 0; a < stdevs.length; a++) {
            for (int b = 0; b < stdevs[0].length; b++) {
                for (int c = 0; c < stdevs[0][0].length; c++) {
                    stdevs[a][b][c] = Math.sqrt(stdevs[a][b][c] / densities[a][b][c]);
                }
            }
        }
    }

    private double[][] summarizRateStdevs(List<Double>[][] individualRates, double[][] meanRates) {
        double[][] rateStdevs = new double[individualRates.length][individualRates[0].length];
        for (int i = 0; i < individualRates.length; i++) {
            for (int j = 0; j < individualRates[0].length; j++) {
                if (individualRates[i][j].isEmpty()) {
//                    sumOfrates = Double.NaN;
                    rateStdevs[i][j] = Double.NaN;
                } else {
                    List<Double> cellRates = individualRates[i][j];
                    double sumOfSquaredDifferences = 0;
                    for (Double element : cellRates) {
                        sumOfSquaredDifferences += Math.pow((element - meanRates[i][j]), 2);
                    }
                    rateStdevs[i][j] = Math.sqrt(sumOfSquaredDifferences / (cellRates.size()));
                }

            }
        }

        return rateStdevs;
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

    public static void reverse(double[] array) {
        if (array == null) {
            return;
        }
        int i = 0;
        int j = array.length - 1;
        double tmp;
        while (j > i) {
            tmp = array[j];
            array[j] = array[i];
            array[i] = tmp;
            j--;
            i++;
        }
    }

    public static double[] extractUnique(double[] array) {
        Set<Double> tmp = new LinkedHashSet<Double>();
        for (Double each : array) {
            tmp.add(each);
        }
        double[] output = new double[tmp.size()];
        int i = 0;
        for (Double each : tmp) {
            output[i++] = each;
        }
        return output;
    }


    public static <T extends Enum<T>> String[] enumNamesToStringArray(T[] values) {
        int i = 0;
        String[] result = new String[values.length];
        for (T value : values) {
            result[i++] = value.name();
        }
        return result;
    }

    public static double getMaxMatrix(double[][] gridDensities) {
        double max = 0;
        for (int x = 0; x < gridDensities.length; x++) {
            for (int y = 0; y < gridDensities[0].length; y++) {
                if (gridDensities[x][y] > max) {
                    max = gridDensities[x][y];
                }
            }

        }
        return max;
    }

    private void setPosteriorCutoff() {
// no discrete traits
//        for (int a = 0; a < densities.length; a++){
//            for (int b = 0; b < densities[0].length; b++){
//                for (int c = 0; c < densities[0][0].length; c++){
//                    if (densities[a][b][c] < (posteriorCutoff*treesAnalyzed)){
//                        densities[a][b][c] = 0;
//                        rates[a][b][c] = Double.NaN;
//                        stdevs[a][b][c] = Double.NaN;
//                    }
//                }
//            }
//        }

        for (int a = 0; a < densitiesByDTrait.length; a++) {
            for (int b = 0; b < densitiesByDTrait[0].length; b++) {
                for (int c = 0; c < densitiesByDTrait[0][0].length; c++) {
                    for (int d = 0; d < densitiesByDTrait[0][0][0].length; d++) {
                        if (summarizeByDiscreteTrait) {
                            if (densitiesByDTrait[a][b][c][d] < (posteriorCutoff * treesAnalyzed)) {
                                densitiesByDTrait[a][b][c][d] = 0;
                            }
                        }
                        if (densities[a][c][d] < (posteriorCutoff * treesAnalyzed)) {
                            densities[a][c][d] = 0;
                            rates[a][c][d] = Double.NaN;
                            if (getStdevs){
                                stdevs[a][c][d] = Double.NaN;
                            }
                        }
                    }
                }
            }
        }
    }

    public History setUpHistory(String nodeState, double timeLow, double timeUp) {
        double[] heights = new double[]{timeUp, timeLow};
        String[] states = new String[]{nodeState};
        return new History(heights, states);
    }

    public History setUpHistory(Object[] historyObject, String nodeState, String parentNodeState, double timeLow, double timeUp) {
        double[] heights;
        String[] states;

        Object[] histories = (Object[]) historyObject;
//        System.out.println(histories.length);
        Object[] testHistory = (Object[]) histories[0];
        String[][] jumpStrings = new String[histories.length][testHistory.length];


        for (int q = 0; q < histories.length; q++) {
            Object[] singleHistory = (Object[]) histories[q];
            for (int r = 0; r < singleHistory.length; r++) {
//                System.out.println(singleHistory[r]);
                jumpStrings[q][r] = singleHistory[r].toString();
            }

        }

        //sorting jumpStrings not necessary: jumps are in order of their occurrence
        //fill heights and states
        heights = new double[histories.length + 2];
        states = new String[histories.length + 1];
        for (int b = 0; b < histories.length; b++) {
            states[b] = jumpStrings[b][1];
            heights[b + 1] = Double.valueOf(jumpStrings[b][0]);
        }

        //sanity check
        if (!jumpStrings[0][1].equals(parentNodeState)) {
            System.out.println(jumpStrings[0][1] + "\t" + parentNodeState);
            System.err.println("mismatch in jump history and parent node state");
            System.exit(-1);
        }

        //sanity check
        states[histories.length] = jumpStrings[histories.length - 1][2];
        if (!jumpStrings[histories.length - 1][2].equals(nodeState)) {
            System.err.println("mismatch in jump history and node state");
            System.exit(-1);
        }

        heights[0] = timeUp;
        heights[histories.length + 1] = timeLow;


        return new History(heights, states);
    }

    private String getState(int stateInt) {
        String returnString = null;
        try {
//            returnString = markovJumpLikelihood.formattedState(new int[] {stateInt}).replaceAll("\"","");
        } catch (IndexOutOfBoundsException iobe) {
            System.err.println("no state found for int = " + stateInt + "...");
            System.exit(-1);
        }
        return returnString;
    }

    private int treesRead = 0;
    private int treesAnalyzed = 0;
    private double latMin;
    private double latMax;
    private double longMin;
    private double longMax;
    private int gridXcells;
    private int gridYcells;
    private double cellXWidth;
    private double cellYHeight;
    private String rateAttributeString = "rate";
    //    private double[][] gridMeanRates;
//    private int[][] gridHits;
    private ArrayList treeLengths = new ArrayList();
    private boolean outputRateWarning = true;
    private double[][][] densities;
    private double[][][][] densitiesByDTrait;
    private double[][][] rates;
    private double[][][] stdevs;
    // by discrete trait
//    private double[][][][] densities;
//    private double[][][][] rates;
//    private double[][][][] stdevs;
    private int sliceCount;
    private double[] sliceHeights;
    //    private double mostRecentSamplingDate;
    private boolean printedBar = false;
    private boolean getStdevs = true;
    private double posteriorCutoff;
    private String discreteTrait;
    private String[] discreteTraitStates;
    private boolean summarizeByDiscreteTrait = false;
    private String historyAnnotation;

    //    // Messages to stderr, output to stdout
    private static PrintStream progressStream = System.err;
    private PrintStream resultsStream;
    //
    private final static Version version = new BeastVersion();
    //
    private static final String commandName = "Branch2dRateToGrid";


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
        centreLine("branchGrid2Drate " + version.getVersionString() + ", " + version.getDateString(), 60);
        centreLine("MCMC Output analysis", 60);
        centreLine("by", 60);
        centreLine("Philippe Lemey, Marc A. Suchard", 60);
        progressStream.println();
        centreLine("Rega Institute for Medical Research", 60);
        centreLine("KU Leuven", 60);
        centreLine("philippe.lemey@gmail.com", 60);
        progressStream.println();
        centreLine("Department of Biomathematics", 60);
        centreLine("University of California, Los Angeles", 60);
        centreLine("msuchard@ucla.edu", 60);
        progressStream.println();
        progressStream.println();
    }


    public static void main(String[] args) throws IOException {

        String inputFileName = null;
        String outputFileName = null;
        String traitName = "location";
        boolean trait2DNoise = true;
        boolean rate2DNoise = false;
        Normalization normalize = Normalization.LENGTH;
        int burnin = -1;
        int skipEvery = 1;
        double latMax = 90;
        double latMin = -90;
        double longMax = 180;
        double longMin = -180;
        int gridXcells = 360;
        int gridYcells = 180;
        String rateString = "rate";
        double maxPathLength = Double.MAX_VALUE;
        double posteriorCutoff = 0;
        //slicing
        double[] sliceHeights = null;
        double mrsd = 0;
        boolean getStdevs = true;
        String discreteTraitName = null;
        String[] discreteTraitStates = null;
        String historyAnnotation = "history";

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.IntegerOption(BURNIN, "the number of states to be considered as 'burn-in' [default = 0]"),
                        new Arguments.IntegerOption(SKIP, "skip every i'th tree [default = 0]"),
                        new Arguments.StringOption(TRAIT, "trait_name", "specifies an attribute to use to summarize the 2D trait info [default = location]"),
                        new Arguments.StringOption(RATE_ATTRIBUTE, "rate_attribute", "specifies the trait rate attribute string [default=rate]"),
                        new Arguments.Option(HELP, "option to print this message"),
                        new Arguments.StringOption(TRAIT_NOISE, falseTrue, false,
                                "add true noise to 2D traits [default = true])"),
                        new Arguments.StringOption(RATE_NOISE, falseTrue, false,
                                "add true noise to rates [default = true])"),
                        new Arguments.StringOption(NORMALIZATION, enumNamesToStringArray(Normalization.values()), false,
                                "tree normalization [default = length"),
                        new Arguments.RealOption(LATMAX, "specifies the maximum latitude for the grid [default=90]"),
                        new Arguments.RealOption(LATMIN, "specifies the minimum latitude for the grid [default=-90]"),
                        new Arguments.RealOption(LONGMAX, "specifies the maximum longitude for the grid [default=180]"),
                        new Arguments.RealOption(LONGMIN, "specifies the minimum longitude for the grid [default=-180]"),
                        new Arguments.IntegerOption(GRIDXCELLS, "the number of cells along the x-axis (longitude) of the grid [default = 360]"),
                        new Arguments.IntegerOption(GRIDYCELLS, "the number of cells along the y-axis (latitude) of the grid [default = 180]"),
                        new Arguments.RealOption(MAXPATHLENGTH, "specifies the maximum (time) length a branch can go in one direction before adding Brownian noise to it [default=MAX_VALUE (no noise)]"),
                        new Arguments.StringOption(SLICE_TIMES, "slice_times", "specifies a slice time-list [default=none]"),
                        new Arguments.StringOption(SLICE_HEIGHTS, "slice_heights", "specifies a slice height-list [default=none]"),
                        new Arguments.StringOption(SLICE_FILE_HEIGHTS, "heights_file", "specifies a file with a slice heights-list, is overwritten by command-line specification of slice heights [default=none]"),
                        new Arguments.StringOption(SLICE_FILE_TIMES, "Times_file", "specifies a file with a slice Times-list, is overwritten by command-line specification of slice times [default=none]"),
                        new Arguments.StringOption(STDEVS, falseTrue, true,
                                "get standard deviations for the rates [default = true])"),
                        new Arguments.RealOption(CUTOFF, "specifies the posterior cut-off for summarize grid values [default=-0]"),
                        new Arguments.StringOption(DISCRETE_TRAIT_NAME, "discrete_traitName", "specifies the name for a discrete trait that is annotated to the tree nodes and by which states the grid needs to be summarized [default=none]"),
                        new Arguments.StringOption(DISCRETE_TRAIT_STATES, "discrete_traitStates", "specifies the state of a discrete trait by which the grid needs to be summarized [default=none]"),
                        new Arguments.StringOption(HISTORY_ANNOTATION, "history_annotation", "specifies the name for the history annotation for the discrete trait [default=history]"),

                }
        );

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

        try {

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

            if (arguments.hasOption(MAXPATHLENGTH)) {
                maxPathLength = arguments.getRealOption(MAXPATHLENGTH);
            }

            if (arguments.hasOption(GRIDXCELLS)) {
                gridXcells = arguments.getIntegerOption(GRIDXCELLS);
            }

            if (arguments.hasOption(GRIDYCELLS)) {
                gridYcells = arguments.getIntegerOption(GRIDYCELLS);
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
                System.err.println("Skipping every " + skipEvery + " is not possible, no trees will be skipped");
                skipEvery = 1;
            }

            String traitString = arguments.getStringOption(TRAIT);
            if (traitString != null) {
                traitName = traitString;
            }

            String trait2DnoiseString = arguments.getStringOption(TRAIT_NOISE);
            if (trait2DnoiseString != null && trait2DnoiseString.compareToIgnoreCase("false") == 0)
                trait2DNoise = false;

            String rate2DNoiseString = arguments.getStringOption(TRAIT_NOISE);
            if (rate2DNoiseString != null && rate2DNoiseString.compareToIgnoreCase("true") == 0)
                rate2DNoise = true;

            String normalizeString = arguments.getStringOption(NORMALIZATION);
            if (normalizeString != null) {
                try {
                    normalize = Normalization.valueOf(normalizeString.toUpperCase());
                } catch (IllegalArgumentException iae) {
                    System.err.println("Unrecognized normalization mode: " + normalizeString);
                }

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

            String stdevString = arguments.getStringOption(STDEVS);
            if (stdevString != null && stdevString.compareToIgnoreCase("false") == 0)
                getStdevs = false;

            if (arguments.hasOption(CUTOFF)) {
                posteriorCutoff = arguments.getRealOption(CUTOFF);
            }

            String discreteTraitNameString = arguments.getStringOption(DISCRETE_TRAIT_NAME);
            if (discreteTraitNameString != null) {
                discreteTraitName = discreteTraitNameString;
                if (!arguments.hasOption(DISCRETE_TRAIT_STATES)) {
                    System.err.print("a discrete trait name is specified (" + discreteTraitName + "), but no associated states of that trait are specified");
                    System.exit(-1);
                }
            }

            String discreteTraitStatesString = arguments.getStringOption(DISCRETE_TRAIT_STATES);
            if (discreteTraitStatesString != null) {
                discreteTraitStates = parseVariableLengthStringArray(discreteTraitStatesString);
                if (discreteTraitName == null) {
                    System.err.print("states for a discrete trait are specified, but no associated trait name for trait annotation is specified");
                    System.exit(-1);
                }
            }

            String historyString = arguments.getStringOption(HISTORY_ANNOTATION);
            if (historyString != null) {
                historyAnnotation = historyString;
            }


            //slicing

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
                double startTime;
                if (arguments.hasOption(START_TIME)) {
                    startTime = arguments.getRealOption(START_TIME);
                } else {
                    progressStream.println("slice count specified, but no associated start time?");
                    System.exit(-1);
                    startTime = 0;
                }
                double delta;
                if (mrsd != 0) {
                    if (sliceCount == 1) {
                        delta = (mrsd - startTime);
                    } else {
                        delta = (mrsd - startTime) / (sliceCount - 1);
                    }
                } else {
                    if (sliceCount == 1) {
                        delta = startTime;
                    } else {
                        delta = startTime / (sliceCount - 1);
                    }
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

            if ((sliceTimeString == null) && (sliceHeightString == null) && (!arguments.hasOption(SLICE_COUNT))) {
                sliceHeights = new double[]{0};
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

        Branch2dRateToGrid grid = new Branch2dRateToGrid(inputFileName, burnin, skipEvery, traitName, trait2DNoise, rate2DNoise, rateString, maxPathLength, normalize, latMin, latMax, longMin, longMax, gridXcells, gridYcells, sliceHeights, getStdevs, posteriorCutoff, discreteTraitName, discreteTraitStates, historyAnnotation);
        grid.output(outputFileName);

        System.exit(0);

    }

    private class History {

        private double[] historyHeights;
        private String[] historyStates;
        private double[][] historyTraits;

        public History(double historyHeights[], String historyStates[]) {
            this.historyHeights = historyHeights;
            this.historyStates = historyStates;
        }

        public String getStateForHeight(double height) {
            String returnState = null;
            if ((height > historyHeights[0]) || (height < historyHeights[(historyHeights.length - 1)])) {
                System.err.print("height " + height + " is outside the history range!!");
                System.exit(-1);
            }
            for (int a = 0; a < (historyHeights.length - 1); a++) {
                if ((height <= historyHeights[a]) && (height >= historyHeights[a + 1])) {
                    returnState = historyStates[a];
                }
            }

            return returnState;
        }

        public double[] getHeights() {
            return historyHeights;
        }

        public void truncateUpper(double time) {
            int cutFrom = -1;
            for (int a = 0; a < (historyHeights.length - 1); a++) {
                if ((time < historyHeights[a]) && (time > historyHeights[a + 1])) {
                    cutFrom = a;
                }
            }

            if (cutFrom < 0) {
                System.err.println("no upper truncation of discrete trait history on branch possible");
                System.exit(0);
            }

            double[] tempHeights = new double[historyHeights.length - cutFrom];
            String[] tempStates = new String[historyStates.length - cutFrom];

            tempHeights = Arrays.copyOfRange(historyHeights, cutFrom, historyHeights.length);
            tempHeights[0] = time;

            tempStates = Arrays.copyOfRange(historyStates, cutFrom, historyStates.length);

            historyHeights = tempHeights;
            historyStates = tempStates;
        }

        public void truncateLower(double time) {
            int cutTo = -1;

            for (int a = (historyHeights.length - 1); a > 0; a--) {
                if ((time > historyHeights[a]) && (time < historyHeights[a - 1])) {
                    cutTo = a;
                }
            }

            if (cutTo < 0) {
                System.err.println("no lower truncation of discrete trait history on branch possible");
                System.exit(0);
            }

            double[] tempHeights = new double[cutTo + 1];
            String[] tempStates = new String[cutTo];

            tempHeights = Arrays.copyOfRange(historyHeights, 0, cutTo + 1);
            tempHeights[(tempHeights.length - 1)] = time;

            tempStates = Arrays.copyOfRange(historyStates, 0, cutTo);

            historyHeights = tempHeights;
            historyStates = tempStates;
        }

        public double getStateTime(String state) {
            double time = 0;
            for (int x = 0; x < historyStates.length; x++) {
                if (state.equals(historyStates[x])) {
                    time += (historyHeights[x] - historyHeights[x + 1]);
                }
            }
            return time;
        }

        private void setTraitsforHeights(double[] traitUp, double[] traitLow, double[] precisionArray, double rate, boolean trueNoise) {
            historyTraits = new double[historyHeights.length][2];
            for (int x = 0; x < historyHeights.length; x++) {
                if (x == 0) {
                    historyTraits[x] = traitUp;
                } else if (x == (historyTraits.length - 1)) {
                    historyTraits[x] = traitLow;
                } else {
                    historyTraits[x] = imputeValue(traitUp, traitLow, historyHeights[x], historyHeights[(historyHeights.length - 1)], historyHeights[0], precisionArray, rate, trueNoise);
                }

            }
        }

        public double getStateGreatCircleDistance(String state) {
            double distance = 0;
            for (int x = 0; x < historyStates.length; x++) {
                if (state.equals(historyStates[x])) {
                    distance += getKilometerGreatCircleDistance(historyTraits[x], historyTraits[x + 1]);
                }
            }
            return distance;
        }

        public double getStateDifferenceInGreatCircleDistanceFromRoot(String state, double[] rootTrait) {
            double distance = 0;
            for (int x = 0; x < historyStates.length; x++) {
                if (state.equals(historyStates[x])) {
                    distance += (getKilometerGreatCircleDistance(historyTraits[x + 1], rootTrait) - getKilometerGreatCircleDistance(historyTraits[x], rootTrait));
                }
            }
            return distance;
        }

        public double getStateNativeDistance(String state) {
            double distance = 0;
            for (int x = 0; x < historyStates.length; x++) {
                if (state.equals(historyStates[x])) {
                    distance += getNativeDistance(historyTraits[x], historyTraits[x + 1]);
                }
            }
            return distance;
        }

        public double getStateDifferenceInNativeDistanceFromRoot(String state, double[] rootTrait) {
            double distance = 0;
            for (int x = 0; x < historyStates.length; x++) {
                if (state.equals(historyStates[x])) {
                    distance += (getNativeDistance(historyTraits[x + 1], rootTrait) - getNativeDistance(historyTraits[x], rootTrait));
                }
            }
            return distance;
        }

    }

    class Trait {

        Trait(Object obj) {
            this.obj = obj;
            if (obj instanceof Object[]) {
                isMultivariate = true;
                array = (Object[]) obj;
            }
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
            for (int i = 0; i < dim; i++)
                result[i] = (Double) array[i];
            return result;
        }

        private Object obj;
        private Object[] array;
        private boolean isMultivariate = false;

        public String toString() {
            if (!isMultivariate)
                return obj.toString();
            StringBuffer sb = new StringBuffer(array[0].toString());
            for (int i = 1; i < array.length; i++)
                sb.append("\t").append(array[i]);
            return sb.toString();
        }
    }


}

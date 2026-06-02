/*
 * AnnotationAction.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

package dr.app.tools.treeannotator;

import dr.evolution.tree.MutableTree;
import dr.evolution.tree.NodeRef;
import dr.evolution.tree.Tree;
import dr.geo.contouring.ContourMaker;
import dr.geo.contouring.ContourPath;
import dr.geo.contouring.ContourWithSynder;
import dr.stats.DiscreteStatistics;
import dr.util.HeapSort;

import java.util.*;

public class AnnotationAction implements CladeAction {
    private static final boolean USE_R = false;

    private static final List<TreeAnnotationPlugin> plugins = new ArrayList<>();

    private final TreeAnnotator.HeightsSummary heightsOption;
    private double posteriorLimit;
    private double countLimit;
    double[] hpd2D = {0.80};
    Set<String> attributeNames = new HashSet<String>();
    private boolean forceIntegerToDiscrete = false;
    private boolean computeESS;

    private final String location1Attribute = "longLat1";
    private final String location2Attribute = "longLat2";
    private final String locationOutputAttribute = "location";

    private final static boolean PROCESS_BIVARIATE_ATTRIBUTES = true;

    AnnotationAction(TreeAnnotator.HeightsSummary heightsOption,
                     double posteriorLimit,
                     int countLimit,
                     double[] hpd2D,
                     boolean computeESS,
                     boolean forceIntegerToDiscrete) {
        this.heightsOption = heightsOption;
        this.posteriorLimit = posteriorLimit;
        this.countLimit = countLimit;
        this.hpd2D = hpd2D;
        this.forceIntegerToDiscrete = forceIntegerToDiscrete;
    }

    public void addAttributeName(String attributeName) {
        this.attributeNames.add(attributeName);
    }
    public void addAttributeNames(Collection<String> attributeNames) {
        this.attributeNames.addAll(attributeNames);
    }


    @Override
    public void actOnClade(Clade clade, Tree tree, NodeRef node) {
        assert tree instanceof MutableTree;
        annotateNode((MutableTree)tree, node, clade);
    }

    @Override
    public boolean expectAllClades() {
        return false;
    }

    private void annotateNode(MutableTree tree, NodeRef node, Clade clade) {
        boolean filter = false;
        assert clade != null;

        if (!tree.isExternal(node)) {
            final double posterior = clade.getCredibility();
            tree.setNodeAttribute(node, "posterior", posterior);
            if (posterior < posteriorLimit || clade.getCount() < countLimit) {
                filter = true;
            }
        }

        int i = 0;
        for (String attributeName : attributeNames) {
            if (attributeName.equals("height")) {
                if (!filter) {
                    tree.setNodeAttribute(node, "height_mean", ((BiClade) clade).getMeanHeight());
                    tree.setNodeAttribute(node, "height_median", ((BiClade) clade).getMedianHeight());
                    if (((BiClade) clade).getHeightHPDs() != null){
                        tree.setNodeAttribute(node, "height_95%_HPD", ((BiClade) clade).getHeightHPDs());
                    }
                    if (((BiClade) clade).getHeightRange() != null){
                        tree.setNodeAttribute(node, "height_range", ((BiClade) clade).getHeightRange());
                    }
                }
                if (heightsOption == TreeAnnotator.HeightsSummary.MEAN_HEIGHTS) {
                    tree.setNodeHeight(node, ((BiClade) clade).getMeanHeight());
                } else if (heightsOption == TreeAnnotator.HeightsSummary.MEDIAN_HEIGHTS) {
                    tree.setNodeHeight(node, ((BiClade) clade).getMedianHeight());
                } else {
                    // keep the existing height
                }
//                assert tree.isExternal(node) || (tree.getNodeHeight(node) - tree.getNodeHeight(tree.getChild(node, 0))) >= 0.0;
//                assert tree.isExternal(node) || (tree.getNodeHeight(node) - tree.getNodeHeight(tree.getChild(node, 1))) >= 0.0;
            } else {
                if (clade.getAttributeValues() != null && !clade.getAttributeValues().isEmpty()) {
                    double[] values = new double[clade.getAttributeValues().size()];

                    HashMap<Object, Integer> hashMap = new HashMap<>();

                    Object[] v = clade.getAttributeValues().get(0);
                    if (v[i] != null) {
                        boolean isBoolean = v[i] instanceof Boolean;

                        boolean isDiscrete = v[i] instanceof String;

                        if (forceIntegerToDiscrete && v[i] instanceof Integer) isDiscrete = true;

                        double minValue = Double.MAX_VALUE;
                        double maxValue = -Double.MAX_VALUE;

                        final boolean isArray = v[i] instanceof Object[];
                        boolean isDoubleArray = isArray && ((Object[]) v[i])[0] instanceof Double;
                        // This is Java, friends - first value type does not imply all.
                        if (isDoubleArray) {
                            for (Object n : (Object[]) v[i]) {
                                if (!(n instanceof Double)) {
                                    isDoubleArray = false;
                                    break;
                                }
                            }
                        }
                        // todo Handle other types of arrays

                        double[][] valuesArray = null;
                        double[] minValueArray = null;
                        double[] maxValueArray = null;
                        int lenArray = 0;

                        if (isDoubleArray) {
                            lenArray = ((Object[]) v[i]).length;

                            valuesArray = new double[lenArray][clade.getAttributeValues().size()];
                            minValueArray = new double[lenArray];
                            maxValueArray = new double[lenArray];

                            for (int k = 0; k < lenArray; k++) {
                                minValueArray[k] = Double.MAX_VALUE;
                                maxValueArray[k] = -Double.MAX_VALUE;
                            }
                        }

                        for (int j = 0; j < clade.getAttributeValues().size(); j++) {
                            Object value = clade.getAttributeValues().get(j)[i];
                            if (isDiscrete) {
                                final Object s = value;
                                if (hashMap.containsKey(s)) {
                                    hashMap.put(s, hashMap.get(s) + 1);
                                } else {
                                    hashMap.put(s, 1);
                                }
                            } else if (isBoolean) {
                                values[j] = (((Boolean) value) ? 1.0 : 0.0);
                            } else if (isDoubleArray) {
                                // Forcing to Double[] causes a cast exception. MAS
                                Object[] array = (Object[]) value;
                                for (int k = 0; k < lenArray; k++) {
                                    valuesArray[k][j] = ((Double) array[k]);
                                    if (valuesArray[k][j] < minValueArray[k]) minValueArray[k] = valuesArray[k][j];
                                    if (valuesArray[k][j] > maxValueArray[k]) maxValueArray[k] = valuesArray[k][j];
                                }
                            } else {
                                // Ignore other (unknown) types
                                if (value instanceof Number) {
                                    values[j] = ((Number) value).doubleValue();
                                    if (values[j] < minValue) minValue = values[j];
                                    if (values[j] > maxValue) maxValue = values[j];
                                }
                            }
                        }
                        if (!filter) {
                            boolean processed = false;
                            for (TreeAnnotationPlugin plugin : plugins) {
                                if (plugin.handleAttribute(tree, node, attributeName, values)) {
                                    processed = true;
                                }
                            }

                            if (!processed) {
                                if (!isDiscrete) {
                                    if (!isDoubleArray)
                                        annotateMeanAttribute(tree, node, attributeName, values);
                                    else {
                                        for (int k = 0; k < lenArray; k++) {
                                            annotateMeanAttribute(tree, node, attributeName + (k + 1), valuesArray[k]);
                                        }
                                    }
                                } else {
                                    annotateModeAttribute(tree, node, attributeName, hashMap);
                                    annotateFrequencyAttribute(tree, node, attributeName, hashMap);
                                }
                                if (!isBoolean && minValue < maxValue && !isDiscrete && !isDoubleArray) {
                                    // Basically, if it is a boolean (0, 1) then we don't need the distribution information
                                    // Likewise if it doesn't vary.
                                    annotateMedianAttribute(tree, node, attributeName + "_median", values);
                                    annotateHPDAttribute(tree, node, attributeName + "_95%_HPD", 0.95, values);
                                    annotateRangeAttribute(tree, node, attributeName + "_range", values);
                                    annotateSignAttribute(tree, node, attributeName + "_signDistribution", values);
                                    if (computeESS == true) {
                                        annotateESSAttribute(tree, node, attributeName + "_ESS", values);
                                    }
                                }

                                if (isDoubleArray) {
                                    String name = attributeName;
                                    // todo
//                                    if (name.equals(location1Attribute)) {
//                                        name = locationOutputAttribute;
//                                    }
                                    boolean want2d = PROCESS_BIVARIATE_ATTRIBUTES && lenArray == 2;
                                    if (name.equals("dmv")) {  // terrible hack
                                        want2d = false;
                                    }
                                    for (int k = 0; k < lenArray; k++) {
                                        if (minValueArray[k] < maxValueArray[k]) {
                                            annotateMedianAttribute(tree, node, name + (k + 1) + "_median", valuesArray[k]);
                                            annotateRangeAttribute(tree, node, name + (k + 1) + "_range", valuesArray[k]);
                                            annotatePositiveProbability(tree, node, name + (k + 1) + "_positiveProb", valuesArray[k]);
                                            if (!want2d)
                                                annotateHPDAttribute(tree, node, name + (k + 1) + "_95%_HPD", 0.95, valuesArray[k]);
                                        }
                                    }
                                    // 2D contours
                                    if (want2d) {

                                        boolean variationInFirst = (minValueArray[0] < maxValueArray[0]);
                                        boolean variationInSecond = (minValueArray[1] < maxValueArray[1]);

                                        if (variationInFirst && !variationInSecond)
                                            annotateHPDAttribute(tree, node, name + "1" + "_95%_HPD", 0.95, valuesArray[0]);

                                        if (variationInSecond && !variationInFirst)
                                            annotateHPDAttribute(tree, node, name + "2" + "_95%_HPD", 0.95, valuesArray[1]);

                                        if (variationInFirst && variationInSecond) {

                                            for (int l = 0; l < hpd2D.length; l++) {

                                                if (hpd2D[l] > 1) {
                                                    System.err.println("no HPD for proportion > 1 (" + hpd2D[l] + ")");
                                                } else if (hpd2D[l] < 0) {
                                                    System.err.println("no HPD for proportion < 0 (" + hpd2D[l] + ")");
                                                } else {
                                                    annotate2DHPDAttribute(tree, node, name, "_" + (int) (100 * hpd2D[l]) + "%HPD", hpd2D[l], valuesArray);
                                                }

                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                i++;
            }
        }
    }

    private void annotateMeanAttribute(MutableTree tree, NodeRef node, String label, double[] values) {
        double mean = DiscreteStatistics.mean(values);
        tree.setNodeAttribute(node, label, mean);
    }

    private void annotateMedianAttribute(MutableTree tree, NodeRef node, String label, double[] values) {
        double median = DiscreteStatistics.median(values);
        tree.setNodeAttribute(node, label, median);
    }

    private void annotateModeAttribute(MutableTree tree, NodeRef node, String label, HashMap<Object, Integer> values) {
        Object mode = null;
        int maxCount = 0;
        int totalCount = 0;
        int countInMode = 1;

        for (Object key : values.keySet()) {
            int thisCount = values.get(key);
            if (thisCount == maxCount) {
                // I hope this is the intention
                mode = mode.toString().concat("+" + key);
                countInMode++;
            } else if (thisCount > maxCount) {
                mode = key;
                maxCount = thisCount;
                countInMode = 1;
            }
            totalCount += thisCount;
        }
        double freq = (double) maxCount / (double) totalCount * countInMode;
        tree.setNodeAttribute(node, label, mode);
        tree.setNodeAttribute(node, label + ".prob", freq);
    }

    private void annotateFrequencyAttribute(MutableTree tree, NodeRef node, String label, HashMap<Object, Integer> values) {
        double totalCount = 0;
        Set keySet = values.keySet();
        int length = keySet.size();
        String[] name = new String[length];
        Double[] freq = new Double[length];
        int index = 0;
        for (Object key : values.keySet()) {
            name[index] = key.toString();
            freq[index] = (double) values.get(key);
            totalCount += freq[index];
            index++;
        }
        for (int i = 0; i < length; i++)
            freq[i] /= totalCount;

        tree.setNodeAttribute(node, label + ".set", name);
        tree.setNodeAttribute(node, label + ".set.prob", freq);
    }

    private void annotateSignAttribute(MutableTree tree, NodeRef node, String label, double[] values) {
        double negativePortion = DiscreteStatistics.negativeProbability(values);
        double positivePortion = 1 - negativePortion;
        tree.setNodeAttribute(node, label, new Object[]{negativePortion, positivePortion});
    }

    private void annotatePositiveProbability(MutableTree tree, NodeRef node, String label, double[] values) {
        double negativePortion = DiscreteStatistics.negativeProbability(values);
        double positivePortion = 1 - negativePortion;
        tree.setNodeAttribute(node, label, positivePortion);
    }

    private void annotateRangeAttribute(MutableTree tree, NodeRef node, String label, double[] values) {
        double min = DiscreteStatistics.min(values);
        double max = DiscreteStatistics.max(values);
        tree.setNodeAttribute(node, label, new Object[]{min, max});
    }

    private void annotateHPDAttribute(MutableTree tree, NodeRef node, String label, double hpd, double[] values) {
        int[] indices = new int[values.length];
        HeapSort.sort(values, indices);

        double minRange = Double.MAX_VALUE;
        int hpdIndex = 0;

        int diff = (int) Math.round(hpd * (double) values.length);
        for (int i = 0; i <= (values.length - diff); i++) {
            double minValue = values[indices[i]];
            double maxValue = values[indices[i + diff - 1]];
            double range = Math.abs(maxValue - minValue);
            if (range < minRange) {
                minRange = range;
                hpdIndex = i;
            }
        }
        double lower = values[indices[hpdIndex]];
        double upper = values[indices[hpdIndex + diff - 1]];
        tree.setNodeAttribute(node, label, new Object[]{lower, upper});
    }

    private void annotateESSAttribute(MutableTree tree, NodeRef node, String label, double[] values) {
        throw new UnsupportedOperationException("Not fully implemented");
        // array --> list (to construct traceCorrelation obj)
//        List<Double> values2 = new ArrayList<Double>(0);
//        for (int i = 0; i < values.length; i++) {
//            values2.add(values[i]);
//        }
//
//        TraceType traceType = TraceType.REAL;
//        // maxState / totalTrees = stepSize for ESS
//        int logStep = (int) (maxState / totalTrees);
//        TraceCorrelation traceCorrelation = new TraceCorrelation(values2, traceType, logStep);
//
//        double ESS = traceCorrelation.getESS();
//        tree.setNodeAttribute(node, label, ESS);
    }

    // todo Move rEngine to outer class; create once.
//        Rengine rEngine = null;

    private final String[] rArgs = {"--no-save"};

//	    private int called = 0;

    private final String[] rBootCommands = {
            "library(MASS)",
            "makeContour = function(var1, var2, prob=0.95, n=50, h=c(1,1)) {" +
                    "post1 = kde2d(var1, var2, n = n, h=h); " +    // This had h=h in argument
                    "dx = diff(post1$x[1:2]); " +
                    "dy = diff(post1$y[1:2]); " +
                    "sz = sort(post1$z); " +
                    "c1 = cumsum(sz) * dx * dy; " +
                    "levels = sapply(prob, function(x) { approx(c1, sz, xout = 1 - x)$y }); " +
                    "line = contourLines(post1$x, post1$y, post1$z, level = levels); " +
                    "return(line) }"
    };

    private String makeRString(double[] values) {
        StringBuilder sb = new StringBuilder("c(");
        sb.append(values[0]);
        for (int i = 1; i < values.length; i++) {
            sb.append(",");
            sb.append(values[i]);
        }
        sb.append(")");
        return sb.toString();
    }

    public static final String CORDINATE = "cordinates";

//		private String formattedLocation(double loc1, double loc2) {
//			return formattedLocation(loc1) + "," + formattedLocation(loc2);
//		}

    private String formattedLocation(double x) {
        return String.format("%5.8f", x);
    }

    private void annotate2DHPDAttribute(MutableTree tree, NodeRef node, String preLabel, String postLabel,
                                        double hpd, double[][] values) {
        int N = 50;
        if (USE_R) {
//
//                // Uses R-Java interface, and the HPD routines from 'emdbook' and 'coda'
//
//                if (rEngine == null) {
//
//                    if (!Rengine.versionCheck()) {
//                        throw new RuntimeException("JRI library version mismatch");
//                    }
//
//                    rEngine = new Rengine(rArgs, false, null);
//
//                    if (!rEngine.waitForR()) {
//                        throw new RuntimeException("Cannot load R");
//                    }
//
//                    for (String command : rBootCommands) {
//                        rEngine.eval(command);
//                    }
//                }
//
//                // todo Need a good method to pick grid size
//
//
//                REXP x = rEngine.eval("makeContour(" +
//                        makeRString(values[0]) + "," +
//                        makeRString(values[1]) + "," +
//                        hpd + "," +
//                        N + ")");
//
//                RVector contourList = x.asVector();
//                int numberContours = contourList.size();
//
//                if (numberContours > 1) {
//                    System.err.println("Warning: a node has a disjoint " + 100 * hpd + "% HPD region.  This may be an artifact!");
//                    System.err.println("Try decreasing the enclosed mass or increasing the number of samples.");
//                }
//
//
//                tree.setNodeAttribute(node, preLabel + postLabel + "_modality", numberContours);
//
//                StringBuffer output = new StringBuffer();
//                for (int i = 0; i < numberContours; i++) {
//                    output.append("\n<" + CORDINATE + ">\n");
//                    RVector oneContour = contourList.at(i).asVector();
//                    double[] xList = oneContour.at(1).asDoubleArray();
//                    double[] yList = oneContour.at(2).asDoubleArray();
//                    StringBuffer xString = new StringBuffer("{");
//                    StringBuffer yString = new StringBuffer("{");
//                    for (int k = 0; k < xList.length; k++) {
//                        xString.append(formattedLocation(xList[k])).append(",");
//                        yString.append(formattedLocation(yList[k])).append(",");
//                    }
//                    xString.append(formattedLocation(xList[0])).append("}");
//                    yString.append(formattedLocation(yList[0])).append("}");
//
//                    tree.setNodeAttribute(node, preLabel + "1" + postLabel + "_" + (i + 1), xString);
//                    tree.setNodeAttribute(node, preLabel + "2" + postLabel + "_" + (i + 1), yString);
//                }
//
//
        } else { // do not use R


//                KernelDensityEstimator2D kde = new KernelDensityEstimator2D(values[0], values[1], N);
            //ContourMaker kde = new ContourWithSynder(values[0], values[1], N);
            boolean bandwidthLimit = false;

            ContourMaker kde = new ContourWithSynder(values[0], values[1], bandwidthLimit);

            ContourPath[] paths = kde.getContourPaths(hpd);

            tree.setNodeAttribute(node, preLabel + postLabel + "_modality", paths.length);

            if (paths.length > 1) {
                System.err.println("Warning: a node has a disjoint " + 100 * hpd + "% HPD region.  This may be an artifact!");
                System.err.println("Try decreasing the enclosed mass or increasing the number of samples.");
            }

            StringBuilder output = new StringBuilder();
            int i = 0;
            for (ContourPath p : paths) {
                output.append("\n<").append(CORDINATE).append(">\n");
                double[] xList = p.getAllX();
                double[] yList = p.getAllY();
                StringBuilder xString = new StringBuilder("{");
                StringBuilder yString = new StringBuilder("{");
                for (int k = 0; k < xList.length; k++) {
                    xString.append(formattedLocation(xList[k])).append(",");
                    yString.append(formattedLocation(yList[k])).append(",");
                }
                xString.append(formattedLocation(xList[0])).append("}");
                yString.append(formattedLocation(yList[0])).append("}");

                tree.setNodeAttribute(node, preLabel + "1" + postLabel + "_" + (i + 1), xString);
                tree.setNodeAttribute(node, preLabel + "2" + postLabel + "_" + (i + 1), yString);
                i++;

            }
        }
    }

    public interface TreeAnnotationPlugin {
        Set<String> setAttributeNames(Set<String> attributeNames);

        boolean handleAttribute(Tree tree, NodeRef node, String attributeName, double[] values);
    }

}

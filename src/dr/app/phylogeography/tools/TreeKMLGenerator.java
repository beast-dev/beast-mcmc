/*
 * TreeKMLGenerator.java
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

import dr.app.phylogeography.tools.kml.*;
import jebl.evolution.trees.RootedTree;
import jebl.evolution.graphs.Node;
import jebl.evolution.io.TreeImporter;
import jebl.evolution.io.NexusImporter;
import jebl.evolution.io.ImportException;
import org.jdom.*;
import org.jdom.output.XMLOutputter;
import org.jdom.output.Format;

import java.util.*;
import java.util.List;
import java.awt.*;
import java.io.FileReader;
import java.io.IOException;
import java.io.BufferedWriter;
import java.io.FileWriter;

/**
 * @author Philippe Lemey
 * @author Andrew Rambaut
 * @author Marc A. Suchard
 */
public class TreeKMLGenerator {

    // used to calibrate the color range for the branches
    double minRate;
    double maxRate;
    double medianRate;

    double maxHeight;
    double maxBranchLength;

    String altitudeMode = "clampToGround";

    private final RootedTree tree;
    private final Map<String, Location> locationMap;


    public TreeKMLGenerator(RootedTree tree) {
        this(tree, null);
    }

    public TreeKMLGenerator(RootedTree tree, List<Location> locationList) {
        this.tree = tree;

//        double[] rates = new double[(tree.getNodes().size() - 1)];
//
//        maxBranchLength = 0.0;
//        int counter = 0;
//        for (Node node : tree.getNodes()) {
//            if (!tree.isRoot(node)) {
//                rates[counter] = (Double)node.getAttribute("rate");
//                counter ++;
//
//                if (tree.getLength(node) > maxBranchLength) {
//                    maxBranchLength = tree.getLength(node);
//                }
//            }
//        }
//
//        minRate = DiscreteStatistics.min(rates);
//        maxRate = DiscreteStatistics.max(rates);
//        medianRate = DiscreteStatistics.quantile(0.5, rates);

        maxHeight = tree.getHeight(tree.getRootNode());
        maxBranchLength = 0.0;
        for (Node node : tree.getNodes()) {
            if (!tree.isRoot(node)) {
                if (tree.getLength(node) > maxBranchLength) {
                    maxBranchLength = tree.getLength(node);
                }
            }
        }

        if (locationList != null) {
            double minLat = Double.MAX_VALUE;
            double maxLat = -Double.MAX_VALUE;
            double minLong = Double.MAX_VALUE;
            double maxLong = -Double.MAX_VALUE;

            locationMap = new HashMap<String, Location>();

            for (Location location : locationList) {
                if (location.getLatitude() < minLat)  {
                    minLat = location.getLatitude();
                }
                if (location.getLatitude() > maxLat)  {
                    maxLat = location.getLatitude();
                }

                if (location.getLongitude() < minLong)  {
                    minLong = location.getLongitude();
                }
                if (location.getLongitude() > maxLong)  {
                    maxLong = location.getLongitude();
                }
                locationMap.put(location.getState(), location);
            }

//            radius = 100*Math.abs(maxLat-minLat)*Math.abs(maxLong-minLong);
//            radius = 200000;

//            if (mostRecentDate - treeToExport.getHeight(treeToExport.getRootNode()) < 0) {
//                ancient = true;
//            }

        } else {
            locationMap = null;
        }
    }

    public Element generate(String documentName, Settings settings) {

        Element root = new Element("kml");

        Element doc = new Element("Document");
        doc.addContent(generateElement("name", documentName));

        List<Element> schema = new ArrayList<Element>();
        Element branchSchema = new Element("Schema");
        branchSchema.setAttribute("id", "Branch_Schema");
        branchSchema.addContent(new Element("SimpleField")
                .setAttribute("name", "Name")
                .setAttribute("type", "string")
                .addContent(new Element("displayName").addContent("Name")));
        branchSchema.addContent(new Element("SimpleField")
                .setAttribute("name", "Height")
                .setAttribute("type", "double")
                .addContent(new Element("displayName").addContent("Height")));
        branchSchema.addContent(new Element("SimpleField")
                .setAttribute("name", "StartTime")
                .setAttribute("type", "double")
                .addContent(new Element("displayName").addContent("StartTime")));
        branchSchema.addContent(new Element("SimpleField")
                .setAttribute("name", "FinishTime")
                .setAttribute("type", "double")
                .addContent(new Element("displayName").addContent("FinishTime")));
        branchSchema.addContent(new Element("SimpleField")
                .setAttribute("name", "Support")
                .setAttribute("type", "double")
                .addContent(new Element("displayName").addContent("Support")));
        branchSchema.addContent(new Element("SimpleField")
                .setAttribute("name", "Rate")
                .setAttribute("type", "double")
                .addContent(new Element("displayName").addContent("Rate")));
        schema.add(branchSchema);

        Element nodeSchema = new Element("Schema");
        branchSchema.setAttribute("id", "Node_Schema");
        branchSchema.addContent(new Element("SimpleField")
                .setAttribute("name", "Label")
                .setAttribute("type", "string")
                .addContent(new Element("displayName").addContent("Label")));
        branchSchema.addContent(new Element("SimpleField")
                .setAttribute("name", "Height")
                .setAttribute("type", "double")
                .addContent(new Element("displayName").addContent("Height")));
        branchSchema.addContent(new Element("SimpleField")
                .setAttribute("name", "Time")
                .setAttribute("type", "double")
                .addContent(new Element("displayName").addContent("Time")));
        schema.add(nodeSchema);

        List<Element> styles = new ArrayList<Element>();
        List<Element> trees = new ArrayList<Element>();
        List<Element> projections = new ArrayList<Element>();
        List<Element> taxonLabels = new ArrayList<Element>();
        List<Element> locationLabels = new ArrayList<Element>();
        List<Element> contours = new ArrayList<Element>();

//        styles.add(generateLineStyle("projections", 1.0, Color.lightGray));
//        styles.add(generatePolyStyle("contours", new Color(64, 255, 0, 128)));

        if (settings.getAnalysisType() == AnalysisType.CONTINUOUS) {
//            trees.add(generateTree(tree, settings, settings.getAltitudeTreeSettings(),  styles));
            trees.add(generateTree(tree, settings, settings.getGroundTreeSettings(),  styles));

//            double scaleFactor = settings.getPlotAltitude()/tree.getHeight(tree.getRootNode());
//
//            int nodeNumber = 0;
//            for (Node node : tree.getNodes()) {
//                nodeNumber++;
//
//                double latitude = getDoubleNodeAttribute(node, settings.getLatitudeName());
//                double longitude = getDoubleNodeAttribute(node, settings.getLongitudeName());
//                double altitude = (tree.getHeight(node)*scaleFactor);
//
//
//                if (tree.isExternal(node)) {
//                    projections.add(generateProjection(altitude, latitude, longitude));
//
//                    taxonLabels.add(generateTaxonLabel(tree, node, latitude, longitude));
//                } else {
//                    double date = settings.getMostRecentDate() - tree.getHeight(node);
//                    if (settings.getAgeCutOff() > 0.0 && date > settings.getAgeCutOff()) {
//
//                        contours.add(generateContour(
//                                settings.getGroundContours(), tree, node,  nodeNumber, 0.0, date,
//                                settings.getTraitName(),  settings.getLatitudeName(), settings.getLongitudeName(),
//                                true, styles));
//                    }
//                }
//            }

        } else if (settings.getAnalysisType() == AnalysisType.DISCRETE) {

            for (Location location : locationMap.values()) {
                locationLabels.add(generateLocation(location));
            }

        } else {
            throw new IllegalArgumentException("Unknown AnalysisType");
        }

        doc.addContent(schema);
//        doc.addContent(styles);
        doc.addContent(trees);

        if (projections.size() > 0) {
            Element placeMark = generateContainer("Placemark", "Projections", "projections from tips to surface", true, "projections");
            placeMark.addContent(new Element("MultiGeometry").addContent(projections));
            doc.addContent(placeMark);
        }

        if (taxonLabels.size() > 0) {
            Element folder = generateContainer("Folder", "Taxon Labels", "Taxon labels", settings.getTaxonLabels().isVisible(), null);
            folder.addContent(taxonLabels);
//        placeMark.addContent(new Element("MultiGeometry").addContent(labels));
            doc.addContent(folder);
        }

        if (locationLabels.size() > 0) {
            Element folder = generateContainer("Folder", "Location Labels", "Location labels", settings.getLocationLabels().isVisible(), null);
            folder.addContent(locationLabels);
//        placeMark.addContent(new Element("MultiGeometry").addContent(labels));
            doc.addContent(folder);
        }

        if (contours.size() > 0) {
            Element folder = generateContainer("Folder", "Surface Contours", "Contours overlaid on the surface to represent location", true, null);
            folder.addContent(contours);
            doc.addContent(folder);
        }

        root.addContent(doc);

        return root;
    }

    private Element generateTree(RootedTree tree, Settings settings, TreeSettings treeSettings, List<Element> styles) {

        BranchStyle branchStyle = treeSettings.getBranchStyle();

        Element element = generateContainer("Folder", treeSettings.getName(), treeSettings.getDescription(), true, null);

        double scaleFactor = settings.getPlotAltitude() / tree.getHeight(tree.getRootNode());

        int nodeNumber = 0;
        for (Node node : tree.getNodes()) {
            nodeNumber++;

            double latitude;
            double longitude;
            Object[] location = getArrayNodeAttribute(node, settings.getTraitName());
            if (location != null) {
                latitude =  (Double)location[0];
                longitude =  (Double)location[0];
            } else {
             latitude = getDoubleNodeAttribute(node, settings.getLatitudeName());
             longitude = getDoubleNodeAttribute(node, settings.getLongitudeName());
            }
            double altitude = (tree.getHeight(node)*scaleFactor);
            double date = settings.getMostRecentDate() - tree.getHeight(node);

            if (tree.isExternal(node)) {
                System.out.println(tree.getTaxon(node).getName() + "\t" + latitude + "\t" + longitude);
            }

            if (!tree.isRoot(node)) {

                // Create each branch of the tree..
                Element branch;

                String nodeName = treeSettings.getName() + "_node" + nodeNumber;
                if (tree.isExternal(node)) {
                    Element tip = generateTaxonLabel(tree, node, tree.getTaxon(node).getName(), tree.getHeight(node), date, latitude, longitude);
                    element.addContent(tip);
                }

                Node parentNode = tree.getParent(node);
                double parentLatitude = getDoubleNodeAttribute(parentNode, settings.getLatitudeName());
                double parentLongitude = getDoubleNodeAttribute(parentNode, settings.getLongitudeName());
                double parentAltitude = (tree.getHeight(parentNode)*scaleFactor);
                double parentDate = settings.getMostRecentDate() - tree.getHeight(parentNode);

                double width = branchStyle.getWidthBase();
                if (branchStyle.getWidthProperty() != null) {
                    double property = getDoubleNodeAttribute(node, branchStyle.getWidthProperty(), 0.0);
                    width += ((property - branchStyle.getWidthPropertyMinimum()) /
                            (branchStyle.getWidthPropertyMaximum() - branchStyle.getWidthPropertyMinimum())) * branchStyle.getWidthScale();
                }

                Color color = branchStyle.getColorStart();
                if (branchStyle.getColorProperty() != null) {
                    double property;
                    if (branchStyle.getColorProperty().equalsIgnoreCase("height")) {
                        property = tree.getHeight(node) / maxHeight;
                    } else {
                        property = getDoubleNodeAttribute(node, branchStyle.getColorProperty());
                        property = ((property - branchStyle.getWidthPropertyMinimum()) /
                                (branchStyle.getWidthPropertyMaximum() - branchStyle.getWidthPropertyMinimum())) * branchStyle.getWidthScale();
                    }
                    color = getBlendedColor((float)property, branchStyle.getColorStart(), branchStyle.getColorFinish());
                }

                LineStyle lineStyle = new LineStyle(width, color);

                if (treeSettings.getTreeType() == TreeType.RECTANGLE_TREE ||
                        treeSettings.getTreeType() == TreeType.TRIANGLE_TREE) {

                    branch = generateBranch(treeSettings.getTreeType(), tree, node, nodeName,
                            parentLatitude, latitude,
                            parentLongitude, longitude,
                            parentAltitude, altitude,
                            0.0, parentAltitude,
                            date, date,
                            settings.getAgeCutOff(),
                            lineStyle, null,
                            0,
                            styles);
                } else if (treeSettings.getTreeType() == TreeType.ARC_TREE) {

                    LineStyle lineStyle2 = null;

                    // If we are coloring by height and we are using subdivided branches (ARC_TREE) then we can give the
                    // color range from the parent height to the node height.
                    if (branchStyle.getColorProperty() != null && branchStyle.getColorProperty().equalsIgnoreCase("height")) {
                        double divisor = (settings.getAgeCutOff() > 0 ? settings.getMostRecentDate() - settings.getAgeCutOff(): maxHeight);
                        color = getBlendedColor((float)(tree.getHeight(parentNode) / divisor), branchStyle.getColorStart(), branchStyle.getColorFinish());
                        lineStyle = new LineStyle(width, color);
                        color = getBlendedColor((float)(tree.getHeight(node) / divisor), branchStyle.getColorStart(), branchStyle.getColorFinish());
                        lineStyle2 = new LineStyle(width, color);
                    }

                    branch = generateBranch(treeSettings.getTreeType(), tree, node, nodeName,
                            parentLatitude, latitude,
                            parentLongitude, longitude,
                            parentAltitude, altitude,
                            0.0, parentAltitude,
                            parentDate, date,
                            settings.getAgeCutOff(),
                            lineStyle, lineStyle2,
                            50,
                            styles);
                } else if (treeSettings.getTreeType() == TreeType.SURFACE_TREE) {
                    double maxAltitude = settings.getPlotAltitude() * (tree.getHeight(parentNode) - tree.getHeight(node)) / maxBranchLength;

                    LineStyle lineStyle2 = null;

                    // If we are coloring by height and we are using subdivided branches (ARC_TREE) then we can give the
                    // color range from the parent height to the node height.
                    if (branchStyle.getColorProperty() != null && branchStyle.getColorProperty().equalsIgnoreCase("height")) {
                        double divisor = (settings.getAgeCutOff() > 0 ? settings.getMostRecentDate() - settings.getAgeCutOff(): maxHeight);
                        color = getBlendedColor((float)(tree.getHeight(parentNode) / divisor), branchStyle.getColorStart(), branchStyle.getColorFinish());
                        lineStyle = new LineStyle(width, color);
                        color = getBlendedColor((float)(tree.getHeight(node) / divisor), branchStyle.getColorStart(), branchStyle.getColorFinish());
                        lineStyle2 = new LineStyle(width, color);
                    }

                    branch = generateBranch(treeSettings.getTreeType(), tree, node, nodeName,
                            parentLatitude, latitude,
                            parentLongitude, longitude,
                            0.0, 0.0,
                            0.5, maxAltitude,
                            parentDate, date,
                            settings.getAgeCutOff(),
                            lineStyle, lineStyle2,
                            settings.getTimeDivisionCount(),
                            styles);


                } else {
                    throw new IllegalArgumentException("Unknown TreeType");
                }
                element.addContent(branch);
            }
        }

        return element;
    }

    /**
     * Create a branch that goes from [startLatitude, startLongitude, startAltitude] to
     * [finishLatitude, finishLongitude, finishAltitude], optionally arching through
     * peakAltitude. If divisionCount > 0 then the branch is subdivided into segments,
     * this is required if the branch is arcing or has a colour gradient.
     * @param treeType the tree type
     * @param tree the tree
     * @param node the node
     * @param nodeName the node's name
     * @param startLatitude the starting (parental) latitude
     * @param finishLatitude the finishing (descendent) latitude
     * @param startLongitude the starting (parental) longitude
     * @param finishLongitude the finishing (descendent) longitude
     * @param startAltitude the starting (parental) altitude
     * @param finishAltitude the finishing (descendent) altitude
     * @param peakPosition the position of the peak (0 = beginning, 0.5 = midpoint)
     * @param peakAltitude the peak altitude
     * @param startDate the starting date
     * @param finishDate the finishing date
     * @param divisionCount the number of divisions
     * @param styles an array into which style elements are put
     * @return the branch element
     */
    private Element generateBranch(TreeType treeType, RootedTree tree, Node node,
                                   String nodeName,
                                   double startLatitude, double finishLatitude,
                                   double startLongitude, double finishLongitude,
                                   double startAltitude, double finishAltitude,
                                   double peakPosition, double peakAltitude,
                                   double startDate, double finishDate,
                                   double ageCutOff,
                                   LineStyle startStyle, LineStyle finishStyle,
                                   int divisionCount,
                                   List<Element> styles) {

        Element element;

        Double rate = (Double)node.getAttribute("rate");
        Double support = (Double)node.getAttribute("posterior");
        double height = tree.getHeight(node);

        boolean hasDivisionStyles = startStyle != null && finishStyle != null;

        if (divisionCount > 0) {
            double latDiff = finishLatitude - startLatitude;
            double latDelta = latDiff / divisionCount;

            double longDiff = finishLongitude - startLongitude;
            double longDelta = longDiff / divisionCount;

            double altDiff = finishAltitude - startAltitude;
            double altDelta = altDiff / divisionCount;

            double dateDiff = finishDate - startDate;
            double dateDelta = dateDiff / divisionCount;

            double lastLatitude = startLatitude;
            double latitude = startLatitude + latDelta;

            double lastLongitude = startLongitude;
            double longitude = startLongitude + longDelta;

            // x goes 0 to 1, offset by peak Position
            double x = -peakPosition;
            double xDelta = 1.0 / divisionCount;
            // assume a parabolic curve that peaks at peakAltitude
            double a = peakAltitude - finishAltitude;
            if (peakPosition == 0.5) {
                a *= 4;
            }
            double altitude = peakAltitude - (a * (x * x));
            double lastAltitude = altitude;
            x += xDelta;

            double date = startDate;

            String styleName = null;
            if (!hasDivisionStyles && startStyle != null) {
                styleName = nodeName + "_style";

                styles.add(generateLineStyle(styleName, startStyle.getWidth(), startStyle.getColor()));
            }

            element = generateContainer("Folder", nodeName, null, (styleName != null ? "#" + styleName : null));

            for (int division = 0; division < divisionCount; division++) {
                if (ageCutOff == 0.0 || date > ageCutOff) {
                    String partName = nodeName + "_part" + (division + 1);
                    styleName = null;
                    if (hasDivisionStyles) {
                        styleName = partName + "_style";
                    }
                    Element placeMark = generateContainer("Placemark", partName, null, (styleName != null ? "#" + styleName : null));

                    annotateBranch(placeMark, height, startDate, finishDate, rate, support);

                    if (hasDivisionStyles) {
                        // Create a style for this branch segment
                        double width = startStyle.getWidth();
                        Color color = getBlendedColor(((float)division) / (divisionCount - 1), startStyle.getColor(), finishStyle.getColor());
                        styles.add(generateLineStyle(styleName, width, color));
                    }

                    if (dateDiff > 0.0) {
                        Element timeSpan = new Element("TimeSpan");

                        //convert height of the branch segment to a real date (based on th date for the most recent sample)
                        timeSpan.addContent(generateElement("begin", getKMLDate(date)));
                        placeMark.addContent(timeSpan);

                    }

                    Element lineString = new Element("LineString");
                    Element coordinates = new Element("coordinates");
                    if (altDiff > 0.0 || peakAltitude > 0.0) {
                        lineString.addContent(generateElement("altitudeMode", altitudeMode));

                        altitude = peakAltitude - (a * (x * x));

                        coordinates.addContent("" +lastLongitude + "," +lastLatitude + "," + lastAltitude + "\r");
                        coordinates.addContent("" +longitude + "," +latitude + "," + altitude + "\r");
                    } else {
                        lineString.addContent(generateElement("altitudeMode", "clampToGround"));
                        lineString.addContent(generateElement("tessellate", true));

                        coordinates.addContent("" +lastLongitude + "," +lastLatitude + "\r");
                        coordinates.addContent("" +longitude + "," +latitude + "\r");
                    }

                    lastLatitude = latitude;
                    latitude += latDelta;

                    lastLongitude = longitude;
                    longitude += longDelta;

                    lastAltitude = altitude;
                    x += xDelta;

                    date += dateDelta;

                    lineString.addContent(coordinates);

                    placeMark.addContent(lineString);

                    element.addContent(placeMark);
                }

            }

        } else {
            String styleName = null;
            if (startStyle != null) {
                styleName = nodeName + "_style";

                double width = startStyle.getWidth();
                styles.add(generateLineStyle(styleName, width, startStyle.getColor()));
            }

            element = generateContainer("Placemark", nodeName, null, (styleName != null ? "#" + styleName : null));

            annotateBranch(element, height, startDate, finishDate, rate, support);

            Element lineString = new Element("LineString");
            lineString.addContent(generateElement("altitudeMode", altitudeMode));

            Element coordinates = new Element("coordinates");
            if (treeType == TreeType.RECTANGLE_TREE) {
                coordinates.addContent(""+finishLongitude+","+finishLatitude+","+finishAltitude+"\r");
                coordinates.addContent(""+finishLongitude+","+finishLatitude+","+startAltitude+"\r");
                coordinates.addContent(""+startLongitude+","+startLatitude+","+startAltitude+"\r");
            } else { // TRIANGLE_TREE
                coordinates.addContent(""+finishLongitude+","+finishLatitude+","+finishAltitude+"\r");
                coordinates.addContent(""+startLongitude+","+startLatitude+","+startAltitude+"\r");

            }
            lineString.addContent(coordinates);

            element.addContent(lineString);

        }

        return element;
    }

    private void annotateBranch(final Element placeMark, final double height, final double startDate, final double finishDate, final Double rate, final Double support) {
        Element data = new Element("ExtendedData");
        Element schemaData = new Element("SchemaData");
        schemaData.setAttribute("schemaUrl", "#Branch_Schema");
        schemaData.addContent(new Element("SimpleData").setAttribute("name", "Height").addContent(Double.toString(height)));
        schemaData.addContent(new Element("SimpleData").setAttribute("name", "StartTime").addContent(Double.toString(startDate)));
        schemaData.addContent(new Element("SimpleData").setAttribute("name", "FinishTime").addContent(Double.toString(finishDate)));
        if (support != null) {
            schemaData.addContent(new Element("SimpleData").setAttribute("name", "Support").addContent(Double.toString(support)));
        }
        if (rate != null) {
            schemaData.addContent(new Element("SimpleData").setAttribute("name", "Rate").addContent(Double.toString(rate)));
        }
        data.addContent(schemaData);
        placeMark.addContent(data);
    }

    private void annotateTip(final Element placeMark, final String label,  final double height, final double date) {
        Element data = new Element("ExtendedData");
        Element schemaData = new Element("SchemaData");
        schemaData.setAttribute("schemaUrl", "#Node_Schema");
        schemaData.addContent(new Element("SimpleData").setAttribute("name", "Label").addContent(label));
        schemaData.addContent(new Element("SimpleData").setAttribute("name", "Height").addContent(Double.toString(height)));
        schemaData.addContent(new Element("SimpleData").setAttribute("name", "Time").addContent(Double.toString(date)));
        data.addContent(schemaData);
        placeMark.addContent(data);
    }

    private Element generateContour(SurfaceDecoration surfaces, RootedTree tree,
                                    Node node, int nodeNumber, double plotHeight, double date,
                                    String latLongName, String latitudeName, String longitudeName,
                                    boolean clampToGround, List<Element> styles) {

        String nodeName = "contour" + nodeNumber;
        String styleName = "contours";
        Element placemark = generateContainer("Placemark", nodeName, null, (styleName != null ? "#" + styleName : null));


        Element timeSpan = new Element("TimeSpan");
        timeSpan.addContent(generateElement("begin", getKMLDate(date)));
        placemark.addContent(timeSpan);

        double altitude = 0;

        if (!clampToGround) {
            double scaleFactor = plotHeight/tree.getHeight(tree.getRootNode());
            altitude = (tree.getHeight(node)*scaleFactor);
        }

        int modality = getIntegerNodeAttribute(node, latLongName+"_95%HPD_modality");

        Element container = placemark;
        if (modality > 1) {
            container = new Element("MultiGeometry");
            placemark.addContent(container);
        }

        for (int x = 0; x < modality; x++) {
            Object[] latitudeHPDs = getArrayNodeAttribute(node, latitudeName+"_95%HPD_"+(x + 1));
            Object[] longitudeHPDs = getArrayNodeAttribute(node, longitudeName+"_95%HPD_"+(x + 1));

            container.addContent(generateContourPolygon(latitudeHPDs, longitudeHPDs, altitude, clampToGround));
        }

//        if (styleName != null) {
//            styles.add(generatePolyStyle(styleName, surfaces.getStartColor()));
//        }

        return placemark;
    }

    private Element generateContourPolygon(Object[] latitudeHPDs, Object[] longitudeHPDs, double altitude, boolean clampToGround) {

        Element polygon = new Element("Polygon");
        if (clampToGround) {
            polygon.addContent(generateElement("altitudeMode", "clampToGround"));
            polygon.addContent(generateElement("tessellate", true));
        } else {
            polygon.addContent(generateElement("altitudeMode", altitudeMode));
        }
        Element outerBoundaryIs = new Element("outerBoundaryIs");
        Element linearRing = new Element("LinearRing");
        Element coordinates = new Element("coordinates");
        if (clampToGround) {
            for (int y = 0; y < longitudeHPDs.length; y++) {
                coordinates.addContent(""+longitudeHPDs[y]+","+latitudeHPDs[y]+"\r");
            }
        } else {
            for (int y = 0; y < longitudeHPDs.length; y++) {
                coordinates.addContent(""+longitudeHPDs[y]+","+latitudeHPDs[y]+","+altitude+"\r");
            }
        }
        linearRing.addContent(coordinates);
        outerBoundaryIs.addContent(linearRing);
        polygon.addContent(outerBoundaryIs);

        return polygon;
    }



    private Element generateProjection(double altitude, double latitude, double longitude) {

        Element lineString = new Element("LineString");
        lineString.addContent(generateElement("altitudeMode", altitudeMode));
        Element coordinates = new Element("coordinates");
        coordinates.addContent(""+longitude+","+latitude+","+altitude+"\r");
        coordinates.addContent(""+longitude+","+latitude+",0\r");
        lineString.addContent(coordinates);

        return lineString;
    }


    private Element generateTaxonLabel(RootedTree tree, Node node, String label, double height, double date, double latitude, double longitude) {

        Element placeMark = generateContainer("Placemark", tree.getTaxon(node).getName(), null, null);

        annotateTip(placeMark, label, height, date);

        Element point = new Element("Point");
        point.addContent(generateElement("altitudeMode", "clampToGround"));
        Element coordinates = new Element("coordinates");
        coordinates.addContent(""+longitude+","+latitude+"\r");
        point.addContent(coordinates);
        placeMark.addContent(point);

        return placeMark;
    }


    private Element generateLocation(Location location) {

        Element placeMark = generateContainer("Placemark", location.getName(), null, null);

        Element point = new Element("Point");
        point.addContent(generateElement("altitudeMode", "clampToGround"));
        Element coordinates = new Element("coordinates");
        coordinates.addContent("" + location.getLongitude() + "," + location.getLatitude() + ",0\r");
        point.addContent(coordinates);
        placeMark.addContent(point);

        return placeMark;
    }

    private Element generateLineStyle(String styleName, double width, Color color) {
        Element style = new Element("Style");
        style.setAttribute("id", styleName);

        Element lineStyle = new Element("LineStyle");
        lineStyle.addContent(generateElement("width", "" + width));
        lineStyle.addContent(generateElement("color", getKMLColor(color)));
        style.addContent(lineStyle);

        return style;
    }

    private Element generatePolyStyle(String styleName, Color color) {
        Element style = new Element("Style");
        style.setAttribute("id", styleName);

        Element polyStyle = new Element("PolyStyle");
        polyStyle.addContent(generateElement("color", getKMLColor(color)));
        polyStyle.addContent(generateElement("outline", false));
        style.addContent(polyStyle);

        return style;
    }

    private Element generatePolyStyle(String styleName, Color color, double outlineWidth, Color outlineColor) {
        Element style = new Element("Style");
        style.setAttribute("id", styleName);

        Element lineStyle = new Element("LineStyle");
        lineStyle.addContent(generateElement("width", "" + outlineWidth));
        lineStyle.addContent(generateElement("color", getKMLColor(outlineColor)));
        style.addContent(lineStyle);

        Element polyStyle = new Element("PolyStyle");
        polyStyle.addContent(generateElement("color", getKMLColor(color)));
        polyStyle.addContent(generateElement("outline", true));
        style.addContent(polyStyle);

        return style;
    }


    private Element generateContainer(String elementTag, String name, String description, String styleURL) {
        Element element = new Element(elementTag);
        if (name != null) {
            element.addContent(generateElement("name", name));
        }
        if (description != null) {
            element.addContent(generateElement("description", description));
        }
        if (styleURL != null) {
            element.addContent(generateElement("styleUrl", styleURL));
        }
        return element;
    }

    private Element generateContainer(String elementTag, String name, String description, boolean visibility, String styleURL) {
        Element element = new Element(elementTag);
        if (name != null) {
            element.addContent(generateElement("name", name));
        }
        if (description != null) {
            element.addContent(generateElement("description", description));
        }
        if (styleURL != null) {
            element.addContent(generateElement("styleUrl", styleURL));
        }
        element.addContent(generateElement("visibility", visibility));
        return element;
    }

    private Element generateElement(String elementName, boolean content) {
        Element e = new Element(elementName);
        e.addContent(content ? "1" : "0");
        return e;
    }

    private Element generateElement(String elementName, String content) {
        Element e = new Element(elementName);
        e.addContent(content);
        return e;
    }

    private static String getKMLDate(double fractionalDate) {

        int year = (int) fractionalDate;
        String yearString;

        if (year < 10) {
            yearString = "000"+year;
        } else if (year < 100) {
            yearString = "00"+year;
        } else if (year < 1000) {
            yearString = "0"+year;
        } else {
            yearString = ""+year;
        }

        double fractionalMonth = fractionalDate - year;

        int month = (int) (12.0 * fractionalMonth);
        String monthString;

        if (month < 10) {
            monthString = "0"+month;
        } else {
            monthString = ""+month;
        }

        int day = (int) Math.round(30*(12*fractionalMonth - month));
        String dayString;

        if (day < 10) {
            dayString = "0"+day;
        } else {
            dayString = ""+day;
        }


        return yearString + "-" + monthString + "-" + dayString;
    }


    private int getIntegerNodeAttribute(Node node, String attributeName) {
        if (node.getAttribute(attributeName) == null) {
            throw new RuntimeException("Attribute, " + attributeName + ", missing from node");
        }
        return (Integer)node.getAttribute(attributeName);
    }

    private int getIntegerNodeAttribute(Node node, String attributeName, int defaultValue) {
        if (node.getAttribute(attributeName) == null) {
            return defaultValue;
        }
        return (Integer)node.getAttribute(attributeName);
    }

    private double getDoubleNodeAttribute(Node node, String attributeName) {
        if (node.getAttribute(attributeName) == null) {
            throw new RuntimeException("Attribute, " + attributeName + ", missing from node");
        }
        return (Double)node.getAttribute(attributeName);
    }

    private double getDoubleNodeAttribute(Node node, String attributeName, double defaultValue) {
        if (node.getAttribute(attributeName) == null) {
            return defaultValue;
        }
        return (Double)node.getAttribute(attributeName);
    }

    private Object[] getArrayNodeAttribute(Node node, String attributeName) {
        if (node.getAttribute(attributeName) == null) {
            return null;
//            throw new RuntimeException("Attribute, " + attributeName + ", missing from node");
        }
        return (Object[])node.getAttribute(attributeName);
    }

    /**
     * converts a Java color into a 4 channel hex color string.
     * @param color
     * @return the color string
     */
    public static String getKMLColor(Color color) {
        String a = Integer.toHexString(color.getAlpha());
        String b = Integer.toHexString(color.getBlue());
        String g = Integer.toHexString(color.getGreen());
        String r = Integer.toHexString(color.getRed());
        return  (a.length() < 2 ? "0" : "") + a +
                (b.length() < 2 ? "0" : "") + b +
                (g.length() < 2 ? "0" : "") + g +
                (r.length() < 2 ? "0" : "") + r;
    }

    /**
     * converts a Java color into a 4 channel hex color string.
     * @param color
     * @return the color string
     */
    public static String getKMLColor(Color color, double opacity) {
        int alpha = (int)(256 * (1.0 - opacity));
        String a = Integer.toHexString(alpha);
        String b = Integer.toHexString(color.getBlue());
        String g = Integer.toHexString(color.getGreen());
        String r = Integer.toHexString(color.getRed());
        return  (a.length() < 2 ? "0" : "") + a +
                (b.length() < 2 ? "0" : "") + b +
                (g.length() < 2 ? "0" : "") + g +
                (r.length() < 2 ? "0" : "") + r;
    }

    public static Color getBlendedColor(float proportion, Color startColor, Color endColor) {
        proportion = Math.max(proportion, 0.0F);
        proportion = Math.min(proportion, 1.0F);
        float[] start = startColor.getRGBColorComponents(null);
        float[] end = endColor.getRGBColorComponents(null);

        float[] color = new float[start.length];
        for (int i = 0; i < start.length; i++) {
            color[i] = start[i] + ((end[i] - start[i]) * proportion);
        }

        return new Color(color[0], color[1], color[2]);
    }

    public static void main(String[] args) {


        String inputTreeFile = args[0];
        RootedTree tree = null;

        try {
            TreeImporter importer = new NexusImporter(new FileReader(inputTreeFile));
            tree = (RootedTree)importer.importNextTree();
        } catch (ImportException e) {
            e.printStackTrace();
            return;
        } catch (IOException e) {
            e.printStackTrace();
            return;
        }

        TreeKMLGenerator generator = new TreeKMLGenerator(tree);
        Settings settings = new Settings(AnalysisType.CONTINUOUS);
//        settings.getAltitudeTreeSettings().setTreeType(TreeType.ARC_TREE);
//        settings.getAltitudeTreeSettings().getBranchStyle().setColorProperty("height");
        settings.getGroundTreeSettings().setTreeType(TreeType.SURFACE_TREE);
        settings.getGroundTreeSettings().getBranchStyle().setColorProperty("height");
        settings.setPlotAltitude(0);
        settings.setMostRecentDate(2003);
        //settings.setAgeCutOff(1995);
        settings.setTimeDivisionCount(0);

        settings.setTraitName("antigenic");
        settings.setLatitudeName("antigenic1");
        settings.setLongitudeName("antigenic2");

        try {

            BufferedWriter out = new BufferedWriter(new FileWriter(args[0]+".kml"));
            Document doc = new Document(generator.generate(args[0], settings));

            try {
                XMLOutputter outputter = new XMLOutputter(Format.getPrettyFormat());
                outputter.output(doc, out);
            }
            catch (IOException e) {
                System.err.println(e);
            }
            out.close();
        } catch (IOException e) {
            e.printStackTrace();
            return;
        }

    }

}
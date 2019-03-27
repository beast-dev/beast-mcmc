/*
 * TreeGeoJSONGenerator.java
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
import jebl.evolution.graphs.Node;
import jebl.evolution.io.ImportException;
import jebl.evolution.io.NexusImporter;
import jebl.evolution.io.TreeImporter;
import jebl.evolution.trees.RootedTree;

import java.awt.*;
import java.io.*;
import java.util.*;
import java.util.List;

/**
 * @author Andrew Rambaut
 */
public class TreeGeoJSONGenerator {
    /*
        { "type": "FeatureCollection",
                "features": [
            { "type": "Feature",
                    "geometry": {"type": "Point", "coordinates": [102.0, 0.5]},
                "properties": {"prop0": "value0"}
            },
            { "type": "Feature",
                    "geometry": {
                "type": "LineString",
                        "coordinates": [
                [102.0, 0.0], [103.0, 1.0], [104.0, 0.0], [105.0, 1.0]
                ]
            },
                "properties": {
                "prop0": "value0",
                        "prop1": 0.0
            }
            },
            { "type": "Feature",
                    "geometry": {
                "type": "Polygon",
                        "coordinates": [
                [ [100.0, 0.0], [101.0, 0.0], [101.0, 1.0],
                [100.0, 1.0], [100.0, 0.0] ]
                ]
            },
                "properties": {
                "prop0": "value0",
                        "prop1": {"this": "that"}
            }
            }
            ]
        }
    */


    private final RootedTree tree;
    private final int locationHPDPercentile;
    private final Map<String, Location> locationMap;
    private final List<GeoPoint> points = new ArrayList<GeoPoint>();
    private final List<GeoLine> lines = new ArrayList<GeoLine>();
    private final List<GeoPoly> polys = new ArrayList<GeoPoly>();

    public TreeGeoJSONGenerator(RootedTree tree, List<Location> locationList, int locationHPDPercentile) {
        this.tree = tree;

        this.locationHPDPercentile = locationHPDPercentile;

        if (locationList != null) {
            this.locationMap = new HashMap<String, Location>();
            double minLat = Double.MAX_VALUE;
            double maxLat = -Double.MAX_VALUE;
            double minLong = Double.MAX_VALUE;
            double maxLong = -Double.MAX_VALUE;

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

        } else {
            this.locationMap = null;
        }

        traverseTree(tree, tree.getRootNode());
    }

    private GeoPoint traverseTree(RootedTree tree, Node node) {
        GeoPoint point = new GeoPoint();
        point.time = tree.getHeight(node);
        if (locationMap != null) {
            point.location = (String) node.getAttribute("location");
//            point.host = (String) node.getAttribute("host");
            point.longitude = locationMap.get(point.location).getLongitude();
            point.latitude = locationMap.get(point.location).getLatitude();
        } else {
            // is continuous
            point.latitude = Double.parseDouble(node.getAttribute("Location1").toString());
            point.longitude = Double.parseDouble(node.getAttribute("Location2").toString());

            GeoPoly poly = new GeoPoly();
            poly.label = "hpd_" + locationHPDPercentile;

            if (node.getAttribute("Location_" + locationHPDPercentile + "%HPD_modality") != null) {
                int modality = Integer.parseInt(node.getAttribute("Location_" + locationHPDPercentile + "%HPD_modality").toString());

                poly.time = point.time;
                poly.latitude = new Double[modality][];
                poly.longitude = new Double[modality][];

                for (int i = 0; i < modality; i++) {
                    Object[] lat = (Object[]) node.getAttribute("Location1_" + locationHPDPercentile + "%HPD_" + (i + 1));
                    poly.latitude[i] = Arrays.copyOf(lat, lat.length, Double[].class);
                    Object[] lon = (Object[]) node.getAttribute("Location2_" + locationHPDPercentile + "%HPD_" + (i + 1));
                    poly.longitude[i] = Arrays.copyOf(lon, lon.length, Double[].class);
                }

                polys.add(poly);
            }
        }

        point.isExternal = tree.isExternal(node);

        if (point.isExternal) {
            point.label = tree.getTaxon(node).getName();
        } else {
            point.label = "";
            point.probability = (Double)node.getAttribute("posterior");

            for (Node child : tree.getChildren(node)) {

                GeoPoint point1 = traverseTree(tree, child);

                GeoLine line = new GeoLine();
                line.label = "";

                line.time0 = point.time;
                line.time1 = point1.time;

                line.location0 = point.location;
                line.location1 = point1.location;

//                line.host0 = point.host;
//                line.host1 = point1.host;

                line.longitude0 = point.longitude;
                line.longitude1 = point1.longitude;

                line.latitude0 = point.latitude;
                line.latitude1 = point1.latitude;

                lines.add(line);
            }
        }

        points.add(point);

        return point;
    }

    public void generate(String documentName, String fileNameStem) {
        if (points.size() > 0) {
            try {
                PrintWriter writer = new PrintWriter(fileNameStem + "_tips.geojson");

                writer.println("{ \"type\": \"FeatureCollection\",");
                writer.println("\t\"features\": [");
                for (GeoPoint point : points) {
                    if (point.isExternal) {
                        generatePoint(writer, point);
                    }
                }
                writer.println("\t]");
                writer.println("}");

                writer.close();
            } catch (FileNotFoundException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }

            try {
                PrintWriter writer = new PrintWriter(fileNameStem + "_nodes.geojson");

                writer.println("{ \"type\": \"FeatureCollection\",");
                writer.println("\t\"features\": [");
                for (GeoPoint point : points) {
                    if (!point.isExternal) {
                        generatePoint(writer, point);
                    }
                }
                writer.println("\t]");
                writer.println("}");

                writer.close();
            } catch (FileNotFoundException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
        }

        if (lines.size() > 0) {
            try {
                PrintWriter writer = new PrintWriter(fileNameStem + "_tree.geojson");
                writer.println("{ \"type\": \"FeatureCollection\",");
                writer.println("\t\"features\": [");
                for (GeoLine line : lines) {
                    generateLine(writer, line);
                }
                writer.println("\t]");
                writer.println("}");
                writer.close();
            } catch (FileNotFoundException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
        }

        if (polys.size() > 0) {
            try {
                PrintWriter writer = new PrintWriter(fileNameStem + "_hpds.geojson");

                writer.println("{ \"type\": \"FeatureCollection\",");
                writer.println("\t\"features\": [");
                for (GeoPoly poly : polys) {
                    generatePoly(writer, poly);
                }
                writer.println("\t]");
                writer.println("}");
                writer.close();
            } catch (FileNotFoundException e) {
                e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
            }
        }

    }

    private void generatePoint(PrintWriter writer, GeoPoint point) {
        writer.println("\t\t{ \"type\": \"Feature\",");
        writer.println("\t\t\t\"geometry\": {\"type\": \"Point\", \"coordinates\": [" + point.longitude + ", " + point.latitude + "]},");
        if (locationMap != null) {
            if (point.isExternal) {
                writer.println("\t\t\t\"properties\": {\"label\": \"" + point.label + "\", \"time\": \"" + point.time + "\", \"probability\": \"" + point.probability +
                        "\", \"location\": \"" + point.location + "\"}");
            } else {
                writer.println("\t\t\t\"properties\": {\"label\": \"" + point.label + "\", \"time\": \"" + point.time + "\", \"probability\": \"" + point.probability +
                        "\", \"location\": \"" + point.location + "\"}");
            }
        } else {
            if (point.isExternal) {
                writer.println("\t\t\t\"properties\": {\"label\": \"" + point.label + "\", \"time\": \"" + point.time + "\"}");
            } else {
                writer.println("\t\t\t\"properties\": {\"label\": \"" + point.label + "\", \"time\": \"" + point.time + "\", \"probability\": \"" + point.probability + "\"}");
            }
        }
        writer.println("\t\t},");
    }

    private void generateLine(PrintWriter writer, GeoLine line) {
        writer.println("\t\t{ \"type\": \"Feature\",");
        writer.println("\t\t\t\"geometry\": {\"type\": \"LineString\", \"coordinates\": [[" + line.longitude0 + ", " + line.latitude0 + "], [" + line.longitude1 + ", " + line.latitude1 + "]]},");
        if (locationMap != null) {
            writer.println("\t\t\t\"properties\": {\"label\": \"" + line.label + "\", \"time0\": \"" + line.time0 + "\", \"time0\": \"" + line.time1 +
                    "\", \"location0\": \"" + line.location0 + "\", \"location1\": \"" + line.location1 +
                    "\", \"host\": \"" + line.host0 + "\", \"host1\": \"" + line.host1 + "\"}");
        } else {
            writer.println("\t\t\t\"properties\": {\"label\": \"" + line.label + "\", \"time0\": \"" + line.time0 + "\", \"time0\": \"" + line.time1 + "\"}");
        }
        writer.println("},");
    }

    private void generatePoly(PrintWriter writer, GeoPoly poly) {
        writer.println("\t\t{ \"type\": \"Feature\",");
        writer.print("\t\t\t\"geometry\": {\"type\": \"Polygon\", \"coordinates\": [");
        for (int i = 0; i < poly.latitude.length; i++) {
            if (i > 0) {
                writer.print(", ");
            }
            writer.print("[");
            for (int j = 0; j < poly.latitude[i].length; j++) {
                if (j > 0) {
                    writer.print(", ");
                }
                writer.print("[" + poly.longitude[i][j] + ", " + poly.latitude[i][j] + "]");
            }
            writer.print("]");
        }
        writer.println("]},");
        writer.println("\t\t\t\"properties\": {\"label\": \"" + poly.label + "\", \"time\": \"" + poly.time + "\"}");
        writer.println("\t\t},");
    }

    class GeoPoint {
        String label;
        double time;
        double probability;
        String location;
        String host;
        double longitude;
        double latitude;
        boolean isExternal;
    }

    class GeoLine {
        String label;
        double time0;
        double time1;
        String location0;
        String location1;
        String host0;
        String host1;
        double longitude0;
        double latitude0;
        double longitude1;
        double latitude1;
    }

    class GeoPoly {
        String label;
        double time;
        Double[][] longitude;
        Double[][] latitude;
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

        TreeGeoJSONGenerator generator = new TreeGeoJSONGenerator(tree, null, 80);
        generator.generate("", args[1]);
    }

}
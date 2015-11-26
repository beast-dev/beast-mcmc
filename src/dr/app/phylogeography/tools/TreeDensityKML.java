/*
 * TreeDensityKML.java
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

import dr.app.beast.BeastVersion;
import dr.app.phylogeography.tools.kml.Location;
import dr.app.util.Arguments;
import dr.util.Version;
import jebl.evolution.graphs.Node;
import jebl.evolution.io.ImportException;
import jebl.evolution.io.NexusImporter;
import jebl.evolution.io.TreeImporter;
import jebl.evolution.trees.RootedTree;
import jebl.math.Random;
import org.jdom.Document;
import org.jdom.Element;
import org.jdom.output.Format;
import org.jdom.output.XMLOutputter;

import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/*
 * @author Andrew Rambaut
 */

public class TreeDensityKML {
    private static final String STATE_ATTRIBUTE_NAME = "states";

    private final static Version version = new BeastVersion();

    private final Map<String, Location> locationMap;

    public TreeDensityKML(int burnin,
                          int skipEvery,
                          String inputFileName,
                          String outputFileName) throws IOException {

        locationMap = readCoordinates("locations.txt");

        System.out.println("Reading trees...");

        Element root = new Element("kml");

        Element doc = new Element("Document");
        doc.addContent(generateElement("name", inputFileName));

        root.addContent(doc);

//        List<Element> schema = new ArrayList<Element>();
//        Element branchSchema = new Element("Schema");
//        branchSchema.setAttribute("id", "Branch_Schema");
//        branchSchema.addContent(new Element("SimpleField")
//                .setAttribute("name", "Name")
//                .setAttribute("type", "string")
//                .addContent(new Element("displayName").addContent("Name")));
//        branchSchema.addContent(new Element("SimpleField")
//                .setAttribute("name", "StartTime")
//                .setAttribute("type", "double")
//                .addContent(new Element("displayName").addContent("StartTime")));
//        branchSchema.addContent(new Element("SimpleField")
//                .setAttribute("name", "FinishTime")
//                .setAttribute("type", "double")
//                .addContent(new Element("displayName").addContent("FinishTime")));
//        branchSchema.addContent(new Element("SimpleField")
//                .setAttribute("name", "Rate")
//                .setAttribute("type", "double")
//                .addContent(new Element("displayName").addContent("Rate")));
//        schema.add(branchSchema);

        List<Element> trees = new ArrayList<Element>();

        boolean firstTree = true;
        FileReader fileReader = new FileReader(inputFileName);
        TreeImporter importer = new NexusImporter(fileReader);
        try {
            while (importer.hasTree()) {
                RootedTree tree = (RootedTree)importer.importNextTree();

                if (firstTree) {
                    firstTree = false;
                }

                if (totalTrees >= burnin) {
                    if (totalTrees % skipEvery == 0) {
                        trees.add(generateKMLTree("Tree_" + totalTrees, tree));

                        if (totalTrees % 50 == 0) {
                            System.out.print(".");
                        }
                        totalTreesUsed += 1;
                    }
                }
                totalTrees += 1;

            }

            System.out.println();

        } catch (ImportException e) {
            System.err.println("Error Parsing Input Tree: " + e.getMessage());
            return;
        }
        fileReader.close();

//        doc.addContent(schema);
//        doc.addContent(styles);
        doc.addContent(trees);

        try {

            BufferedWriter out = new BufferedWriter(new FileWriter(outputFileName));
            Document document = new Document(root);

            try {
                XMLOutputter outputter = new XMLOutputter(Format.getPrettyFormat());
                outputter.output(document, out);
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

    private Element generateKMLTree(String name, RootedTree tree) {


        Element element = generateContainer("Folder", name, null, null);

        double longNoise = Random.nextGaussian() * 0.5;
        double latNoise = Random.nextGaussian() * 0.5;

        int nodeNumber = 0;
        for (Node node : tree.getNodes()) {
            nodeNumber++;

            if (!tree.isRoot(node)) {
                String state = (String)node.getAttribute(STATE_ATTRIBUTE_NAME);
                Location location = locationMap.get(state);
                if (location == null) {
                    throw new RuntimeException("No location called " + state + " in location list");
                }

                // Create each branch of the tree..

                String nodeName = name + "_node" + nodeNumber;


                Node parentNode = tree.getParent(node);
                String parentState = (String)parentNode.getAttribute(STATE_ATTRIBUTE_NAME);

                Location parentLocation = locationMap.get(parentState);
                if (parentLocation == null) {
                    throw new RuntimeException("No location called " + parentState + " in location list");
                }

                Element branch = generateContainer("Placemark", nodeName, null, null);

//                annotateBranch(element, height, startDate, finishDate, rate, support);

                Element lineString = new Element("LineString");
                lineString.addContent(generateElement("altitudeMode", "clampToGround"));

                Element coordinates = new Element("coordinates");
                    coordinates.addContent(""+ (parentLocation.getLongitude() + longNoise)+","+ (parentLocation.getLatitude()  + latNoise)+"\r");
                    coordinates.addContent(""+ (location.getLongitude() + longNoise) +","+ (location.getLatitude() + latNoise) +"\r");
                lineString.addContent(coordinates);

                branch.addContent(lineString);

                element.addContent(branch);
            }
        }

        return element;
    }


//    private void annotateBranch(final Element placeMark, final double height, final double startDate, final double finishDate, final Double rate, final Double support) {
//        Element data = new Element("ExtendedData");
//        Element schemaData = new Element("SchemaData");
//        schemaData.setAttribute("schemaUrl", "#Branch_Schema");
//        schemaData.addContent(new Element("SimpleData").setAttribute("name", "Height").addContent(Double.toString(height)));
//        schemaData.addContent(new Element("SimpleData").setAttribute("name", "StartTime").addContent(Double.toString(startDate)));
//        schemaData.addContent(new Element("SimpleData").setAttribute("name", "FinishTime").addContent(Double.toString(finishDate)));
//        if (support != null) {
//            schemaData.addContent(new Element("SimpleData").setAttribute("name", "Support").addContent(Double.toString(support)));
//        }
//        if (rate != null) {
//            schemaData.addContent(new Element("SimpleData").setAttribute("name", "Rate").addContent(Double.toString(rate)));
//        }
//        data.addContent(schemaData);
//        placeMark.addContent(data);
//    }

    private int getIntegerNodeAttribute(Node node, String attributeName) {
        if (node.getAttribute(attributeName) == null) {
            throw new RuntimeException("Attribute, " + attributeName + ", missing from node");
        }
        return (Integer)node.getAttribute(attributeName);
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


    private Element generateElement(String elementName, String content) {
        Element e = new Element(elementName);
        e.addContent(content);
        return e;
    }

    private static Map<String, Location> readCoordinates(String fileName) {
        Map<String, Location> locationMap = new HashMap<String, Location>();

        try {
            BufferedReader reader = new BufferedReader(new FileReader(fileName));
            String line = reader.readLine();
            while (line != null && line.trim().length() > 0) {
                String[] parts = line.split("\t");
                Location location;
                if (parts.length == 4) {
                     location = new Location(parts[0], parts[1], Double.parseDouble(parts[2]), Double.parseDouble(parts[3]));
                } else if (parts.length == 3) {
                    location = new Location(parts[0], parts[0], Double.parseDouble(parts[1]), Double.parseDouble(parts[2]));
                } else {
                    throw new RuntimeException("Wrong number of columns in coordinates file");
                }

                locationMap.put(location.getState(), location);
                line = reader.readLine();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return locationMap;
    }

    int totalTrees = 0;
    int totalTreesUsed = 0;

    public static void printTitle() {
        System.out.println();
        centreLine("TreeDensityKML " + version.getVersionString() + ", " + version.getDateString(), 60);
        centreLine("by", 60);
        centreLine("Andrew Rambaut", 60);
        System.out.println();
        System.out.println();
    }

    public static void centreLine(String line, int pageWidth) {
        int n = pageWidth - line.length();
        int n1 = n / 2;
        for (int i = 0; i < n1; i++) {
            System.out.print(" ");
        }
        System.out.println(line);
    }


    public static void printUsage(Arguments arguments) {

        arguments.printUsage("TreeDensityKML", "<input-file-name> [<output-file-name>]");
        System.out.println();
        System.out.println("  Example: TreeDensityKML -burnin 100 test.trees test.kml");
        System.out.println();
    }

    //Main method
    public static void main(String[] args) throws IOException {

        String inputFileName = null;
        String outputFileName = null;

        printTitle();

        Arguments arguments = new Arguments(
                new Arguments.Option[]{
                        new Arguments.IntegerOption("burnin", "the number of states to be considered as 'burn-in' [default = 0]"),
                        new Arguments.IntegerOption("skip", "skip over this many trees per sample [default = 1]"),
                        new Arguments.Option("help", "option to print this message")
                });

        try {
            arguments.parseArguments(args);
        } catch (Arguments.ArgumentException ae) {
            System.out.println(ae);
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

        int skipEvery = 1;
        if (arguments.hasOption("skip")) {
            skipEvery = arguments.getIntegerOption("skip");
        }

        String[] args2 = arguments.getLeftoverArguments();

        if (args2.length > 2) {
            System.err.println("Unknown option: " + args2[2]);
            System.err.println();
            printUsage(arguments);
            System.exit(1);
        }

        if (args2.length == 2) {
            inputFileName = args2[0];
            outputFileName = args2[1];
        } else {
            System.err.println("Missing input or output file name");
            printUsage(arguments);
            System.exit(1);
        }

        new TreeDensityKML(burnin,
                skipEvery,
                inputFileName,
                outputFileName
        );

        System.exit(0);
    }

}
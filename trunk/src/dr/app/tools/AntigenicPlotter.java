/*
 * TreeAnnotator.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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
import dr.geo.KMLCoordinates;
import dr.geo.KernelDensityEstimator2D;
import dr.geo.contouring.*;
import dr.inference.trace.*;
import dr.util.Version;
import org.jdom.Element;
import org.jdom.output.Format;
import org.jdom.output.XMLOutputter;

import java.io.*;
import java.util.Calendar;

/*
 * @author Marc Suchard and Andrew Rambaut
 */

public class AntigenicPlotter {

    private final static Version version = new BeastVersion();

    public final static ContourMode CONTOUR_MODE = ContourMode.SNYDER;
    public final static double HPD_VALUE = 0.95;
    public static final int GRIDSIZE = 200;

    private Element rootElement;
    private Element documentElement;
    private Element contourFolderElement;
    private Element traceFolderElement;

    public AntigenicPlotter(int burnin,
                            String inputFileName,
                            String outputFileName
    ) throws IOException {

        System.out.println("Reading log file...");

        FileReader fileReader = new FileReader(inputFileName);
        try {
            File file = new File(inputFileName);

            LogFileTraces traces = new LogFileTraces(inputFileName, file);
            traces.loadTraces();


            if (burnin == -1) {
                burnin = traces.getMaxState() / 10;
            }

            traces.setBurnIn(burnin);

            System.out.println();
            System.out.println("burnIn   <= " + burnin);
            System.out.println("maxState  = " + traces.getMaxState());
            System.out.println();

            writeKML(outputFileName, traces);

        } catch (Exception e) {
            System.err.println("Error Parsing Input File: " + e.getMessage());
            return;
        }
        fileReader.close();

    }

    private void writeKML(String fileName, LogFileTraces traces) {
        Element hpdSchema = new Element("Schema");
        hpdSchema.setAttribute("schemaUrl", "HPD_Schema");
        hpdSchema.addContent(new Element("SimpleField")
                .setAttribute("name", "Name")
                .setAttribute("type", "string")
                .addContent(new Element("displayName").addContent("Name")));
        hpdSchema.addContent(new Element("SimpleField")
                .setAttribute("name", "Point")
                .setAttribute("type", "integer")
                .addContent(new Element("displayName").addContent("Point")));
        hpdSchema.addContent(new Element("SimpleField")
                .setAttribute("name", "HPD")
                .setAttribute("type", "double")
                .addContent(new Element("displayName").addContent("HPD")));

        Element nodeSchema = new Element("Schema");
        nodeSchema.setAttribute("schemaUrl", "Trace_Schema");
        nodeSchema.addContent(new Element("SimpleField")
                .setAttribute("name", "State")
                .setAttribute("type", "integer")
                .addContent(new Element("displayName").addContent("State")));

        contourFolderElement = new Element("Folder");
        Element contourFolderNameElement = new Element("name");
        contourFolderNameElement.addContent("HPDs");
        contourFolderElement.addContent(contourFolderNameElement);

        traceFolderElement = new Element("Folder");
        Element traceFolderNameElement = new Element("name");
        traceFolderNameElement.addContent("traces");
        traceFolderElement.addContent(traceFolderNameElement);

        Element documentNameElement = new Element("name");
        String documentName = fileName;
        if (documentName.endsWith(".kml"))
            documentName = documentName.replace(".kml", "");
        documentNameElement.addContent(documentName);

        documentElement = new Element("Document");
        documentElement.addContent(documentNameElement);
//        documentElement.addContent(hpdSchema);
//        documentElement.addContent(nodeSchema);
        documentElement.addContent(contourFolderElement);
        documentElement.addContent(traceFolderElement);

        rootElement = new Element("kml");
        rootElement.addContent(documentElement);

        int pointCount = traces.getTraceCount() / 2;

        int stateCount = traces.getStateCount();
        double[][] xy = new double[2][stateCount];

        for (int i = 0; i < pointCount; i++) {

            for (int j = 0; j < stateCount; j++) {
                xy[0][j] = traces.getStateValue(i * 2, j) / 1000.0;
                xy[1][j] = traces.getStateValue((i * 2) + 1, j) / 1000.0;
            }

            Element traceElement = generateTraceElement(i + 1, xy);
            traceFolderElement.addContent(traceElement);

            ContourMaker contourMaker;
            if (CONTOUR_MODE == ContourMode.JAVA)
                contourMaker = new KernelDensityEstimator2D(xy[0], xy[1], GRIDSIZE);
            else if (CONTOUR_MODE == ContourMode.R)
                contourMaker = new ContourWithR(xy[0], xy[1], GRIDSIZE);
            else if (CONTOUR_MODE == ContourMode.SNYDER)
                contourMaker = new ContourWithSynder(xy[0], xy[1], GRIDSIZE);
            else
                throw new RuntimeException("Unimplemented ContourModel!");

            ContourPath[] paths = contourMaker.getContourPaths(HPD_VALUE);
            int pathCounter = 1;
            for (ContourPath path : paths) {

                KMLCoordinates coords = new KMLCoordinates(path.getAllX(), path.getAllY());

                //because KML polygons require long,lat,alt we need to switch lat and long first
                coords.switchXY();
                Element placemarkElement = generatePlacemarkElementWithPolygon(HPD_VALUE, Double.NaN, coords, -1, pathCounter);
                //testing how many points are within the polygon
                contourFolderElement.addContent(placemarkElement);

                pathCounter ++;
            }


        }

        PrintStream resultsStream;

        try {
            resultsStream = new PrintStream(new File(fileName));
            XMLOutputter xmlOutputter = new XMLOutputter(Format.getPrettyFormat().setTextMode(Format.TextMode.PRESERVE));
            xmlOutputter.output(rootElement, resultsStream);

        } catch (IOException e) {
            System.err.println("Error opening file: " + fileName);
            System.exit(-1);
        }
        
    }

    private Element generateTraceElement(int pointNumber, double[][] points) {
        Element traceElement = new Element("Folder");
        Element nameKDEElement = new Element("name");
        String name = "trace_" + pointNumber;

        nameKDEElement.addContent(name);
        traceElement.addContent(nameKDEElement);

        for (int a = 0; a < points[0].length; a++)  {
            Element placemarkElement = new Element("Placemark");

//            placemarkElement.addContent(generateTraceData(a));


            Element pointElement = new Element("Point");
            Element coordinates = new Element("coordinates");
            coordinates.addContent(points[1][a]+","+points[0][a]+",0");
            pointElement.addContent(coordinates);
            placemarkElement.addContent(pointElement);

            traceElement.addContent(placemarkElement);

        }

        return traceElement;
    }

    private Element generateTraceData(int state) {
        Element data = new Element("ExtendedData");
        Element schemaData = new Element("SchemaData");
        schemaData.setAttribute("schemaUrl", "Trace_Schema");
        schemaData.addContent(new Element("SimpleData").setAttribute("state", "State").addContent(Integer.toString(state)));
        data.addContent(schemaData);
        return data;
    }

    private Element generatePlacemarkElementWithPolygon(double hpdValue, double sliceValue, KMLCoordinates coords, int pointNumber, int pathCounter) {
        Element placemarkElement = new Element("Placemark");

        String name;
        Element placemarkNameElement = new Element("name");
        name = "hpdRegion_" + sliceValue + "_path_" + pathCounter;

        placemarkNameElement.addContent(name);
        placemarkElement.addContent(placemarkNameElement);

//        placemarkElement.addContent(generateContourData(name, pointNumber, hpdValue));

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

    private Element generateContourData(String name, int point, double hpd) {
        Element data = new Element("ExtendedData");
        Element schemaData = new Element("SchemaData");
        schemaData.setAttribute("schemaUrl", "HPD_Schema");
        schemaData.addContent(new Element("SimpleData").setAttribute("name", "Name").addContent(name));
        schemaData.addContent(new Element("SimpleData").setAttribute("point", "Point").addContent(Integer.toString(point)));
        if (hpd > 0) {
            schemaData.addContent(new Element("SimpleData").setAttribute("name", "HPD").addContent(Double.toString(hpd)));
        }
        data.addContent(schemaData);
        return data;
    }

    public static void printTitle() {
        System.out.println();
        centreLine("AntigenicPlotter " + version.getVersionString() + ", " + version.getDateString(), 60);
        centreLine("BEAST time vs. parameter density analysis", 60);
        centreLine("by", 60);
        centreLine("Andrew Rambaut and Marc A. Suchard", 60);
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

        arguments.printUsage("antigenicplotter", "<input-file-name> [<output-file-name>]");
        System.out.println();
        System.out.println("  Example: antigenicplotter -burnin 100 locations.log locations.kml");
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

        new AntigenicPlotter(burnin,
                inputFileName,
                outputFileName
        );

        System.exit(0);
    }

}
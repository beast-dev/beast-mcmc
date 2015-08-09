/*
 * TopographicalMap.java
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

package dr.evomodel.continuous;

import dr.math.SparseMatrixExponential;
import dr.math.matrixAlgebra.Vector;
import org.jdom.Element;
import org.jdom.output.XMLOutputter;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

/**
 * @author Marc A. Suchard
 */
public class TopographicalMap {


    public static final String FORMAT = "%3.2f";


    public TopographicalMap(double xStart, double xEnd, double yStart, double yEnd, int xDim, int yDim) {
        this.xStart = xStart;
        this.xEnd = xEnd;
        this.yStart = yStart;
        this.yEnd = yEnd;

        this.xDim = xDim;
        this.yDim = yDim;
        map = new double[xDim][yDim];

        xDelta = (xEnd - xStart) / xDim;
        yDelta = (yEnd - yStart) / yDim;
    }

    private double xStart;
    private double xEnd;
    private double yStart;
    private double yEnd;
    private double xDelta;
    private double yDelta;

//        public TopographicalMap(int xDim, int yDim) {
//            TopographicalMap(0,1,0,1,xDim,yDim);
//        }

    public TopographicalMap(double[][] map) {
        xDim = map.length;
        yDim = map[0].length;
        this.map = map;
        order = setUpIndices();
        setUpRandomWalkMatrix();
    }

    public int getXDim() {
        return xDim;
    }

    public int getYDim() {
        return yDim;
    }

    public int getOrder() {
        return order;
    }

    public int getNonZeroSize() {
        return nonZeroElements;
    }


    private int setUpIndices() {
        int index = 0;
        indexByXY = new int[xDim * yDim];
        int[] tmpXYByIndex = new int[xDim * yDim];
        for (int x = 0; x < xDim; x++) {
            for (int y = 0; y < yDim; y++) {
                final int xy = x * yDim + y;
                if (!Double.isNaN(map[x][y])) {
                    tmpXYByIndex[index] = xy;
                    indexByXY[xy] = index;
                    index++;
                } else
                    indexByXY[xy] = -1;
            }
        }
        xyByIndex = new int[index];
        xyByIndex = tmpXYByIndex;
        return index;
    }

    public double getCTMCProbability(double[] start, double[] stop, double time) {
        int startIndex = getIndex((int) start[0], (int) start[1]);
        int stopIndex = getIndex((int) stop[0], (int) stop[1]);

//		startIndex = stopIndex = 0;

        if (startIndex == -1 || stopIndex == -1)
            return 0;
//		System.err.println("Indices: "+startIndex+" -> "+stopIndex);
//		System.err.println("Time = "+time);
        double prob = matrixExp.getExponentialEntry(startIndex, stopIndex, time);
//		System.err.println("prob = "+prob);
        return prob;
    }


    public SparseMatrixExponential getMatrix() {
        return matrixExp;
    }

    public boolean isValidPoint(int x, int y) {
        return (getIndex(x, y) != -1);
    }

    public int getIndex(int x, int y) {
        if (x < 0 || x >= xDim || y < 0 || y >= yDim)
            return -1;
        return indexByXY[x * yDim + y];
    }

    public int getX(int index) {
        return xyByIndex[index] / yDim;
    }

    public int getY(int index) {
        final int xy = xyByIndex[index];
        final int dim = yDim;
        final int mod = xy % dim;
        return mod;
    }

    private void setUpRandomWalkMatrix() {

        List<GraphEdge> edgeList = new ArrayList<GraphEdge>();
        for (int index = 0; index < order; index++) {
            int x = getX(index);
            int y = getY(index);
            int dest;

            dest = getIndex(x - 1, y - 1);
            if (dest != -1)
                edgeList.add(new GraphEdge(
                        index, dest, getWeight(x, y, x - 1, y - 1)
                ));

            dest = getIndex(x - 1, y);
            if (dest != -1)
                edgeList.add(new GraphEdge(
                        index, dest, getWeight(x, y, x - 1, y)
                ));

            dest = getIndex(x - 1, y + 1);
            if (dest != -1)
                edgeList.add(new GraphEdge(
                        index, dest, getWeight(x, y, x - 1, y + 1)
                ));

            dest = getIndex(x, y - 1);
            if (dest != -1)
                edgeList.add(new GraphEdge(
                        index, dest, getWeight(x, y, x, y - 1)
                ));

            dest = getIndex(x, y + 1);
            if (dest != -1)
                edgeList.add(new GraphEdge(
                        index, dest, getWeight(x, y, x, y + 1)
                ));

            dest = getIndex(x + 1, y - 1);
            if (dest != -1)
                edgeList.add(new GraphEdge(
                        index, dest, getWeight(x, y, x + 1, y - 1)
                ));

            dest = getIndex(x + 1, y);
            if (dest != -1)
                edgeList.add(new GraphEdge(
                        index, dest, getWeight(x, y, x + 1, y)
                ));

            dest = getIndex(x + 1, y + 1);
            if (dest != -1)
                edgeList.add(new GraphEdge(
                        index, dest, getWeight(x, y, x + 1, y + 1)
                ));


        }


        double[] rowTotal = new double[order];
        nonZeroElements = edgeList.size() + order;  // don't forget the diagonals

        matrixExp = new SparseMatrixExponential(order, nonZeroElements);

//		int index = 0;
        for (GraphEdge edge : edgeList) {
            int start = edge.getStart();
            int stop = edge.getStop();
            double weight = edge.getWeight();
            rowTotal[start] -= weight;
//			bandMatrix.set(start, stop, weight);
            matrixExp.addEntry(start, stop, weight);
        }

        double norm = 0;
        for (int i = 0; i < order; i++) {
            if (-2 * rowTotal[i] > norm)
                norm = -2 * rowTotal[i];
//			bandMatrix.set(i, i, rowTotal[i]);
            matrixExp.addEntry(i, i, rowTotal[i]);

        }

        matrixExp.setNorm(norm);
    }


    public void doStuff() {
        System.err.println(matrixExp.getExponentialEntry(0, 1, 1000));
        System.err.println(matrixExp.getExponentialEntry(0, 0, 1000));


    }

    public static void writeEigenvectors(String file, double[][] mat) throws IOException {

        PrintWriter writer;
        if (file.endsWith("gz"))
            writer = new PrintWriter(
                    new OutputStreamWriter
                            (new GZIPOutputStream
                                    (new FileOutputStream(file))));
        else
            writer = new PrintWriter(new FileWriter(file));

        final int dim = mat.length;
        writer.println(dim);
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++)
                writer.println(mat[i][j]);
        }
        writer.close();

    }

    public static void writeEigenvalues(String file, double[] vec) throws IOException {

        PrintWriter writer;
        if (file.endsWith("gz"))
            writer = new PrintWriter(
                    new OutputStreamWriter
                            (new GZIPOutputStream
                                    (new FileOutputStream(file))));
        else
            writer = new PrintWriter(new FileWriter(file));

        writer.println(vec.length);
        for (double d : vec)
            writer.println(d);
        writer.close();

    }

    public static int extractInt(String line) {
        StringTokenizer st = new StringTokenizer(line);
        st.nextToken();
        return Integer.parseInt(st.nextToken());
    }

    public static final String defaultInvalidString = "*";

    public static double[][] readGRASSAscii(String file) throws IOException {
        return readGRASSAscii(file, defaultInvalidString);
    }

    public static double[][] readGRASSAscii(String file, String invalidString) throws IOException {
        double[][] map;
        BufferedReader reader = getReader(file);
        String northLine = reader.readLine();
        String southLine = reader.readLine();
        String eastLine = reader.readLine();
        String westLine = reader.readLine();
        int numberRows = extractInt(reader.readLine());
        int numberCols = extractInt(reader.readLine());
        map = new double[numberRows][numberCols];
        for (int i = 0; i < numberRows; i++) {
            String line = reader.readLine();
            StringTokenizer st = new StringTokenizer(line);
            for (int j = 0; j < numberCols; j++) {
                String item = st.nextToken();
                if (item.compareTo(invalidString) == 0) {
                    map[i][j] = Double.NaN;
                } else {
                    map[i][j] = Double.parseDouble(item);
                }
            }
        }
        return map;
    }

    public String toString() {
        return toCartogram();
    }

    public String toCartogram() {
        StringBuffer sb = new StringBuffer();
        for (int i = 0; i < xDim; i++) {
            for (int j = 0; j < yDim; j++) {
                sb.append(String.format(FORMAT, map[i][j]));
                if (j < (yDim - 1))
                    sb.append(" ");
            }
            sb.append("\n");
        }
        return sb.toString();
    }

    public static BufferedReader getReader(String file) throws IOException {
        BufferedReader reader;
        if (file.endsWith("gz"))
            reader = new BufferedReader
                    (new InputStreamReader
                            (new GZIPInputStream
                                    (new FileInputStream(file))));
        else
            reader = new BufferedReader(new FileReader(file));
        return reader;
    }

    public static BufferedWriter getWriter(String file) throws IOException {
        BufferedWriter writer;
        if (file.endsWith("gz"))
            writer = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(file))));
        else
            writer = new BufferedWriter(new FileWriter(file));
        return writer;
    }

    public static double[] readEigenvalues(String file) throws IOException {

        BufferedReader reader = getReader(file);

        final int dim = Integer.parseInt(reader.readLine());
        double[] vec = new double[dim];
        for (int i = 0; i < dim; i++)
            vec[i] = Double.parseDouble(reader.readLine());
        reader.close();

        return vec;

    }

    public static double[][] readEigenvectors(String file) throws IOException {

        BufferedReader reader = getReader(file);

        final int dim = Integer.parseInt(reader.readLine());
        double[][] mat = new double[dim][dim];
        for (int i = 0; i < dim; i++) {
            for (int j = 0; j < dim; j++)
                mat[i][j] = Double.parseDouble(reader.readLine());
        }
        reader.close();

        return mat;

    }

    /* Used for converting GRASS files to PathFinding XML */

    public static final String TERRAIN = "terrain";
    public static final String TYPE_TYPE = "tileTypes";
    public static final String TYPE = "type";

    public static String mapToXMLString(double[][] map, int numLevels) {
        XMLOutputter outputter = new XMLOutputter(org.jdom.output.Format.getPrettyFormat());
        return outputter.outputString(mapToXML(map, numLevels));
//				root.toString();
    }

//	static


    public static Element mapToXML(double[][] map, int numLevels) {

        Element root = new Element(TERRAIN);
        Element tileTypes = new Element(TYPE_TYPE);
        Element size = new Element("size");
        Element content = new Element("content");
        root.addContent(tileTypes);
        root.addContent(size);
        root.addContent(content);

        // Make sizes
        size.addContent(new Element("rows").addContent(Integer.toString(map.length)));
        size.addContent(new Element("columns").addContent(Integer.toString(map[0].length)));
        size.addContent(new Element("height").addContent(Integer.toString(map.length)));
        size.addContent(new Element("width").addContent(Integer.toString(map[0].length)));

        // Make cutpoints
        double min = Double.MAX_VALUE;
        double max = Double.MIN_VALUE;

        for (int i = 0; i < map.length; i++) {
            for (int j = 0; j < map[0].length; j++) {
                if (!Double.isNaN(map[i][j])) {
//					System.err.println("die");
//					System.exit(-1);
//				}
                    if (map[i][j] < min)
                        min = map[i][j];
                    if (map[i][j] > max)
                        max = map[i][j];
                }
            }
        }

//		System.err.println("min = ");

        double[] cuts = new double[numLevels];
//		cuts[numLevels-1] = max;
        double range = max - min;
        double delta = range / numLevels;
        for (int i = 0; i < numLevels; i++)
            cuts[i] = min + delta * (i + 1);

        int deltaColor = 255 / numLevels;

        // Make tile types
        double average = min + delta / 2;
        int grey = 255 - deltaColor;
        for (int i = 0; i < numLevels; i++) {
            Element type = new Element("type");
            type.addContent(new Element("name").addContent("L" + Integer.toString(i)));
            type.addContent(new Element("cost").addContent(Integer.toString((int) average)));
            Element color = new Element("color");
            color.setAttribute("r", Integer.toString(grey));
            color.setAttribute("g", Integer.toString(grey));
            color.setAttribute("b", Integer.toString(grey));
            type.addContent(color);
            tileTypes.addContent(type);
            average += delta;
            grey -= deltaColor;
        }
        Element blocked = new Element("type");
        blocked.addContent(new Element("name").addContent("B"));
        blocked.addContent(new Element("blocked"));
        Element color = new Element("color");
        color.setAttribute("r", Integer.toString(255));
        color.setAttribute("g", Integer.toString(255));
        color.setAttribute("b", Integer.toString(255));
        blocked.addContent(color);
        tileTypes.addContent(blocked);

        int count = 0;

        // Make content
        content.addContent(new Element("default").addContent("B"));
        StringBuffer ascii = new StringBuffer();
        boolean newLine = true;

        for (int r = 0; r < map.length; r++) {
            for (int c = 0; c < map[0].length; c++) {
                if (!Double.isNaN(map[r][c])) {
                    count++;
                    double value = map[r][c];
//				if( value != Double.NaN ) {
                    int cut = 0;
                    try {
                        while (value > cuts[cut])
                            cut++;
                    } catch (Exception e) {
                        System.err.println("Error: " + e);
                        System.err.println("value = " + value);
                        System.err.println("cuts = " + new Vector(cuts));
                        System.err.println("min = " + min);
                        System.err.println("max = " + max);
                        System.exit(-1);
                    }
                    Element cell = new Element("column");
                    cell.setAttribute("r", Integer.toString(r));
                    cell.setAttribute("c", Integer.toString(c));
                    cell.setAttribute("length", "1");
                    cell.addContent("L" + Integer.toString(cut));
//					if( count < 10)
                    content.addContent(cell);
                    ascii.append((cut + 1) + " ");

//				}  else {
//					System.err.println("about time");
//					System.exit(-1);
                } else {
                    ascii.append("NA ");
                }

            }
            ascii.append("\n");
        }

        System.out.println(ascii);

        return root;
    }


    public static void writeXML(String file) throws IOException {

    }

//	private double getProbability(int I, int J, double time) {
//		double probability = 0;
//		for (int k = 0; k < order; k++) {
//			probability += eVec[I][k] * Math.exp(time * eVal[k]) * eVec[J][k];
//		}
//		return probability;
//
//	}


    public class GraphEdge {
        private int start;
        private int stop;
        private double weight;

        public GraphEdge(int start, int stop, double weight) {
            this.start = start;
            this.stop = stop;
            this.weight = weight;
        }

        public int getStart() {
            return start;
        }

        public void setStart(int start) {
            this.start = start;
        }

        public int getStop() {
            return stop;
        }

        public void setStop(int stop) {
            this.stop = stop;
        }

        public double getWeight() {
            return weight;
        }

        public void setWeight(double weight) {
            this.weight = weight;
        }
    }

    public double getWeight(int x0, int y0, int x1, int y1) {
        return 1.0 / (100 * Math.abs(map[x0][y0] - map[x1][y1]) + 1.0);
    }


    public static void main(String[] args) {
        double[][] map = null;
        try {
//		try {
            map = readGRASSAscii(args[1]);
//		} catch (IOException e) {
//			e.printStackTrace();
//		}
            System.err.println("Read in GRASS file.");
            System.err.println("Writing XML file.");
            String xml = mapToXMLString(map, Integer.parseInt(args[0]));
            PrintWriter writer = new PrintWriter(new FileWriter(args[2]));
            writer.println(xml);
            writer.close();
        } catch (Exception e) {
            System.err.println("Command-line error.");
            System.err.println("USAGE: program_name <# of levels> <GRASS file> <XML file>");
        }

//		new TopographicalMap(map).doStuff();
    }

    private double[][] map;
    private int[] indexByXY;
    private int[] xyByIndex;
    SparseMatrixExponential matrixExp;

    private int xDim;
    private int yDim;
    private int order;
    private int nonZeroElements;

    public static double[][] sillyMap = {
            {0.0, 100.0, 0.0, 100.0, 0.0},
            {0.0, 0.0, 0.0, 100.0, 0.0},
            {100.0, 0.0, 100.0, 100.0, 0.0},
            {100.0, 0.0, 100.0, 100.0, 0.0},
            {100.0, 0.0, 0.0, 0.0, 0.0}
    };

}

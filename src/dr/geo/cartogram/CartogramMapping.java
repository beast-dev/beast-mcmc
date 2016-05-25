/*
 * CartogramMapping.java
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

package dr.geo.cartogram;

import dr.util.FileHelpers;

import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.io.*;
import java.util.StringTokenizer;
import java.util.zip.GZIPInputStream;

/**
 * Provides basic functionality for output from Mark Newman's 'cartogram'
 * - Reads resulting flat text file
 * - Uses bilinear interpolation to map (real space x, real space y) -> (cartogram x, cartogram y)
 *
 * @author Marc A. Suchard
 */
public class CartogramMapping {

    public static final String CARTOGRAM_OUTPUT = "cartogram";


    public CartogramMapping(int gridXSize, int gridYSize, Rectangle2D boundingBox) {

        this.gridXSize = gridXSize;
        this.gridYSize = gridYSize;
        this.boundingBox = boundingBox;

        gridPt = new Point2D[gridXSize + 1][gridYSize + 1];

        dX = (boundingBox.getMaxX() - boundingBox.getMinX()) / gridXSize;
        dY = (boundingBox.getMaxY() - boundingBox.getMinY()) / gridYSize;

    }

    public double getAverageDensity() {
        return averageDensity;
    }

    public void setAverageDensity(double d) {
        averageDensity = d;
    }

    public void readCartogramOutput(String fileName) throws IOException {
        readCartogramOutput(FileHelpers.getFile(fileName));
    }

    public void readCartogramOutput(File file) throws IOException {
        BufferedReader reader;
        if (file.getName().endsWith("gz"))
            reader = new BufferedReader
                    (new InputStreamReader
                            (new GZIPInputStream
                                    (new FileInputStream(file))));
        else
            reader = new BufferedReader(new FileReader(file));

        // Newman stores values in column-major
        for (int y = 0; y <= gridYSize; y++) {
            for (int x = 0; x <= gridXSize; x++) {
                String line = reader.readLine();
                if (line == null)
                    throw new IOException("Premature end of file in '" + file + "'");
                StringTokenizer st = new StringTokenizer(line);
                try {
                    double mappedX = Double.parseDouble(st.nextToken());
                    double mappedY = Double.parseDouble(st.nextToken());
                    gridPt[x][y] = new Point2D.Double(mappedX, mappedY);
                } catch (NumberFormatException e) {
                    throw new IOException("Unable to parse line: " + line + " in '" + file + "'");
                }
            }
        }
        reader.close();
        loaded = true;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("grid: [x=" + gridXSize + ",y=" + gridYSize + "], boundingBox: " + boundingBox.toString());
        return sb.toString();
    }

    public Point2D map(Point2D inPt) {

        if (!boundingBox.contains(inPt))
            return null;

        if (!loaded)
            return inPt;

        final double offsetX = (inPt.getX() - boundingBox.getMinX()) / dX;
        final double offsetY = (inPt.getY() - boundingBox.getMinY()) / dY;

        final int iX = (int) offsetX;
        final int iY = (int) offsetY;

        final double rX = offsetX - iX;
        final double rY = offsetY - iY;

        assert (rX >= 0 && rY < 1.0 && rY >= 0 && rY < 1.0);

        final Point2D gridiXiY = gridPt[iX][iY];
        final Point2D gridiX1iY = gridPt[iX + 1][iY];
        final Point2D gridiXiY1 = gridPt[iX][iY + 1];
        final Point2D gridiX1iY1 = gridPt[iX + 1][iY + 1];

        final double outX = (1 - rX) * (1 - rY) * gridiXiY.getX() +
                rX * (1 - rY) * gridiX1iY.getX() +
                (1 - rX) * rY * gridiXiY1.getX() +
                rX * rY * gridiX1iY1.getX();


        final double outY = (1 - rX) * (1 - rY) * gridiXiY.getY() +
                rX * (1 - rY) * gridiX1iY.getY() +
                (1 - rX) * rY * gridiXiY1.getY() +
                rX * rY * gridiX1iY1.getY();

        return new Point2D.Double(outX, outY);

    }

    /*
     * This should replicate the behavoir of Newman's 'interp' program
     */
    public static void main(String[] args) {

        int gridXSize = Integer.parseInt(args[0]);
        int gridYSize = Integer.parseInt(args[1]);
        String fileName = args[2];

        CartogramMapping mapping = new CartogramMapping(gridXSize, gridYSize,
                new Rectangle2D.Double(0, 0, gridXSize, gridYSize));
        try {
            mapping.readCartogramOutput(fileName);

            BufferedReader reader = new BufferedReader(new InputStreamReader(System.in));
            while (true) {

                String line = reader.readLine();
                if (line == null) {
                    // end of file
                    break;
                } else if (line.length() == 0) {
                    break;
                } else {
                    // not a blank line
                    StringTokenizer st = new StringTokenizer(line);
                    Point2D inPt = new Point2D.Double(
                            Double.parseDouble(st.nextToken()),
                            Double.parseDouble(st.nextToken())
                    );
                    System.out.println(mapping.map(inPt));
                }

            }
        } catch (IOException e) {
            System.err.println(e);
            System.exit(-1);
        }
    }

    private Rectangle2D boundingBox;
    private int gridXSize, gridYSize;
    private double dX, dY;
    private Point2D[][] gridPt;
    private boolean loaded = false;
    private double averageDensity = 1.0;
}

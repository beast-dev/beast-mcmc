/*
 * GTOPO30Tile.java
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

package dr.geo;

import dr.app.gui.ColorFunction;

import javax.imageio.ImageIO;
import javax.imageio.stream.ImageInputStream;
import javax.swing.*;
import java.awt.*;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.*;
import java.nio.ByteOrder;
import java.util.HashMap;
import java.util.Map;
import java.util.StringTokenizer;

/**
 * @author Alexei Drummond
 * @author Remco Bouckaert
 */
public class GTOPO30Tile extends JComponent {

    static final int NROWS = 6000;
    static final int NCOLS = 4800;
    static final int NODATA = -9999;

    short[][] height;
    BufferedImage image;
    ColorFunction colorFunction;

    double ulxmap, ulymap, xdim, ydim;

    Rectangle2D bounds;

    boolean createImage = true;


    public GTOPO30Tile(String filename, ColorFunction function) throws IOException {

        String headerFileName = filename.substring(0, filename.lastIndexOf('.')) + ".HDR";

        File headerFile = new File(headerFileName);

        if (headerFile.isFile()) {
            try {
                readHeader(headerFile);
            } catch (FileNotFoundException e) {
                System.err.println("No header file named " + headerFile.toString() + " found");
            }
        }

        height = new short[NROWS][NCOLS];

        colorFunction = function;

        read(filename);
    }

    private void readHeader(File headerFile) throws IOException {

        BufferedReader reader = new BufferedReader(new FileReader(headerFile));
        Map<String, String> properties = new HashMap<String, String>();

        String line = reader.readLine();
        while (line != null) {
            StringTokenizer tokens = new StringTokenizer(line);
            if (tokens.countTokens() == 2) {
                properties.put(tokens.nextToken(), tokens.nextToken());
            }
            line = reader.readLine();
        }
        ulxmap = Double.parseDouble(properties.get("ULXMAP"));
        ulymap = Double.parseDouble(properties.get("ULYMAP"));
        xdim = Double.parseDouble(properties.get("XDIM"));
        ydim = Double.parseDouble(properties.get("YDIM"));

        bounds = new Rectangle2D.Double(ulxmap - xdim / 2.0, ulymap + ydim / 2.0 - ydim * NROWS, xdim * NCOLS, ydim * NROWS);

        if (NROWS != Integer.parseInt(properties.get("NROWS"))) {
            throw new RuntimeException("Expected NROWS to be " + NROWS + " in header file " + headerFile);
        }
        if (NCOLS != Integer.parseInt(properties.get("NCOLS"))) {
            throw new RuntimeException("Expected NCOLS to be " + NCOLS + " in header file " + headerFile);
        }
    }

    void read(String sFile) {

        if (createImage) {
            image = new BufferedImage(NCOLS, NROWS, BufferedImage.TYPE_INT_RGB);
        }

        try {
            ImageInputStream iis = ImageIO.createImageInputStream(new File(sFile));
            iis.setByteOrder(ByteOrder.BIG_ENDIAN);

            for (int y = 0; y < NROWS; y++) {
                iis.readFully(height[y], 0, NCOLS);

                if (createImage) {

                    for (int x = 0; x < NCOLS; x++) {
                        int color = height[y][x];//iis.readShort();
                        if (color == -9999) {
                            image.setRGB(x, y, Color.blue.darker().getRGB());
                        } else {
                            image.setRGB(x, y, colorFunction.getColor((float) color).getRGB());
                        }
                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        System.err.println();
    } // read

    protected void paintComponent(Graphics g) {
        int nW = getWidth();
        int nH = getHeight();
        g.drawImage(image, 0, 0, nW, nH, 0, 0, NCOLS, NROWS, null);
    }

    /**
     * @param y latitudinal pixel, increasing south from north edge
     * @param x longitudinal pixel index, increasing east from west edge
     * @return the height of the pixel at y,x
     */
    public short getHeight(int y, int x) {
        return height[y][x];
    }

    boolean contains(double longitude, double latitude) {
        return bounds.contains(longitude, latitude);
    }

    public int getMinLongitude() {
        return (int) Math.round(bounds.getMinX());
    }

    public int getMaxLongitude() {
        return (int) Math.round(bounds.getMaxX());
    }

    public int getMinLatitude() {
        return (int) Math.round(bounds.getMinY());
    }

    public int getMaxLatitude() {
        return (int) Math.round(bounds.getMaxY());
    }
}

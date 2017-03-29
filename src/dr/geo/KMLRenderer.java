/*
 * KMLRenderer.java
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

import java.awt.*;
import java.awt.geom.AffineTransform;
import java.awt.geom.GeneralPath;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.awt.image.Raster;
import java.util.ArrayList;

/**
 * @author Alexei Drummond
 */
public class KMLRenderer implements Lattice {

    BufferedImage image;
    int[][] lattice;
    Rectangle2D bounds;

    java.util.List<AbstractPolygon2D> polygons;
    java.util.List<Shape> shapes;

    ViewTransform viewTransform;

    Color shapeColor;
    Color background;

    public KMLRenderer(String kmlFileName, Color shapeColor, Color background) {

        polygons = Polygon2D.readKMLFile(kmlFileName);

        this.shapeColor = shapeColor;
        this.background = background;

        System.out.println("Read " + polygons.size() + " polygons");

        System.out.println("Converting polygons to shapes");

        double minX = Double.MAX_VALUE, maxX = -Double.MAX_VALUE, minY = Double.MAX_VALUE, maxY = -Double.MAX_VALUE;

        shapes = new ArrayList<Shape>();

        for (AbstractPolygon2D p : polygons) {

            Shape shape = p.getShape();

            bounds = shape.getBounds();

            if (bounds.getMinX() < minX) minX = bounds.getMinX();
            if (bounds.getMaxX() > maxX) maxX = bounds.getMaxX();
            if (bounds.getMinY() < minY) minY = bounds.getMinY();
            if (bounds.getMaxY() > maxY) maxY = bounds.getMaxY();

            shapes.add(shape);
            System.out.print(".");
            System.out.flush();
        }

        bounds = new Rectangle2D.Double(minX, minY, maxX - minX, maxY - minY);

        System.out.println();

    }

    public Rectangle2D getBounds() {
        return bounds;
    }

    public BufferedImage render(int size) {
        int width;
        int height;
        if (bounds.getHeight() > bounds.getWidth()) {
            height = size;
            width = (int) (height * bounds.getWidth() / bounds.getHeight());
        } else {
            width = size;
            height = (int) (width * bounds.getHeight() / bounds.getWidth());
        }
        return render(width, height);
    }

    public BufferedImage render(int width, int height) {
        image = new BufferedImage(width, height, BufferedImage.TYPE_INT_ARGB);
        render(image);

        Raster raster = image.getData();

        lattice = new int[width][height];
        int[] pixel = new int[4];
        for (int i = 0; i < width; i++) {
            for (int j = 0; j < height; j++) {
                raster.getPixel(i, j, pixel);
                if (colorDistanceSquared(pixel, shapeColor) < colorDistanceSquared(pixel, background)) {
                    lattice[i][j] = 1;
                } else {
                    lattice[i][j] = 0;
                }
            }
        }

        return image;
    }

    private double colorDistanceSquared(int[] pixel, Color color) {

        double[] argb = new double[4];

        argb[0] = Math.abs(pixel[0] - color.getAlpha());
        argb[1] = Math.abs(pixel[1] - color.getRed());
        argb[2] = Math.abs(pixel[2] - color.getGreen());
        argb[3] = Math.abs(pixel[3] - color.getBlue());

        double dist = 0;
        for (double a : argb) {
            dist += a * a;
        }
        return dist;
    }

    public void render(BufferedImage image) {
        Graphics2D g2d = image.createGraphics();

        g2d.setColor(background);
        g2d.fillRect(0, 0, image.getWidth(), image.getHeight());

        viewTransform = new ViewTransform(bounds, image.getWidth(), image.getHeight());

        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        g2d.setColor(shapeColor);
        AffineTransform transform = viewTransform.getTransform();
        for (Shape s : shapes) {
            GeneralPath path = new GeneralPath(s);
            path.transform(transform);
            g2d.fill(path);
        }
    }

    public int latticeWidth() {
        return lattice.length;
    }

    public int latticeHeight() {
        return lattice[0].length;
    }

    public void setState(int i, int j, int state) {
        lattice[i][j] = state;
    }

    public int getState(int i, int j) {
        return lattice[i][j];
    }

    public void paintLattice(Graphics g) {
        g.drawImage(image, 0, 0, null);
    }

    public void setBounds(Rectangle2D bounds) {
        this.bounds = bounds;
    }
}

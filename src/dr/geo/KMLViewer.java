/*
 * VisualizeKMLBrownianBridge.java
 *
 * Copyright (C) 2002-2010 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.geo;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.*;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexei Drummond
 */
public class KMLViewer extends JComponent {

    List<Polygon2D> polygons;

    List<Shape> shapes;

    Paint shapeColor = Color.GRAY;

    Point2D topLeft;
    Point2D bottomRight;

    double scaleX;
    double scaleY;


    public KMLViewer(String kmlFileName) {

        polygons = Polygon2D.readKMLFile(kmlFileName);

        System.out.println("Read " + polygons.size() + " polygons");

        System.out.println("Converting polygons to shapes");

        double minX = Double.MAX_VALUE, maxX = -Double.MAX_VALUE, minY = Double.MAX_VALUE, maxY = -Double.MAX_VALUE;

        shapes = new ArrayList<Shape>();

        for (Polygon2D p : polygons) {

            Shape shape = getShape(p);

            Rectangle2D bounds = shape.getBounds();

            if (bounds.getMinX() < minX) minX = bounds.getMinX();
            if (bounds.getMaxX() > maxX) maxX = bounds.getMaxX();
            if (bounds.getMinY() < minY) minY = bounds.getMinY();
            if (bounds.getMaxY() > maxY) maxY = bounds.getMaxY();

            shapes.add(shape);
            System.out.print(".");
            System.out.flush();
        }

        topLeft = new Point2D.Double(minX, minY);
        bottomRight = new Point2D.Double(maxX, maxY);

        System.out.println();
    }

    public Rectangle2D getViewport() {
        return new Rectangle2D.Double(topLeft.getX(), topLeft.getY(),
                bottomRight.getX() - topLeft.getX(), bottomRight.getY() - topLeft.getY());
    }

    void computeScales() {
        scaleX = getWidth() / (bottomRight.getX() - topLeft.getX());
        scaleY = getHeight() / (bottomRight.getY() - topLeft.getY());
    }

    public void paintComponent(Graphics g) {

        System.out.println("entering paintComponent()");

        computeScales();

        Graphics2D g2d = (Graphics2D) g;

        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        g2d.setStroke(new BasicStroke(1.5f));

        System.out.println("Painting shapes");

        AffineTransform transform = getFullTransform();
        for (Shape s : shapes) {
            System.out.print(".");
            System.out.flush();
            GeneralPath path = new GeneralPath(s);
            path.transform(transform);
            g2d.setPaint(shapeColor);
            g2d.fill(path);
        }

        g2d.setColor(Color.yellow);
        drawGrid(5, 5, g2d, transform);

        System.out.println("leaving paintComponent()");

    }

    private void drawGrid(int dLat, int dLong, Graphics2D g2d, AffineTransform transform) {

        for (double longitude = -180; longitude < 180; longitude += dLong) {
            Line2D line = new Line2D.Double(longitude, -90, longitude, 90);
            g2d.draw(transform.createTransformedShape(line));
        }
        for (double lat = -90; lat < 90; lat += dLat) {
            Line2D line = new Line2D.Double(-180, lat, 180, lat);
            g2d.draw(transform.createTransformedShape(line));
        }
    }

    Shape getShape(Polygon2D poly) {
        GeneralPath path = new GeneralPath();

        List<Point2D> points = poly.point2Ds;
        path.moveTo((float) points.get(0).getX(), (float) points.get(0).getY());

        System.out.println("x=" + points.get(0).getX() + ", y=" + points.get(0).getY());

        for (int i = 1; i < points.size(); i++) {
            path.lineTo((float) points.get(i).getX(), (float) points.get(i).getY());
        }
        path.closePath();
        return path;
    }

    AffineTransform getTranslate() {
        return AffineTransform.getTranslateInstance(-topLeft.getX(), -bottomRight.getY());
    }

    AffineTransform getScale() {
        return AffineTransform.getScaleInstance(scaleX, -scaleY);
    }

    AffineTransform getFullTransform() {
        AffineTransform transform = getScale();
        transform.concatenate(getTranslate());
        return transform;
    }

    public static void main(String[] args) {

        String filename = args[0];

        JFrame frame = new JFrame("KMLViewer - " + filename);

        KMLViewer viewer = new KMLViewer(filename);
        Rectangle2D viewport = viewer.getViewport();

        frame.getContentPane().add(BorderLayout.CENTER, viewer);

        int width;
        int height;
        if (viewport.getHeight() > viewport.getWidth()) {
            height = 900;
            width = (int) (height * viewport.getWidth() / viewport.getHeight());
        } else {
            width = 900;
            height = (int) (width * viewport.getHeight() / viewport.getWidth());
        }
        System.out.println("Height = " + height + ", Width = " + width);

        frame.setSize(width, height);
        frame.setVisible(true);
    }

}
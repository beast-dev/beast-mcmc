/*
 * DensityRectangle.java
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

package dr.inference.trace;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexei Drummond
 * @version $Id: DensityRectangle.java,v 1.4 2006/08/26 19:24:41 rambaut Exp $
 */
public class DensityRectangle {

    Rectangle2D rectangle;
    List<Point2D> points;
    int pointCount;
    double area;

    static int THRESHOLD = 80;

    DensityRectangle[] quadrants = null;


    public DensityRectangle(double x, double y, double w, double h) {

        rectangle = new Rectangle2D.Double(x, y, w, h);
        area = w * h;
        points = new ArrayList<Point2D>();
    }

    public boolean contains(Point2D point) {

        return rectangle.contains(point);
    }

    public boolean contains(double x, double y) {

        return rectangle.contains(x, y);
    }

    public void addPoint(Point2D point) {
        if (contains(point)) {

            pointCount += 1;

            if (quadrants == null) {

                points.add(point);

                if (points.size() > THRESHOLD) {

                    // partition into quadrants
                    quadrants = new DensityRectangle[4];

                    double x = rectangle.getX();
                    double y = rectangle.getY();
                    double w = rectangle.getWidth() / 2.0;
                    double h = rectangle.getHeight() / 2.0;

                    quadrants[0] = new DensityRectangle(x, y, w, h);
                    quadrants[1] = new DensityRectangle(x + w, y, w, h);
                    quadrants[2] = new DensityRectangle(x, y + h, w, h);
                    quadrants[3] = new DensityRectangle(x + w, y + h, w, h);

                    for (Point2D p : points) {
                        quadrantPoint(p);
                    }
                    points.clear();
                    points = null;
                }
            } else {

                quadrantPoint(point);

            }

        } else {
            // ignore point!
        }
    }

    private void quadrantPoint(Point2D point) {

        boolean found = false;
        for (int i = 0; i < quadrants.length && !found; i++) {
            if (quadrants[i].contains(point)) {
                quadrants[i].addPoint(point);
                found = true;
            }
        }
        if (!found) {
            throw new RuntimeException("Failed to quadrantize point!");
        }
    }

    public final double density() {
        return getPointCount() / getArea();
    }

    public final double maximumDensity() {
        if (quadrants == null) return density();

        double maxDensity = quadrants[0].maximumDensity();
        for (int i = 1; i < quadrants.length; i++) {
            double d = quadrants[i].maximumDensity();
            if (d > maxDensity) maxDensity = d;
        }
        return maxDensity;
    }

    public final int getPointCount() {
        return pointCount;
    }

    public final double getArea() {
        return area;
    }

    public final void paint(Graphics2D g2d, Rectangle2D targetRectangle) {

        paint(g2d, targetRectangle, maximumDensity());
    }

    private void paint(Graphics2D g2d, Rectangle2D targetRectangle, double maxDensity) {

        double density = density();

        float red = (float) (density / maxDensity);
        float blue = 1.0f - red;
        float green = Math.min(red, blue); //(0.5f - Math.abs(red-0.5f));

        Color color = new Color(red, green, blue);
        g2d.setColor(color);
        g2d.fill(targetRectangle);

        if (quadrants != null) {

            double x = targetRectangle.getX();
            double y = targetRectangle.getY();
            double w = targetRectangle.getWidth() / 2.0;
            double h = targetRectangle.getHeight() / 2.0;

            Rectangle2D rect[] = new Rectangle2D[4];

            rect[0] = new Rectangle2D.Double(x, y, w, h);
            rect[1] = new Rectangle2D.Double(x + w, y, w, h);
            rect[2] = new Rectangle2D.Double(x, y + h, w, h);
            rect[3] = new Rectangle2D.Double(x + w, y + h, w, h);

            for (int i = 0; i < quadrants.length; i++) {
                quadrants[i].paint(g2d, rect[i], maxDensity);
            }
        }
    }

    public static void main(String[] args) throws IOException, TraceException {

        int burnin = 0;

        System.out.println("usage: densityPlot <logFileName> <xIndex> <yIndex> <maxX> <maxY>");

        String logFileName = args[0];

        int[] indices = new int[2];

        for (int i = 0; i < 2; i++) {
            indices[i] = Integer.parseInt(args[i + 1]);
        }

        double width = Double.parseDouble(args[3]);
        double height = Double.parseDouble(args[4]);

        System.out.println("  logFileName=" + logFileName);
        System.out.println("  xIndex=" + indices[0]);
        System.out.println("  yIndex=" + indices[1]);
        System.out.println("  maxX=" + width);
        System.out.println("  maxY=" + height);

        File file = new File(logFileName);

        LogFileTraces traces = new LogFileTraces(logFileName, file);
        traces.loadTraces();
        traces.setBurnIn(burnin);

        ArrayList<ArrayList> values = new ArrayList<ArrayList>();
        for (int j = 0; j < 2; j++) {
            values.add(new ArrayList(traces.getValues(indices[j])));
        }
        int totalSize = values.get(0).size();
        System.out.println("total samples = " + totalSize);

        final DensityRectangle densityRectangle = new DensityRectangle(0, 0, width, height);

        for (int i = 0; i < values.get(0).size(); i++) {
            densityRectangle.addPoint(new Point2D.Double((Double) values.get(0).get(i), (Double) values.get(1).get(i)));
        }

        System.out.println("posterior prob visualized = " + ((double) densityRectangle.getPointCount()) / totalSize);


        JFrame frame = new JFrame();
        frame.setSize(400, 400);
        frame.getContentPane().add(
                new JComponent() {
                    /**
                     *
                     */
                    private static final long serialVersionUID = -4045142094590087910L;

                    public void paintComponent(Graphics g) {

                        Rectangle2D bounds = getBounds();

                        densityRectangle.paint((Graphics2D) g, bounds);

                    }
                }, BorderLayout.CENTER
        );
        frame.setVisible(true);
    }

}

/*
 * VisualizeSpaceTimeSim2D.java
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

import dr.math.distributions.MultivariateNormalDistribution;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Random;

/**
 * @author Alexei Drummond
 */
public class VisualizeSpaceTimeSim2D extends JComponent {

    MultivariateNormalDistribution mnd;

    // Starting point
    SpaceTime start;

    SpaceTimeRejector rejector;

    List<Shape> shapes;

    Paint shapeColor = Color.GRAY;

    Point2D topLeft;
    Point2D bottomRight;

    double scaleX;
    double scaleY;

    SpaceTimeSimulator sim;
    Random random = new Random();

    boolean drawPath = false;
    boolean drawPoints = false;
    boolean drawFinalPoint = true;

    boolean paintDensity = true;

    int latticeWidth = 100;
    int latticeHeight = 100;

    public static final int POINT_SIZE = 3;
    public static final int steps = 100;
    public static final double dt = 0.01;

    public VisualizeSpaceTimeSim2D(final Rectangle2D bounds, final List<Shape> obstructions) {

        start = new SpaceTime(0, new double[]{0, 0});

        double x1 = bounds.getMinX();
        double y1 = bounds.getMinY();
        double x2 = bounds.getMaxX();
        double y2 = bounds.getMaxY();

        topLeft = new Point2D.Double(x1, y1);
        bottomRight = new Point2D.Double(x2, y2);

        this.shapes = obstructions;

        rejector = new SpaceTimeRejector() {

            ArrayList<Reject> rejects = new ArrayList<Reject>();

            public boolean reject(double time, double[] space) {
                Point2D p = new Point2D.Double(space[0], space[1]);

                if (!bounds.contains(p)) return true;

                for (Shape s : shapes) {
                    if (s.contains(p)) {

                        rejects.add(new Reject(0, time, space));
                        return true;
                    }
                }
                return false;
            }

            // removes all rejects
            public void reset() {
                rejects.clear();
            }

            public List<Reject> getRejects() {
                return rejects;
            }
        };

        mnd = new MultivariateNormalDistribution(new double[]{0.0}, new double[][]{{10, 0}, {0, 10}});
        sim = new SpaceTimeSimulator(mnd);

    }

    public void setShapeColor(Color c) {
        shapeColor = c;
    }

    void computeScales() {
        scaleX = getWidth() / (bottomRight.getX() - topLeft.getX());
        scaleY = getHeight() / (bottomRight.getY() - topLeft.getY());
    }

    public void paintComponent(Graphics g) {

        System.out.println("entering paintComponent()");

        computeScales();

        Graphics2D g2d = (Graphics2D) g;

        //g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        g2d.setStroke(new BasicStroke(1.5f));

        System.out.println("Painting shapes");
        for (Shape s : shapes) {
            System.out.print(".");
            System.out.flush();
            GeneralPath path = new GeneralPath(s);
            path.transform(getFullTransform());
            g2d.setPaint(shapeColor);
            g2d.fill(path);
        }

        AffineTransform transform = getFullTransform();

        List<SpaceTime> points;
        for (int r = 0; r < getTrials(); r++) {
            Color c = new Color(random.nextInt(256), random.nextInt(256), random.nextInt(256), 128);
            g.setColor(c);

            if (drawPoints || drawPath) {
                points = sim.simulatePath(start, rejector, dt, steps);
            } else {
                points = new ArrayList<SpaceTime>();
                points.add(sim.simulate(start, rejector, dt, steps));
            }

            if (drawPoints) {
                for (SpaceTime s : points) {
                    SpaceTime.paintDot(s, POINT_SIZE, transform, g2d);
                }
            }

            if (drawPath) {
                GeneralPath path = new GeneralPath();
                path.moveTo((float) points.get(0).getX(0), (float) points.get(0).getX(1));
                //System.out.println(points.get(0));
                for (int i = 1; i < points.size(); i++) {
                    path.lineTo((float) points.get(i).getX(0), (float) points.get(i).getX(1));
                    //System.out.println(points.get(i));
                }
                path.transform(getFullTransform());
                g2d.draw(path);
            }

            if (drawFinalPoint) {
                SpaceTime.paintDot(points.get(points.size() - 1), POINT_SIZE, transform, g2d);
            }

        }
        g2d.setPaint(Color.black);
        SpaceTime.paintDot(start, POINT_SIZE, transform, g2d);

        System.out.println("leaving paintComponent()");

    }

    AffineTransform getScale() {
        return AffineTransform.getScaleInstance(scaleX, scaleY);
    }

    AffineTransform getTranslate() {
        return AffineTransform.getTranslateInstance(-topLeft.getX(), -topLeft.getY());
    }

    AffineTransform getFullTransform() {
        AffineTransform transform = getScale();
        transform.concatenate(getTranslate());
        return transform;
    }

    public static void main(String[] args) {

        ArrayList<Shape> shapes = new ArrayList<Shape>();
        shapes.add(new Ellipse2D.Double(0.15, 0.15, 0.3, 0.3));
        shapes.add(new Ellipse2D.Double(0.5, 0.7, 0.2, 0.2));
        shapes.add(new Ellipse2D.Double(0.8, 0.2, 0.15, 0.15));

        Rectangle2D bounds = new Rectangle2D.Double(-0.2, -0.2, 1.4, 1.4);

        JFrame frame = new JFrame("Boulders");
        frame.getContentPane().add(BorderLayout.CENTER, new VisualizeSpaceTimeSim2D(bounds, shapes));
        frame.setSize(700, 700);
        frame.setVisible(true);
    }

    public int getTrials() {

        return 20000;
    }
}
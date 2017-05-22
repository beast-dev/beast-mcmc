/*
 * VisualizeKMLNumericalProbs.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import dr.math.distributions.MultivariateNormalDistribution;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.AffineTransform;
import java.awt.geom.GeneralPath;
import java.awt.geom.Point2D;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexei Drummond
 */
public class VisualizeKMLNumericalProbs extends JComponent {

    List<AbstractPolygon2D> polygons;

    Point2D brussels = new Point2D.Double(4.35, 50.85);
    Point2D amsterdam = new Point2D.Double(4.89, 52.37);
    Point2D berlin = new Point2D.Double(13.41, 52.52);
    Point2D rome = new Point2D.Double(12.48, 41.9);
    Point2D naples = new Point2D.Double(14.283, 40.85);

    Point2D athens = new Point2D.Double(23.72, 37.98);
    Point2D paris = new Point2D.Double(2.35, 48.86);
    Point2D montepelier = new Point2D.Double(3.88, 43.61);
    Point2D munich = new Point2D.Double(11.58, 48.14);
    Point2D bern = new Point2D.Double(7.45, 46.95);

    Point2D start, end;
    Point2D topLeft;
    Point2D bottomRight;
    List<Shape> shapes;
    SpaceTimeRejector rejector;
    MultivariateNormalDistribution D;
    NumericalSpaceTimeProbs2D probs;

    ColorFunction cf;

    double scaleX;
    double scaleY;

    boolean ABSORBING = true;

    public VisualizeKMLNumericalProbs(String kmlFileName) {

        polygons = Polygon2D.readKMLFile(kmlFileName);

        System.out.println("Read " + polygons.size() + " polygons");

        start = bern;
        end = rome;

        topLeft = new Point2D.Double(-5, 28);
        bottomRight = new Point2D.Double(25, 57);

        System.out.println("Converting polygons to shapes");
        shapes = new ArrayList<Shape>();
        for (AbstractPolygon2D p : polygons) {
            shapes.add(p.getShape());
            System.out.print(".");
            System.out.flush();
        }
        System.out.println();

        Rectangle2D bounds = new Rectangle2D.Double(topLeft.getX(), topLeft.getY(), bottomRight.getX() - topLeft.getX(), bottomRight.getY() - topLeft.getY());

        final SpaceTimeRejector boundsRejector = SpaceTimeRejector.Utils.createSimpleBounds2D(bounds);

        rejector = new SpaceTimeRejector() {

//            private boolean stop = false;

            public boolean reject(double time, double[] space) {

                if (boundsRejector.reject(time, space)) return true;

                double x = space[0];
                double y = space[1];
                for (Shape s : shapes) {
                    if (s.contains(x, y)) {
                        return true;
                    }
                }
                return false;
            }

            // removes all rejects
            public void reset() {
            }

            public List<Reject> getRejects() {
                return null;
            }
        };

        cf = new ColorFunction(
                new Color[]{Color.white, Color.blue, Color.magenta, Color.red},
                new float[]{0.0f, 0.10f, 0.20f, 1.0f}
        );

        D = new MultivariateNormalDistribution(new double[]{0.0}, new double[][]{{1, 0}, {0, 1}});

        if (rejector.reject(0.0, new double[]{start.getX(), start.getY()})) {
            throw new RuntimeException("The start position was rejected!");
        }

        //probs = new NumericalSpaceTimeProbs2D(50, 50, 50, 50, 1, bounds, D, boundsRejector);
        probs = new NumericalSpaceTimeProbs2D(50, 50, 50, 50, 1, bounds, D, rejector);

        System.out.println("Populating...");
        if (ABSORBING) {
            int successes = probs.populateAbsorbing(start, 2000000);
            System.out.println(successes + " paths simulated successfully simulated");
        } else {
            probs.populate(start, 25000, false);
        }
        System.out.println("Finished populating...");
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

        int sx = probs.x(start.getX());
        int sy = probs.y(start.getY());
        int t = 49;

        double maxCount1 = probs.maxCount[sx][sy][t];
        System.out.println("start max count = " + maxCount1);


        AffineTransform transform = getFullTransform();

        System.out.println("Painting lattice probs");
        for (int i = 0; i < probs.latticeWidth; i++) {
            for (int j = 0; j < probs.latticeHeight; j++) {

                float I1 = (float) probs.r(sx, sy, i, j, t);

                Rectangle2D rect = new Rectangle2D.Double(i * probs.dx + probs.minx, j * probs.dy + probs.miny, probs.dx, probs.dy);

                g.setColor(cf.getColor(I1));
                g2d.fill(transform.createTransformedShape(rect));
                g.setColor(Color.black);
                g2d.draw(transform.createTransformedShape(rect));
            }
        }


        System.out.println("Painting shapes");
        for (Shape s : shapes) {
            System.out.print(".");
            System.out.flush();
            GeneralPath path = new GeneralPath(s);
            path.transform(transform);
            g2d.setPaint(Color.BLACK);
            g2d.fill(path);
        }

        g2d.setColor(Color.yellow);
        SpaceTime.paintDot(new SpaceTime(0, start), 4, transform, g2d);

        g2d.setColor(Color.green);
        SpaceTime.paintDot(new SpaceTime(0, end), 4, transform, g2d);

        int ex = probs.x(end.getX());
        int ey = probs.y(end.getY());

        String message = "p=" + probs.p(sx, sy, ex, ey, t) + " r=" + probs.r(sx, sy, ex, ey, t) + " c=" + probs.counts[sx][sy][ex][ey][t];

        g2d.setColor(Color.yellow);
        g2d.drawString(message, 20, getHeight() - 20);
    }

    AffineTransform getFullTransform() {
        AffineTransform transform = getScale();
        transform.concatenate(getTranslate());
        return transform;
    }

    AffineTransform getTranslate() {
        return AffineTransform.getTranslateInstance(-topLeft.getX(), -bottomRight.getY());
    }

    AffineTransform getScale() {
        return AffineTransform.getScaleInstance(scaleX, -scaleY);
    }


    public static void main(String[] args) {

        JFrame frame = new JFrame("Europe");
        frame.getContentPane().add(BorderLayout.CENTER, new VisualizeKMLNumericalProbs(args[0]));
        frame.setSize(900, 900);
        frame.setVisible(true);
    }

}
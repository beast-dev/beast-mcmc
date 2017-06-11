/*
 * VisualizeBrownianBridge2D.java
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

import dr.math.MathUtils;
import dr.math.distributions.MultivariateNormalDistribution;

import javax.swing.*;
import java.awt.*;
import java.awt.geom.AffineTransform;
import java.awt.geom.Ellipse2D;
import java.awt.geom.GeneralPath;
import java.awt.geom.Point2D;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexei Drummond
 * @author Marc Suchard
 */
public class VisualizeBrownianBridge2D extends JComponent {

    MultivariateNormalDistribution mnd;

    // some line segments in space-time
    SpaceTime[] start, end;

    SpaceTimeRejector rejector;

    List<Shape> shapes;

    Paint shapeColor = Color.GRAY;

    Point2D topLeft;
    Point2D bottomRight;

    double scaleX;
    double scaleY;


    public VisualizeBrownianBridge2D() {

        start = new SpaceTime[]{new SpaceTime(0, new double[]{0, 0})};
        end = new SpaceTime[]{new SpaceTime(1, new double[]{1, 1})};

        topLeft = new Point2D.Double(-0.2, -0.2);
        bottomRight = new Point2D.Double(1.2, 1.2);

        shapes = new ArrayList<Shape>();
        shapes.add(new Ellipse2D.Double(0.25, 0.25, 0.4, 0.4));
        shapes.add(new Ellipse2D.Double(0.5, 0.7, 0.2, 0.2));
        shapes.add(new Ellipse2D.Double(0.8, 0.2, 0.15, 0.15));

        rejector = new SpaceTimeRejector() {

            ArrayList<Reject> rejects = new ArrayList<Reject>();

            public boolean reject(double time, double[] space) {
                double x = space[0];
                double y = space[1];
                for (Shape s : shapes) {
                    if (s.contains(x, y)) {
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

        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

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

        for (int r = 0; r < getTrials(); r++) {
            Color c = new Color((float) MathUtils.nextDouble(), (float) MathUtils.nextDouble(), (float)  MathUtils.nextDouble());

            for (int s = 0; s < start.length; s++) {
                List<SpaceTime> points = null;
                g2d.setPaint(c);
                int topLevelRejects = -1;
                while (points == null) {
                    topLevelRejects += 1;
                    points = MultivariateBrownianBridge.divideConquerBrownianBridge(mnd, start[s], end[s], getMaxDepth(), getMaxTries(), rejector);
                }

                GeneralPath path = new GeneralPath();
                path.moveTo((float) points.get(0).getX(0), (float) points.get(0).getX(1));
                //System.out.println(points.get(0));
                for (int i = 1; i < points.size(); i++) {
                    path.lineTo((float) points.get(i).getX(0), (float) points.get(i).getX(1));
                    //System.out.println(points.get(i));
                }

                path.transform(getFullTransform());

                g2d.draw(path);

                g2d.setPaint(Color.black);

                SpaceTime.paintDot(start[s], 3, transform, g2d);
                SpaceTime.paintDot(end[s], 3, transform, g2d);

            }

        }


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

        JFrame frame = new JFrame("Boulders");
        frame.getContentPane().add(BorderLayout.CENTER, new VisualizeBrownianBridge2D());
        frame.setSize(600, 600);
        frame.setVisible(true);
    }

    public int getMaxDepth() {
        return 10;
    }

    public int getTrials() {

        return 10;
    }

    public int getMaxTries() {
        return 10;
    }
}

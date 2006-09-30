/*
 * ContrastPlot.java
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

package dr.app.simcoal;

import dr.evolution.tree.Tree;
import dr.evolution.tree.NodeRef;
import dr.evolution.continuous.Continuous;
import dr.evolution.io.NexusImporter;
import dr.evolution.io.Importer;
import dr.stats.DiscreteStatistics;

import javax.swing.*;
import java.awt.geom.*;
import java.awt.*;
import java.util.Map;
import java.util.StringTokenizer;
import java.util.ArrayList;
import java.util.Collections;
import java.io.IOException;
import java.io.FileReader;
import java.io.BufferedReader;

/**
 */
public class ContrastPlot extends JComponent {

    Tree tree;
    String xLabel, yLabel;
    Rectangle2D viewPort;
    Image backgroundImage;
    Color lineColor = Color.yellow;
    boolean lightMode = true;
    boolean paintBranches = true;
    boolean paintInternalNodes = true;
    boolean showImage = true;
    java.util.List<Point2D> points = null;

    Color pointsColor = new Color(1.0f, 0.0f, 0.0f, 0.4f);
    Color tipColor = Color.black;
    Color internalNodeColor = Color.black;

    private double nodeSize = 4.0;
    private double internalNodeSize = 4.0;
    private double pointSize = 5.0;

    private boolean pointBoundary = false;


    public ContrastPlot(Tree tree,
                        String xLabel,
                        String yLabel,
                        Rectangle2D viewPort,
                        Image backgroundImage,
                        java.util.List<Point2D> points) {
        this.tree = tree;
        this.xLabel = xLabel;
        this.yLabel = yLabel;
        this.viewPort = viewPort;
        this.backgroundImage = backgroundImage;
        this.points = points;

        ImageIcon icon = new ImageIcon(backgroundImage); // ensure image loaded

        setPreferredSize(new Dimension(backgroundImage.getWidth(null),
                backgroundImage.getHeight(null)));
    }


    public void paintComponent(Graphics g) {

        if (showImage) g.drawImage(backgroundImage,0,0,getWidth(), getHeight(), null);

        if (lightMode) {
            // dim image a bit
            g.setColor(new Color(1.0f, 1.0f, 1.0f, 0.2f));
            g.fillRect(0,0,getWidth(), getHeight());
            lineColor = Color.black;
        }

        Graphics2D g2d = (Graphics2D)g;

        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
        g2d.setStroke(new BasicStroke(1.0f));

        if (points != null) for (Point2D point : points) {
            paintPoint(point, g2d, pointsColor, pointSize, true, pointBoundary);
        }


        g2d.setColor(lineColor);
        if (paintBranches) {
            for (int i = 0; i < tree.getNodeCount(); i++) {
                paintBranch(tree, tree.getNode(i), g2d);
            }
        }

        if (paintInternalNodes) {
            for (int i = 0; i < tree.getInternalNodeCount(); i++) {
               paintNode(tree.getInternalNode(i), g2d, internalNodeColor, internalNodeSize);
            }
        }

        for (int i = 0; i < tree.getExternalNodeCount(); i++) {
            paintNode(tree.getExternalNode(i), g2d, tipColor, nodeSize);
        }

    }

    private void paintPoint(Point2D point, Graphics2D g2d, Color fill, double width, boolean transform, boolean boundary) {
        double half = width / 2.0;

        if (transform) {
            //System.out.print(point.getX() + " " + point.getY());
            point = transformOriginal(point);
            //System.out.println(" -> " +
            //        point.getX() + " " + point.getY());
        }


        Ellipse2D circle = new Ellipse2D.Double(
                point.getX() - half,
                point.getY() - half, width, width);

        g2d.setColor(fill);
        g2d.fill(circle);

        if (boundary) {
            Stroke oldStroke = g2d.getStroke();
            g2d.setStroke(new BasicStroke(0.5f));
            g2d.setColor(Color.black);
            g2d.draw(circle);
            g2d.setStroke(oldStroke);
        }

    }

    private void paintNode(NodeRef node, Graphics2D g2d, Color fill, double size) {

            paintPoint(getPoint(node), g2d, fill, size, false, true);
    }

    private void paintBranch(Tree tree, NodeRef node, Graphics2D g2d) {

        if (node == null) throw new IllegalArgumentException("What?!");

        if (tree.isRoot(node)) {

        } else {

            NodeRef parentNode = tree.getParent(node);

            if (parentNode == null) throw new IllegalArgumentException("What?!");

            Point2D parentPoint = getPoint(parentNode);
            Point2D childPoint = getPoint(node);

            //System.out.println(parentPoint.getX() + " " + parentPoint.getY());

            g2d.draw(new Line2D.Double(parentPoint, childPoint));
        }
    }

    private Point2D getPoint(NodeRef node) {

        Point2D original = getPointInOriginalCoordinates(node);
        return transformOriginal(original);

    }

    private Point2D transformOriginal(Point2D original) {
        double xScale = getWidth() / viewPort.getWidth();
        double yScale = getHeight() / viewPort.getHeight();

        return new Point2D.Double(
                (original.getX() - viewPort.getX()) * xScale,
                getHeight() - ((original.getY() - viewPort.getY()) * yScale)
        );
    }

    private static double mercator(double latitude) {

        double phi = latitude * Math.PI / 180.0;

        return Math.log(Math.tan(Math.PI/4.0 + phi/2.0));
    }

    private Point2D getPointInOriginalCoordinates(NodeRef node) {

        double x, y;

        Object xObject = tree.getNodeAttribute(node, xLabel);
        if (xObject instanceof Continuous) x = ((Continuous)xObject).getValue();
        else if (xObject instanceof Double) x = (Double) xObject;

        else if (xObject == null) {
            throw new IllegalArgumentException("xObject is null for node: " + node.getNumber());
        } else {
            throw new IllegalArgumentException("xObject is not a continuous or double for node: " + node.getNumber());
        }

        Object yObject = tree.getNodeAttribute(node, yLabel);
        if (yObject instanceof Continuous) y = ((Continuous)yObject).getValue();
        else if (yObject instanceof Double) y = (Double)yObject;
        else throw new IllegalArgumentException();

        y = mercator(y);

        return new Point2D.Double(x, y);
    }

    private static java.util.List<Point2D> credibleSet(java.util.List<Point2D> points, double proportion) {

        int numToRemove = (int)Math.round(points.size() * (1.0-proportion));

        double[] x = new double[points.size()];
        double[] y = new double[points.size()];
        for (int i = 0; i < x.length; i++) {
            x[i] = points.get(i).getX();
            y[i] = points.get(i).getY();
        }

        double meanX = DiscreteStatistics.mean(x);
        double meanY = DiscreteStatistics.mean(y);
        double stdevX = DiscreteStatistics.stdev(x);
        double stdevY = DiscreteStatistics.stdev(y);

        java.util.List<Double> distances = new ArrayList<Double>();
        for (Point2D point : points) {
            distances.add(normalizedEuclideanDistance(point, meanX, meanY, stdevX, stdevY));
        }

        Collections.sort(distances);

        double cutoffDistance = distances.get(distances.size()-numToRemove);
        System.out.println("cutoffDistance = " + cutoffDistance);

        java.util.List<Point2D> goodPoints = new ArrayList<Point2D>();

        for (Point2D point : points) {
            if (normalizedEuclideanDistance(point, meanX, meanY, stdevX, stdevY) < cutoffDistance) {
                goodPoints.add(point);
            }
        }

        return goodPoints;
    }

    private static double normalizedEuclideanDistance(Point2D point, double meanX, double meanY, double stdevX, double stdevY) {

        double xDistance = Math.abs(point.getX() - meanX) / stdevX;
        double yDistance = Math.abs(point.getY() - meanY) / stdevY;

        return Math.sqrt(xDistance * xDistance + (yDistance * yDistance));

    }

    private static Point2D getMean(java.util.List<Point2D> points) {

        double x = 0.0;
        double y = 0.0;
        for (Point2D point : points) {
            x += point.getX();
            y += point.getY();
        }

        return new Point2D.Double(x / points.size(), y / points.size());
    }

    private static java.util.List<Point2D> readPointsFromFile(String fileName) throws IOException {

        BufferedReader reader = new BufferedReader(new FileReader(fileName));

        String line = reader.readLine();

        ArrayList<Point2D> pointList = new ArrayList<Point2D>();
        while (line != null) {
            StringTokenizer tokens = new StringTokenizer(line);
            double latitude = Double.parseDouble(tokens.nextToken());
            double longitude = Double.parseDouble(tokens.nextToken());

            pointList.add(new Point2D.Double(longitude, mercator(latitude)));

            line = reader.readLine();

        }
        return pointList;
    }

    public static void main(String[] args) throws IOException, Importer.ImportException {

        Map<String, SpaceTime> map = Phylogeography.loadBiekFile(args[0]);
        NexusImporter importer = new NexusImporter(new FileReader(args[1]));

        Tree tree = importer.importNextTree();

        tree = Phylogeography.analyzeContinuousTraits(map, tree);

        double minLatitude = 42.90;
        double rangeInLatitude = 11.17;
        double xMin = -119.15;
        double xRange = 12.25;

        /*double minLatitude = 43.0;
        double rangeInLatitude = 11;
        double xMin = -119.0;
        double xRange = 11.0;*/

        double yMin = mercator(minLatitude);
        double yRange = mercator(minLatitude + rangeInLatitude) - yMin;


        Rectangle2D viewPort = new Rectangle2D.Double(xMin,yMin,xRange,yRange);
        ImageIcon icon = new ImageIcon(args[2]);

        java.util.List<Point2D> points = readPointsFromFile(args[3]);
        System.out.println(points.size() + " points read from file.");

        ContrastPlot plot = new ContrastPlot(tree,
                "longitude", "latitude", viewPort, icon.getImage(), null); //credibleSet(points, 0.95));

        JFrame frame = new JFrame();
        frame.getContentPane().add(plot,BorderLayout.CENTER);
        frame.pack();
        frame.setVisible(true);

    }
}

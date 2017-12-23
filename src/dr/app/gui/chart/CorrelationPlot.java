/*
 * CorrelationPlot.java
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

package dr.app.gui.chart;

import dr.stats.DiscreteStatistics;
import dr.stats.Variate;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;

import java.awt.*;
import java.awt.geom.AffineTransform;
import java.awt.geom.Ellipse2D;
import java.util.*;

/**
 * Description:	A correlation ellipse plot.
 *
 * @author Guy Baele
 */
public class CorrelationPlot extends Plot.AbstractPlot implements Citable {

    private static final boolean PRINT_VISUAL_AIDES = false;

    private final double ELLIPSE_WIDTH = 0.65;
    private final double ELLIPSE_HALF_WIDTH = ELLIPSE_WIDTH / 2;

    private final double NEGATIVE_CORRELATION_DEGREE = 0.785398163;
    private final double POSITIVE_CORRELATION_DEGREE = 2.35619449;

    //colors from plotcorr R package
    private final Color[] colorbrewer_RdBu = {new Color(165,15,21),
            new Color(222, 45, 38),
            new Color(251, 106, 74),
            new Color(252, 174, 145),
            new Color(254, 229, 217),
            new Color( 255, 255, 255),
            new Color(239, 243, 255),
            new Color(189, 215, 231),
            new Color(107, 174, 214),
            new Color(49, 130, 189),
            new Color(8, 81, 156)};

    private final Color[] colorbrewer_RdYlBu = {
            new Color(165,0,38),
            new Color(215,48,39),
            new Color(244,109,67),
            new Color(253,174,97),
            new Color(254,224,144),
            new Color(255,255,191),
            new Color(224,243,248),
            new Color(171,217,233),
            new Color(116,173,209),
            new Color(69,117,180),
            new Color(49,54,149)};

    private final Color[] colorbrewer_Spectral = {
            new Color(158,1,66),
            new Color(213,62,79),
            new Color(244,109,67),
            new Color(253,174,97),
            new Color(254,224,139),
            new Color(255,255,191),
            new Color(230,245,152),
            new Color(171,221,164),
            new Color(102,194,165),
            new Color(50,136,189),
            new Color(94,79,162)};

    private final Color[] translucentColors = {new Color(165,15,21, 32),
            new Color(222, 45, 38, 32),
            new Color(251, 106, 74, 32),
            new Color(252, 174, 145, 32),
            new Color(254, 229, 217, 32),
            new Color(255, 255, 255, 32),
            new Color(239, 243, 255, 32),
            new Color(189, 215, 231, 32),
            new Color(107, 174, 214, 32),
            new Color(49, 130, 189, 32),
            new Color(8, 81, 156, 32)};

    private final Color[] colors = colorbrewer_RdYlBu;
    private final boolean showPoints;

    public CorrelationPlot(String name, java.util.List<Double> xData, java.util.List<Double> yData, boolean showPoints) {
        super(xData, yData);
        setName(name);
        this.showPoints = showPoints;
    }

    /**
     * Paint data series
     */
    protected void paintData(Graphics2D g2, Variate.N xData, Variate.N yData) {

        //System.out.println("CorrelationPlot: paintData");
        //System.out.println("PlotNumber = " + plotNumber);
        //System.out.println("TotalPlotCount = " + this.plotCount);

        int xCount = xData.getCount();
        int yCount = yData.getCount();

        double minX = (Double)xData.getMin();
        double maxX = (Double)xData.getMax();
        double minY = (Double)yData.getMin();
        double maxY = (Double)yData.getMax();

        double xMean = (Double)xData.getMean();
        double yMean = (Double)yData.getMean();

        double[] xDataArray = new double[xCount];
        double[] yDataArray = new double[yCount];
        for (int i = 0; i < xCount; i++) {
            xDataArray[i] = (Double) xData.get(i);
            yDataArray[i] = (Double) yData.get(i);
        }

//        if (this.samples) {
//
//            //Only using a subset of available samples
//            int minCount = Math.min(xCount, yCount);
//            int sampleSize = minCount;
//            /*if (sampleSize < 20) {
//                sampleSize = 20;
//            }*/
//            if (sampleSize > 200) {
//                sampleSize = 200;
//            }
//
//            xDataArray = new double[sampleSize];
//            yDataArray = new double[sampleSize];
//
//            int k = 0;
//            for (int i = 0; i < sampleSize; i++) {
//
//                xDataArray[i] = (Double) xData.get(k);
//                yDataArray[i] = (Double) yData.get(k);
//                k += minCount / sampleSize;
//
//                if (xDataArray[i] < minX) {
//                    minX = xDataArray[i];
//                }
//                if (xDataArray[i] > maxX) {
//                    maxX = xDataArray[i];
//                }
//                if (yDataArray[i] < minY) {
//                    minY = yDataArray[i];
//                }
//                if (yDataArray[i] > maxY) {
//                    maxY = yDataArray[i];
//                }
//            }
//
//            xCount = sampleSize;
//            yCount = sampleSize;
//
//
//        }

        double correlation = DiscreteStatistics.covariance(xDataArray, yDataArray);

        //System.out.println("plotNumber: " + plotNumber + " ; correlation = " + correlation);

        Color fillColor = colors[(int) (5 + 5 * correlation)];
        g2.setColor(fillColor);

        //g2.setPaint(linePaint);
        //g2.setStroke(lineStroke);

        double x1 = (1.0 - ELLIPSE_HALF_WIDTH);
        double y1 = (1.0 - ELLIPSE_HALF_WIDTH);
        double x2 = (1.0 + ELLIPSE_HALF_WIDTH);
        double y2 = (1.0 + ELLIPSE_HALF_WIDTH);

        //System.out.println("(" + x1 + "," + y1 + ") (" + x2 + "," + y2 + ")");

        if (PRINT_VISUAL_AIDES) {
            drawRect(g2, x1, y1, x2, y2);
            g2.drawString("plot" + getPlotNumber(), (float) transformX(x1), (float) transformY(y1));
        }

        double rotationDegree = NEGATIVE_CORRELATION_DEGREE;
        if (correlation > 0) {
            rotationDegree = POSITIVE_CORRELATION_DEGREE;
        }

        double selectedHeight;
        double absCorrelation = Math.abs(correlation);
        selectedHeight = 1.0 - absCorrelation;

        //System.out.println("selectedHeight = " + selectedHeight);

        double ellipseWidth = Math.abs(transformX(x2) - transformX(x1));
        //double ellipseHeight = 0.5 * Math.abs(transformY(y2) - transformY(y1));
        double ellipseHeight = selectedHeight * Math.abs(transformY(y2) - transformY(y1));

        AffineTransform oldTransform = g2.getTransform();

        //double drawOffset = Math.abs(transformY(y2) - transformY(y1)) - ellipseHeight;
        //System.out.println("ellipse draw offset = " + drawOffset);

        Shape ellipse = new Ellipse2D.Double(transformX(x1), transformY(y2), ellipseWidth, ellipseHeight);

        //this transformation rotates around the upper right rectangle corner
        //Shape rotatedEllipse = AffineTransform.getRotateInstance(rotationDegree, transformX(x1), transformY(y2)).createTransformedShape(ellipse);

        //rotate around rectangle center
        //Shape rotatedEllipse = AffineTransform.getRotateInstance(rotationDegree, transformX((x1+x2)/2.0), transformY((y1+y2)/2.0)).createTransformedShape(ellipse);
        //rotate around ellipse center
        Shape rotatedEllipse = AffineTransform.getRotateInstance(rotationDegree, transformX((x1 + x2) / 2.0), transformY(y2) + ellipseHeight / 2.0).createTransformedShape(ellipse);

        //center rotated ellipse within rectangle
        Shape translatedEllipse = AffineTransform.getTranslateInstance(0.0, transformY((y1 + y2) / 2.0) - (transformY(y2) + ellipseHeight / 2.0)).createTransformedShape(rotatedEllipse);

        //g2.fill(rotatedEllipse);
        g2.fill(translatedEllipse);

        g2.setColor(Color.BLACK);
        //g2.draw(rotatedEllipse);

        g2.draw(translatedEllipse);

        if (PRINT_VISUAL_AIDES) {
            //draw ellipse center in black
            g2.fill(new Ellipse2D.Double(transformX((x1 + x2) / 2.0), transformY(y2) + ellipseHeight / 2.0, 5, 5));

            //draw rectangle center in black
            ///g2.fill(new Ellipse2D.Double(transformX((x1+x2)/2.0), transformY((y1+y2)/2.0), 5, 5));
            g2.draw(new Ellipse2D.Double(transformX((x1 + x2) / 2.0), transformY((y1 + y2) / 2.0), 10, 10));
        }

        if (showPoints) {

            g2.setColor(new Color(0, 0, 0, 32));

            for (int i = 0; i < xDataArray.length; i++) {
                //1 unit wide, so divide by maximum data value on both axes
                double newX = x1 + ((xDataArray[i] - minX) * ELLIPSE_WIDTH / (maxX - minX));
                double newY = y1 + ((yDataArray[i] - minY) * ELLIPSE_WIDTH / (maxY - minY));

                //System.out.println("(" + xDataArray[i] + "," + yDataArray[i] + ")  >>>  (" + newX + "," + newY + ")");
                g2.fill(new Ellipse2D.Double(transformX(newX), transformY(newY), 2, 2));
            }

        }

        g2.setTransform(oldTransform);

    }


    @Override
    public Citation.Category getCategory() {
        return Citation.Category.MISC;
    }

    @Override
    public String getDescription() {
        return "method for displaying correlation matrices";
    }

    @Override
    public java.util.List<Citation> getCitations() {
        return Collections.singletonList(CITATION);
    }

    public static Citation CITATION = new Citation(
            new Author[]{
                    new Author("DJ", "Murdoch"),
                    new Author("ED", "Chow")
            },
            "A graphical display of large correlation matrices",
            1986,
            "Am. Stat.",
            50,
            178, 180,
            Citation.Status.PUBLISHED
    );

}

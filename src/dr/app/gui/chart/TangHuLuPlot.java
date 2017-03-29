/*
 * TangHuLuPlot.java
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

package dr.app.gui.chart;

import dr.stats.FrequencyCounter;
import dr.stats.Variate;

import java.awt.*;
import java.awt.geom.Rectangle2D;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

public class TangHuLuPlot extends ScatterPlot {

    protected List<Integer> counts;

    public TangHuLuPlot(List<Double> xData, List<Double> yData) {
        super(xData, yData);
    }


    /**
     * Set data
     */
    public void setData(List<Double> xData, List<Double> yData) {

        if (xData.size() != yData.size())
            throw new IllegalArgumentException("xData and yData must have the same size ! " +
                    xData.size() + " != " + yData.size());

        List<XY> xyList = new ArrayList<XY>();
        for (int i = 0; i < xData.size(); i++) {
            Double x = xData.get(i);
            Double y = yData.get(i);
//            if (x > 0 && y == 0)
//              System.out.println(i + "  :  " + x + "  " + y);
            xyList.add(new XY(x,y));
        }
        // find unique pairs
        FrequencyCounter<XY> xyFC = new FrequencyCounter<XY>(xyList);

        List<Double> xUnique = new ArrayList<Double>();
        List<Double> yUnique = new ArrayList<Double>();
        counts = new ArrayList<Integer>();
        for (XY xy : xyFC.uniqueValues()) {
            xUnique.add(xy.x);
            yUnique.add(xy.y);
            counts.add(xyFC.getCount(xy)); // store counts for circle size
        }

        Variate.D xd = new Variate.D(xUnique);
        Variate.D yd = new Variate.D(yUnique);

        setData(xd, yd);
    }


    /**
     * Draw a mark transforming co-ordinates to each axis
     */
    protected void drawMark(Graphics2D g2, float x, float y, Color color) {

        Rectangle2D bounds = mark.getBounds2D();
        float w = (float) bounds.getWidth();
        float h = (float) bounds.getHeight();
        x = x - (w / 2);
        y = y - (h / 2);

        g2.translate(x, y);

        if (color == null) {
            if (markFillPaint != null) {
                g2.setPaint(markFillPaint);
                g2.fill(mark);
            }
        } else {
            g2.setPaint(color);
            g2.fill(mark);
        }

        g2.setPaint(markPaint);
        g2.setStroke(markStroke);
        g2.draw(mark);

        g2.translate(-x, -y);

        Rectangle2D rect = new Rectangle2D.Float(x, y, w, h);
        markBounds.add(rect);
    }


    /**
     * Draw a mark transforming co-ordinates to each axis
     */
    protected void drawMarkHilighted(Graphics2D g2, float x, float y) {

        Rectangle2D bounds = mark.getBounds2D();
        float w = (float) bounds.getWidth();
        float h = (float) bounds.getHeight();
        x = x - (w / 2);
        y = y - (h / 2);

        g2.translate(x, y);

        if (hilightedMarkFillPaint != null) {
            g2.setPaint(hilightedMarkFillPaint);
            g2.fill(mark);
        }

        g2.setPaint(hilightedMarkPaint);
        g2.setStroke(hilightedMarkStroke);
        g2.draw(mark);

        g2.translate(-x, -y);

        Rectangle2D rect = new Rectangle2D.Float(x, y, w, h);
        markBounds.add(rect);
    }

    /**
     * Paint data series
     */
    protected void paintData(Graphics2D g2, Variate.N xData, Variate.N yData) {
        float x, y;

        markBounds = new java.util.Vector<Rectangle2D>();

        // what is selectedPoints?
        Set<Integer> selectedPoints = getSelectedPoints();

        int n = xData.getCount();
        double circleSize = 20;
        if (n > 1) {
            double xGap = Math.abs(xAxis.getMajorTickSpacing() * xScale);
            double yGap = Math.abs(yAxis.getMajorTickSpacing() * yScale);
            circleSize = Math.min(xGap, yGap) * 0.9;
        }

        for (int i = 0; i < n; i++) {
            x = (float) transformX(((Number) xData.get(i)).doubleValue());
            y = (float) transformY(((Number) yData.get(i)).doubleValue());

            if (selectedPoints.contains(i)) {

                setHilightedMarkStyle(new BasicStroke(1),
                        Color.black, Color.blue);

                drawMarkHilighted(g2, x, y);

            } else {
//counts.get(i)
                // cred set colour
                setMarkStyle(CIRCLE_MARK, circleSize, new BasicStroke(1),
                        Color.black, Color.yellow);

                if (colours != null && colours.size() == n)
                    drawMark(g2, x, y, colours.get(i));
                else
                    drawMark(g2, x, y, null);
            }
        }
    }

}



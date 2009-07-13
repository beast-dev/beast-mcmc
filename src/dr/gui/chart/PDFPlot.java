/*
 * PDFPlot.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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

package dr.gui.chart;

import dr.math.distributions.Distribution;
import dr.stats.Variate;

import java.awt.*;

public class PDFPlot extends Plot.AbstractPlot {

    private Distribution distribution = null;
    private double offset;
    private double xMax, xMin;
    private double yMax;
    private int stepCount = 100;

    /**
     * Constructor
     */
    public PDFPlot(Distribution distribution, double offset) {
        this.distribution = distribution;
        this.offset = offset;
    }

    /**
     * Set data
     */
    public void setData(double[] xData, double[] yData) {
        throw new UnsupportedOperationException("Not available");
    }

    /**
     * Set data
     */
    public void setData(Variate xData, Variate yData) {
        throw new UnsupportedOperationException("Not available");
    }

    /**
     * Set up the axis with some data
     */
    public void setupAxis(Axis xAxis, Axis yAxis, Variate xData, Variate yData) {
        if (distribution == null) {
            return;
        }

        xMin = distribution.quantile(0.005);
        xMax = distribution.quantile(0.995);
        if (Double.isInfinite(xMin)) {
            xMin = 0.0;
        }
        if (Double.isInfinite(xMax)) {
            xMax = 1.0;
        }
        if (xMin == xMax) xMax += 1;

        double x = xMin + offset;
        yMax = distribution.pdf(x - offset);
        double step = (xMax - xMin) / stepCount;
        for (int i = 1; i < stepCount; i++) {
            x += step;
            double y = distribution.pdf(x - offset);
            if (y > yMax) yMax = y;
        }

        if (xAxis instanceof LogAxis) {
            throw new IllegalArgumentException("Log axis are not compatible to PDFPlot");
        } else {
            xAxis.setRange(offset + xMin, offset + xMax);
        }
        if (yAxis instanceof LogAxis) {
            throw new IllegalArgumentException("Log axis are not compatible to PDFPlot");
        } else {
            yAxis.setRange(0.0, yMax);
        }
    }

    /**
     * Paint actual plot
     */
    public void paintPlot(Graphics2D g2, double xScale, double yScale,
                          double xOffset, double yOffset) {

        if (distribution == null) {
            return;
        }

        super.paintPlot(g2, xScale, yScale, xOffset, yOffset);

        g2.setPaint(linePaint);
        g2.setStroke(lineStroke);

        double x1 = xMin + offset;
        double y1 = distribution.pdf(x1 - offset);
        double step = (xMax - xMin) / stepCount;
        for (int i = 1; i < stepCount; i++) {
            double x2 = x1 + step;
            double y2 = distribution.pdf(x2 - offset);
            drawLine(g2, x1, y1, x2, y2);
            x1 = x2;
            y1 = y2;
        }
    }

    /**
     * Paint data series
     */
    protected void paintData(Graphics2D g2, Variate xData, Variate yData) {
        // do nothing because paintPlot is overridden
    }
}

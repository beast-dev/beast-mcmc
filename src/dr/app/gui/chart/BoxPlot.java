/*
 * NumericalDensityPlot.java
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

import dr.stats.Variate;

import java.awt.*;
import java.awt.geom.GeneralPath;

/**
 * @author Andrew Rambaut
 */
public class BoxPlot extends Plot.AbstractPlot { //Plot.AbstractPlot {

    private final boolean isVertical;
    private final double boxWidth;
    private final boolean showTails;
    private final double upper, lower, mean, lowerTail, upperTail;

    protected Stroke meanLineStroke = new BasicStroke(1.5f);
    protected Paint meanLinePaint = Color.black;

    public BoxPlot(boolean isVertical, double boxWidth, double lower, double upper, double lowerTail, double upperTail, double mean) {
        super();

        this.isVertical = isVertical;
        this.boxWidth = boxWidth;
        this.showTails = true;

        this.lower = lower;
        this.upper = upper;
        this.mean = mean;

        this.lowerTail = lowerTail;
        this.upperTail = upperTail;
    }

    /**
     * Set up the axis with some data
     */
    public void setupAxis(Axis xAxis, Axis yAxis, Variate xData, Variate yData) {
        Axis valueAxis = (isVertical ? yAxis : xAxis);
        double minValue = lower;
        double maxValue = upper;
        if (showTails) {
            minValue = lowerTail;
            maxValue = upperTail;
        }

        valueAxis.addRange(minValue, maxValue);
    }

    /**
     * Set line style
     */
    public void setMeanLineStyle(Stroke lineStroke, Paint linePaint) {
        this.meanLineStroke = lineStroke;
        this.meanLinePaint = linePaint;
    }

    protected boolean hasData() {
        return true;
    }

    /**
     * Paint data series
     */
    protected void paintData(Graphics2D g2, Variate.N xData, Variate.N yData) {

        GeneralPath boxPath = new GeneralPath();
        GeneralPath tailPath = new GeneralPath();
        GeneralPath meanPath = new GeneralPath();

        if (isVertical) {
            // rotate -90
            float y1 = (float) transformY(lower);
            float x1 = (float) transformX(getPlotNumber() + 1 - (boxWidth / 2));
            float y2 = (float) transformY(upper);
            float x2 = (float) transformX(getPlotNumber() + 1 + (boxWidth / 2));

            boxPath.moveTo(x1, y1);
            boxPath.lineTo(x2, y1);
            boxPath.lineTo(x2, y2);
            boxPath.lineTo(x1, y2);
            boxPath.lineTo(x1, y1);
            boxPath.closePath();

            if (showTails) {
                float x = (float) transformX(getPlotNumber() + 1);
                float y3 = (float) transformY(lowerTail);
                float y4 = (float) transformY(upperTail);
                tailPath.moveTo(x, y3);
                tailPath.lineTo(x, y1);
                tailPath.moveTo(x, y2);
                tailPath.lineTo(x, y4);
            }

            float y = (float) transformY(mean);
            meanPath.moveTo(x1, y);
            meanPath.lineTo(x2, y);


        } else {
            float x1 = (float) transformX(lower);
            float y1 = (float) transformY(getPlotNumber() + 1 - (boxWidth / 2));
            float x2 = (float) transformX(upper);
            float y2 = (float) transformY(getPlotNumber() + 1 + (boxWidth / 2));

            boxPath.moveTo(x1, y1);
            boxPath.lineTo(x2, y1);
            boxPath.lineTo(x2, y2);
            boxPath.lineTo(x1, y2);
            boxPath.lineTo(x1, y1);
            boxPath.closePath();

            if (showTails) {
                float y = (float) transformY(getPlotNumber() + 1);
                float x3 = (float) transformX(lowerTail);
                float x4 = (float) transformX(upperTail);
                tailPath.moveTo(x3, y);
                tailPath.lineTo(x1, y);
                tailPath.moveTo(x2, y);
                tailPath.lineTo(x4, y);
            }

            float x = (float) transformX(mean);
            meanPath.moveTo(x, y1);
            meanPath.lineTo(x, y2);


        }

        boolean solid = true;
        if (solid) {
            Paint fillPaint = new Color(
                    ((Color) linePaint).getRed(),
                    ((Color) linePaint).getGreen(),
                    ((Color) linePaint).getBlue(), 32);
            g2.setPaint(fillPaint);

            g2.fill(boxPath);
        }

        g2.setPaint(meanLinePaint);
        g2.setStroke(meanLineStroke);
        g2.draw(meanPath);

        g2.setStroke(lineStroke);
        g2.setPaint(linePaint);

        if (showTails) {
            g2.draw(tailPath);
        }

        g2.draw(boxPath);

    }


}

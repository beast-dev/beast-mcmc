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

import dr.stats.FrequencyCounter;
import dr.stats.Variate;

import java.awt.*;

/**
 * @author Andrew Rambaut
 */
public class IntegerViolinPlot extends Plot.AbstractPlot { //Plot.AbstractPlot {

    private Paint barPaint = new Color(124, 164, 221);
    private Paint intervalTailPaint = new Color(232, 114, 103);

    private boolean hasIntervals = false;
    private double upperInterval = 0.0;
    private double lowerInterval = 0.0;
    private final boolean isVertical;
    private final boolean showTails;

    private double violinWidth = 0.9;

    public IntegerViolinPlot(boolean isVertical, double violinWidth, double lower, double upper, boolean showTails, FrequencyCounter<Integer> frequencyCounter) {
        super();
        this.isVertical = isVertical;
        this.violinWidth = violinWidth;
        this.showTails = showTails;

        this.lowerInterval = lower;
        this.upperInterval = upper;

        setData(frequencyCounter);
    }

    private void setData(FrequencyCounter<Integer> frequencyCounter) {
        double scale = (0.5 * violinWidth) / frequencyCounter.getMaxFrequency();
        Variate.D xd = new Variate.D();
        Variate.D yd = new Variate.D();

        java.util.List<Integer> values = frequencyCounter.getUniqueValues();

        int i = 0;
        for (int value : values) {
            xd.add((double) value);
            yd.add((double)frequencyCounter.getFrequency(value) * scale);
        }

        setData(xd, yd);
    }

    /**
     * Set up the axis with some data
     */
    public void setupAxis(Axis xAxis, Axis yAxis, Variate xData, Variate yData) {
        Axis valueAxis = (isVertical ? yAxis : xAxis);
        Axis densityAxis = (isVertical ? xAxis : yAxis);
        Variate valueData = xData;
        Variate densityData = yData;

        if (valueData != null) {
            if (valueAxis instanceof LogAxis) {
                double minValue = java.lang.Double.POSITIVE_INFINITY;

                for (int i = 0; i < valueData.getCount(); i++) {
                    double value = (Double) valueData.get(i);
                    if (value > 0.0 && value < minValue)
                        minValue = value;
                }

                valueAxis.addRange(minValue, (Double) valueData.getMax());
            } else {
                valueAxis.addRange((Double) valueData.getMin(), (Double) valueData.getMax());
            }
        }

        if (densityData != null) {
            if (densityAxis instanceof LogAxis) {
                throw new IllegalArgumentException("The density axis shouldn't be log scale");
            } else {
                // do nothing - the range will be set up by the chart
            }
        }
    }

    /**
     * Set paints
     */
    public void setPaints(Paint barPaint, Paint intervalTailPaint) {
        this.barPaint = barPaint;
        this.intervalTailPaint = intervalTailPaint;
    }

    /**
     * Set arbitrary intervals to use (0 for none).
     */
    public void setIntervals(double lower, double upper) {
        hasIntervals = true;
        upperInterval = upper;
        lowerInterval = lower;
    }

    /**
     * Paint data series
     */
    protected void paintData(Graphics2D g2, Variate.N xData, Variate.N yData) {

        double x, x1, y1, x2, y2;
        int n = xData.getCount();

        g2.setStroke(lineStroke);

        for (int i = 0; i < n; i ++) {

            x = (Double)xData.get(i);

            x1 = x - 0.5;
            y1 = getPlotNumber() + 1 - (Double) yData.get(i);
            x2 = x + 0.5;
            y2 = getPlotNumber() + 1 + (Double) yData.get(i);

            if (barPaint != null) {
                if (hasIntervals && (x < lowerInterval || x > upperInterval)) {
                    g2.setPaint(intervalTailPaint);
                } else {
                    g2.setPaint(barPaint);
                }
                if (isVertical) {
                    fillRect(g2, y1, x1, y2, x2);
                } else {
                    fillRect(g2, x1, y1, x2, y2);
                }
            }

            if (lineStroke != null && linePaint != null) {
                g2.setStroke(lineStroke);
                g2.setPaint(linePaint);
                if (isVertical) {
                    drawRect(g2, y1, x1, y2, x2);
                } else {
                    drawRect(g2, x1, y1, x2, y2);
                }
            }
        }
    }

}

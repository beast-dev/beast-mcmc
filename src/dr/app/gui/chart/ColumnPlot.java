/*
 * FrequencyPlot.java
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

import dr.inference.trace.TraceDistribution;
import dr.inference.trace.TraceType;
import dr.stats.FrequencyCounter;
import dr.stats.Variate;
import dr.util.FrequencyDistribution;

import java.awt.*;
import java.util.List;
import java.util.Map;

public class ColumnPlot extends Plot.AbstractPlot {

    //    protected Paint barPaint = Color.blue;
    private Paint barPaint = new Color(124, 164, 221);
    private Paint intervalTailPaint = new Color(232, 114, 103);

    private boolean hasIntervals = false;
    private double upperInterval = 0.0;
    private double lowerInterval = 0.0;

    private double columnWidth = 0.9;
    private boolean columnsOnTicks = true;

    public ColumnPlot(FrequencyCounter<Integer> frequencyCounter, List<Integer> order, boolean showFrequency) {
        super();
        setData(frequencyCounter, order, showFrequency);
    }

    public ColumnPlot(List<Double> xData, List<Double> yData) {
        super();
        setData(xData, yData);
    }

    public ColumnPlot(Variate.D xData, Variate.D yData) {
        super();
        setData(xData, yData);
    }

    private void setData(FrequencyCounter<Integer> frequencyCounter, List<Integer> order, boolean showFrequency) {
        Variate.D xd = new Variate.D();
        Variate.D yd = new Variate.D();

        int i = 0;
        if (order == null) {
            order = frequencyCounter.getUniqueValues();
        }
        for (int value : order) {
            xd.add((double)i);
            yd.add(showFrequency ? frequencyCounter.getFrequency(value) : frequencyCounter.getProbability(value));
            i++;
        }

        setData(xd, yd);
    }

    public void setColumnWidth(double columnWidth) {
        this.columnWidth = columnWidth;
    }

    /**
     * Set paints
     */
    public void setPaints(Paint barPaint, Paint quantilePaint) {
        this.barPaint = barPaint;
        this.intervalTailPaint = quantilePaint;
    }

    /**
     * Set arbitrary intervals to use (0 for none).
     */
    public void setIntervals(double lower, double upper) {
        hasIntervals = (upper > 0.0 || lower > 0.0);
        upperInterval = upper;
        lowerInterval = lower;
    }

    /**
     * Set bar fill style. Use a barPaint of null to not fill bar.
     * Bar outline style is set using setLineStyle
     */
    public void setBarFillStyle(Paint barPaint) {
        this.barPaint = barPaint;
    }

    /**
     * Paint data series
     */
    protected void paintData(Graphics2D g2, Variate.N xData, Variate.N yData) {

        double x, x1, y1, x2, y2;
        int n = xData.getCount();

        double delta = (1.0 - columnWidth) / 2;

        g2.setStroke(lineStroke);
        for (int i = 0; i < n; i ++) {

            x = (Double)xData.get(i);

            x1 = x + delta;
            y1 = 0.0;
            x2 = x + 1 - delta;
            y2 = (Double) yData.get(i);

            if (columnsOnTicks) {
                x1 -= 0.5;
                x2 -= 0.5;
            }

            if (y1 != y2) {
                if (barPaint != null) {
                    if (hasIntervals && (x < lowerInterval || x > upperInterval)) {
                        g2.setPaint(intervalTailPaint);
                    } else {
                        g2.setPaint(barPaint);
                    }
                    fillRect(g2, x1, y1, x2, y2);
                }

                if (lineStroke != null && linePaint != null) {
                    g2.setStroke(lineStroke);
                    g2.setPaint(linePaint);
                    drawRect(g2, x1, y1, x2, y2);
                }
            }
        }
    }

}

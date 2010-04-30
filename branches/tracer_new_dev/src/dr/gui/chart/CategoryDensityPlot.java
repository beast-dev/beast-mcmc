/*
 * CategoryFrequencyPlot.java
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

import dr.inference.trace.TraceDistribution;
import dr.stats.Variate;
import dr.util.FrequencyDistribution;

import java.awt.*;

public class CategoryDensityPlot extends FrequencyPlot {


    public CategoryDensityPlot(String[] data, int minimumBinCount, TraceDistribution traceD) {
        super(traceD);

        categoryData = new String[data.length];
        double[] doubleData = new double[data.length];
        for (int i = 0; i < data.length; i++) {
            doubleData[i] = (double) traceD.credSet.getIndex(data[i]);
            categoryData[i] = data[i];
        }
        setData(doubleData, minimumBinCount);
    }

    public CategoryDensityPlot(int[] data, int minimumBinCount, TraceDistribution traceD) {
        super(traceD);
        double[] doubleData = new double[data.length];
        for (int i = 0; i < data.length; i++) {
            doubleData[i] = (double) data[i];
        }
        setData(doubleData, minimumBinCount);
    }

    /**
     * Set data
     */
    public void setData(Variate data, int minimumBinCount) {
        raw = data;
        FrequencyDistribution frequency = getFrequencyDistribution(data, minimumBinCount);

        Variate.Double xData = new Variate.Double();
        Variate.Double yData = new Variate.Double();

        double x = frequency.getLowerBound();

        for (int i = 0; i < frequency.getBinCount(); i++) {

            xData.add(x);
            yData.add(0.0);

            x += frequency.getBinSize();

            xData.add(x);
            yData.add(frequency.getProb(i));

        }
        setData(xData, yData);
    }


    /**
     * Get the FrequencyDistribution object
     */
    protected FrequencyDistribution getFrequencyDistribution(Variate data, int minimumBinCount) {
        double min = data.getMin();
        double max = data.getMax();

        if (min == max) {
            if (min == 0) {
                min = -1.0;
            } else {
                min -= Math.abs(min / 10.0);
            }
            if (max == 0) {
                max = 1.0;
            } else {
                max += Math.abs(max / 10.0);
            }
        }

        Axis axis = new LinearAxis(Axis.AT_MAJOR_TICK, Axis.AT_MAJOR_TICK);
        axis.setRange(min, max);

        int majorTickCount = axis.getMajorTickCount();
        axis.setPrefNumTicks(majorTickCount, 4);

        double binSize = axis.getMinorTickSpacing();
        int binCount = (int) ((axis.getMaxAxis() - axis.getMinAxis()) / binSize) + 2;

        while (binCount < minimumBinCount) {
            majorTickCount++;
            axis.setPrefNumTicks(majorTickCount, 4);

            binSize = axis.getMinorTickSpacing();
            binCount = (int) ((axis.getMaxAxis() - axis.getMinAxis()) / binSize) + 2; // should +2, otherwise the last bar will lose
        }

        FrequencyDistribution frequency = new FrequencyDistribution(axis.getMinAxis(), binCount, binSize);

        for (int i = 0; i < raw.getCount(); i++) {
            frequency.addValue(raw.get(i));
        }

        return frequency;
    }

    /**
     * Paint data series
     */
    protected void paintData(Graphics2D g2, Variate xData, Variate yData) {

        double x1, y1, x2, y2;

        int n = xData.getCount();


        g2.setStroke(lineStroke);
        for (int i = 0; i < n; i += 2) {

            x1 = xData.get(i);
            y1 = yData.get(i);
            x2 = xData.get(i + 1);
            y2 = yData.get(i + 1);

            if (y1 != y2) {
                if (barPaint != null) {
                        g2.setPaint(barPaint);
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



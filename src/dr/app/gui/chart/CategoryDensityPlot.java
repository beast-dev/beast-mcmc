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

package dr.app.gui.chart;

import dr.inference.trace.TraceDistribution;
import dr.stats.Variate;
import dr.util.FrequencyDistribution;

import java.awt.*;
import java.util.List;

public class CategoryDensityPlot extends FrequencyPlot {
    private int barCount = 0;
    private int barId;
    // for string[], passing the int[] storing the index of string[]

    public CategoryDensityPlot(List<Double> data, int minimumBinCount, TraceDistribution traceDistribution, int barCount, int barId) {
        super(traceDistribution);
        this.barCount = barCount;
        this.barId = barId;

//        Double[] doubleData = new Double[data.length];
//        for (int i = 0; i < data.length; i++) {
//            doubleData[i] = data[i].doubleValue();
//        }
        setData(data, minimumBinCount);
    }

    /**
     * Set data
     */
    public void setData(Variate data, int minimumBinCount) {
        raw = data;
        FrequencyDistribution frequency = getFrequencyDistribution(data, minimumBinCount);

        Variate.D xData = new Variate.D();
        Variate.D yData = new Variate.D();

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
        double min = (Double) data.getMin();
        double max = (Double) data.getMax();

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

        if (minimumBinCount > 0) {
            while (binCount < minimumBinCount) {
                majorTickCount++;
                axis.setPrefNumTicks(majorTickCount, 4);

                binSize = axis.getMinorTickSpacing();
                binCount = (int) ((axis.getMaxAxis() - axis.getMinAxis()) / binSize) + 2; // should +2, otherwise the last bar will lose
            }
        }

        FrequencyDistribution frequency = new FrequencyDistribution(axis.getMinAxis(), binCount, binSize);

        for (int i = 0; i < raw.getCount(); i++) {
            frequency.addValue((Double) raw.get(i));
        }

        return frequency;
    }

    /**
     * Paint data series
     */
    protected void paintData(Graphics2D g2, Variate.D xData, Variate.D yData) {
        double x1, y1, x2, y2, x;

        int n = xData.getCount();

        g2.setStroke(lineStroke);
        for (int i = 0; i < n; i += 2) {


               x1 = xData.get(i);
               x2 = xData.get(i + 1);
               x = x2 - x1;

            if (barCount > 1) {
                x1 = x1 - ((double) (barCount - 1)) * x + 2.0 * ((double) barId) * x;
                x2 = x2 - ((double) (barCount - 1)) * x + 2.0 * ((double) barId) * x;
            }

            y1 = yData.get(i);
            y2 = yData.get(i + 1);


            if (y1 != y2) {
                if (linePaint != null) {
                    Paint fillPaint = new Color(
                    ((Color) linePaint).getRed(),
                    ((Color) linePaint).getGreen(),
                    ((Color) linePaint).getBlue(), 125);
                    g2.setPaint(fillPaint);
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

//    protected void fillRect(Graphics2D g2, double x1, double y1, double x2, double y2) {
//        if (traceD != null && traceD.getTraceType() != TraceFactory.TraceType.DOUBLE && numOfBars > 0) {
//            super.fillRect(g2, x1 - ((double) (numOfBars - barId)) * 2.0 * (x2-x1), y1, x2 - ((double) (numOfBars - barId)) * 2.0 * (x2-x1), y2);
//        } else {
//            super.fillRect(g2, x1, y1, x2, y2);
//        }
//    }
//
//    protected void drawRect(Graphics2D g2, double x1, double y1, double x2, double y2) {
//        if (traceD != null && traceD.getTraceType() != TraceFactory.TraceType.DOUBLE && numOfBars > 0) {
////            super.drawRect(g2, x1-(x2-x1), y1, x2, y2);
//           super.fillRect(g2, x1 - ((double) (numOfBars - barId)) * 2.0 * (x2-x1), y1, x2 - ((double) (numOfBars - barId)) * 2.0 * (x2-x1), y2);
//        } else {
//            super.drawRect(g2, x1, y1, x2, y2);
//        }
//    }
}



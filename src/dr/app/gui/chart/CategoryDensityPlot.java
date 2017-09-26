/*
 * CategoryDensityPlot.java
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
import dr.stats.Variate;
import dr.util.FrequencyDistribution;

import java.awt.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

public class CategoryDensityPlot extends FrequencyPlot {
    private final static double BAR_WIDTH = 0.8;
    private int barId;

    public CategoryDensityPlot(List<Double> data, TraceDistribution traceDistribution, int barId) {
        super(traceDistribution);
        this.barId = barId;

        if (!traceDistribution.getTraceType().isDiscrete())
            throw new IllegalArgumentException("Discrete trace is required for frequency plot.");

        setData(new Variate.D(data), traceDistribution.getCategoryOrderMap());
    }

    /**
     * Set data, all integers
     */
    public void setData(Variate.D data, Map<Integer, Integer> categoryOrderMap) {
        setRawData(data);
        FrequencyDistribution frequency = getDiscreteFrequencyDistribution(data, categoryOrderMap);

        Variate.D xData = new Variate.D();
        Variate.D yData = new Variate.D();

        Map<Integer, Integer> orderMap = traceDistribution.getCategoryOrderMap();
        double x = frequency.getLowerBound();

        for (int i = 0; i < frequency.getBinCount(); i++) {
            int index = i;

            if (orderMap != null && orderMap.size() > 0) {
                if (i >= orderMap.size()) {
                    System.out.println("oops");
                }
                index = orderMap.get(i);
            }
            
            xData.add((double)x);
            yData.add(0.0);

            x += frequency.getBinSize();

            xData.add((double)x);
            yData.add(frequency.getProbability(index));

        }
        setData(xData, yData);
    }

    /**
     * Paint data series
     */
    protected void paintData(Graphics2D g2, Variate.N xData, Variate.N yData) {
        int n = xData.getCount();

        for (int i = 0; i < n; i += 2) {
            double cellWidth = Math.abs(xAxis.getMajorTickSpacing() * xScale) * BAR_WIDTH;

            float x = (float) transformX(((Number) xData.get(i)).doubleValue());
            float x1 = x - (float)cellWidth / 2;
            float x2 = x + (float)cellWidth / 2;

            float y1 = (float) transformY(((Number) yData.get(i)).doubleValue());
            float y2 = (float) transformY(((Number) yData.get(i + 1)).doubleValue());


//            x = x2 - x1;

//            if (barCount > 1) {
//                x1 = x1 - ((double) (barCount - 1)) * x + 2.0 * ((double) barId) * x;
//                x2 = x2 - ((double) (barCount - 1)) * x + 2.0 * ((double) barId) * x;
//            }

//            y1 = (Double) yData.get(i);
//            y2 = (Double) yData.get(i + 1);


            if (y1 != y2) {
                if (linePaint != null) {
                    Paint fillPaint = new Color(
                            ((Color) linePaint).getRed(),
                            ((Color) linePaint).getGreen(),
                            ((Color) linePaint).getBlue(), 125);
                    if (barId==0)
                        fillPaint = new Color(124, 164, 221, 125);
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
}



///*
// * TangHuLuPlot.java
// *
// * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
// *
// * This file is part of BEAST.
// * See the NOTICE file distributed with this work for additional
// * information regarding copyright ownership and licensing.
// *
// * BEAST is free software; you can redistribute it and/or modify
// * it under the terms of the GNU Lesser General Public License as
// * published by the Free Software Foundation; either version 2
// * of the License, or (at your option) any later version.
// *
// *  BEAST is distributed in the hope that it will be useful,
// *  but WITHOUT ANY WARRANTY; without even the implied warranty of
// *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// *  GNU Lesser General Public License for more details.
// *
// * You should have received a copy of the GNU Lesser General Public
// * License along with BEAST; if not, write to the
// * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
// * Boston, MA  02110-1301  USA
// */
//
//package dr.app.gui.chart;
//
//import dr.inference.trace.TraceDistribution;
//import dr.stats.Variate;
//import dr.util.FrequencyDistribution;
//
//import java.awt.*;
//import java.util.HashMap;
//import java.util.List;
//
//public class TangHuLuPlot extends FrequencyPlot {
//     protected final HashMap<Integer, String> categoryDataMap;
//
//    // for string[], passing the int[] storing the index of string[]
//    public TangHuLuPlot(List<Double> data, TraceDistribution traceDistribution, HashMap<Integer, String> categoryDataMap) {
//        super(traceDistribution);
//        if (categoryDataMap == null)
//            throw new IllegalArgumentException("Categorical data map is null !");
//        this.categoryDataMap = categoryDataMap;
//
//        setData(new Variate.D(data));
//    }
//
//    /**
//     * Set data, all integers
//     */
//    public void setData(Variate.D data) {
//        setRawData(data);
//        FrequencyDistribution frequency = getFrequencyDistribution(data);
//
//        Variate.D xData = new Variate.D();
//        Variate.D yData = new Variate.D();
//
//        double x = frequency.getLowerBound();
//
//        for (int i = 0; i < frequency.getBinCount(); i++) {
//
//            xData.add(x);
//            yData.add(0.0);
//
//            x += frequency.getBinSize();
//
//            xData.add(x);
//            yData.add(frequency.getProbability(i));
//
//        }
//        setData(xData, yData);
//    }
//
//    protected FrequencyDistribution getFrequencyDistribution(Variate data) {
//        double min = (Double) data.getMin();
//        double max = (Double) data.getMax();
//
//        if (min == max) {
//            if (min == 0) {
//                min = -1.0;
//            } else {
//                min -= 1;
//            }
//            if (max == 0) {
//                max = 1.0;
//            } else {
//                max += 1;
//            }
//        }
//
//        Axis axis = new LinearAxis(Axis.AT_MAJOR_TICK_PLUS, Axis.AT_MAJOR_TICK_PLUS);
//        axis.setRange(min, max);
//
//        int majorTickCount = axis.getMajorTickCount();
//        axis.setPrefNumTicks(majorTickCount, 4);
//
//        double binSize = axis.getMinorTickSpacing();
//        int binCount = (int) ((axis.getMaxAxis() - axis.getMinAxis()) / binSize) + 2;
//
//        double start = Math.floor(axis.getMinAxis()); // to convert x-axis into integer
//        FrequencyDistribution frequency = new FrequencyDistribution(start, binCount, binSize);
//
//        for (int i = 0; i < raw.getCount(); i++) {
//            frequency.addValue((Double) raw.get(i));
//        }
//
//        return frequency;
//    }
//
//    /**
//     * Paint data series
//     */
//    protected void paintData(Graphics2D g2, Variate.N xData, Variate.N yData) {
//        double x1, y1, x2, y2, x;
//
//        int n = xData.getCount();
//
//        for (int i = 0; i < n; i += 2) {
//            x1 = (Double) xData.get(i);
//            x2 = (Double) xData.get(i + 1);
//            x = x2 - x1;
//
//            if (barCount > 1) {
//                x1 = x1 - ((double) (barCount - 1)) * x + 2.0 * ((double) barId) * x;
//                x2 = x2 - ((double) (barCount - 1)) * x + 2.0 * ((double) barId) * x;
//            }
//
//            y1 = (Double) yData.get(i);
//            y2 = (Double) yData.get(i + 1);
//
//
//            if (y1 != y2) {
//                if (linePaint != null) {
//                    Paint fillPaint = new Color(
//                            ((Color) linePaint).getRed(),
//                            ((Color) linePaint).getGreen(),
//                            ((Color) linePaint).getBlue(), 125);
//                    if (barId==0)
//                        fillPaint = new Color(124, 164, 221, 125);
//                    g2.setPaint(fillPaint);
//                    fillRect(g2, x1, y1, x2, y2);
//                }
//
//                if (lineStroke != null && linePaint != null) {
//                    g2.setStroke(lineStroke);
//                    g2.setPaint(linePaint);
//                    drawRect(g2, x1, y1, x2, y2);
//                }
//            }
//        }
//    }
//
//}
//
//

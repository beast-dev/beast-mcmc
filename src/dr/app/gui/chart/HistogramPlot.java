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
import dr.stats.Variate;
import dr.util.FrequencyDistribution;

import java.awt.*;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class HistogramPlot extends Plot.AbstractPlot {

    protected Variate raw = null;

    //    protected Paint barPaint = Color.blue;
    private Paint barPaint = new Color(124, 164, 221);
    private Paint quantilePaint = new Color(232, 114, 103);

    private boolean hasQuantiles = false;
    private double quantiles = 0;

    private boolean hasIntervals = false;
    private double upperInterval = 0.0;
    private double lowerInterval = 0.0;

    private Set<Integer> incredibleSet = Collections.EMPTY_SET;

//    private boolean hasIncredibleSet = false;
//    private TraceDistribution.CredibleSetAnalysis credSet;


    public HistogramPlot(List<Double> data, int minimumBinCount) {
        super();
        setData(data, minimumBinCount);
    }

    public HistogramPlot(Variate.D data, int minimumBinCount) {
        super();
        setData(data, minimumBinCount);
    }

    public HistogramPlot(List<Double> data, int minimumBinCount, TraceDistribution traceDistribution) {
        super();
        setData(data, minimumBinCount);
    }

//    public FrequencyPlot(List<Integer> data, TraceType traceType) {
//        this(traceType);
//        if (!traceType.isCategorical()) {
//            throw new IllegalArgumentException("Categorical value is required for frequency plot.");
//        }
//
////        List<Double> intData = traceDistribution.indexingData(data);
//        // set data by index of unique categorical values
//        setIntegerData(data, 1);
//    }

    /**
     * Set data
     */
    public void setData(List<Double> data, int minimumBinCount) {
        Variate.D d = new Variate.D(data);
        setData(d, minimumBinCount);
    }

    /**
     * Set paints
     */
    public void setPaints(Paint barPaint, Paint quantilePaint) {
        this.barPaint = barPaint;
        this.quantilePaint = quantilePaint;
    }

    /**
     * Set data
     */
    public void setData(Variate.D data, int minimumBinCount) {

        setRawData(data);
        FrequencyDistribution frequency = getFrequencyDistribution(data, minimumBinCount);

        Variate.D xData = new Variate.D();
        Variate.D yData = new Variate.D();

        double x = frequency.getLowerBound();

        for (int i = 0; i < frequency.getBinCount(); i++) {

            xData.add(x);
            yData.add(0.0);

            x += frequency.getBinSize();

            xData.add(x);
            yData.add((double) frequency.getFrequency(i));

        }
        setData(xData, yData);
    }

    /**
     * Set data
     */
//    public void setIntegerData(Variate.I data, int minimumBinCount) {
//
//        setRawData(data);
//        FrequencyDistribution frequency = getFrequencyDistribution(data, minimumBinCount);
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
//            yData.add((double) frequency.getFrequency(i));
//
//        }
//        setData(xData, yData);
//    }

    protected Variate getRawData() {
        return raw;
    }

    protected void setRawData(Variate raw) {
        this.raw = raw;
    }

    /**
     * Get the FrequencyDistribution object
     */
    protected FrequencyDistribution getFrequencyDistribution(Variate data, int minimumBinCount) {
        double min = ((Number) data.getMin()).doubleValue();
        double max = ((Number) data.getMax()).doubleValue();

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
        if (minimumBinCount <= 0) {
            axis = new LinearAxis(Axis.AT_MAJOR_TICK_PLUS, Axis.AT_MAJOR_TICK_PLUS);
        }
        axis.setRange(min, max);

        int majorTickCount = axis.getMajorTickCount();
        axis.setPrefNumTicks(majorTickCount, 4);

        double binSize = axis.getMinorTickSpacing();
        int binCount = (int) ((axis.getMaxAxis() - axis.getMinAxis()) / binSize) + 2;

        if (minimumBinCount > 0) {
            // avoid dead loop
            while (binCount < minimumBinCount && majorTickCount < 1000) {
                majorTickCount++;
                axis.setPrefNumTicks(majorTickCount, 4);

                binSize = axis.getMinorTickSpacing();
                // should +2, otherwise the last bar will lose
                binCount = (int) ((axis.getMaxAxis() - axis.getMinAxis()) / binSize) + 2;
            }
        } else if (binSize > 1) {
            // getMinorTickSpacing() returns 1.25, if the min integer slightly bigger than 0
            binSize = 0.5;
            binCount = (int) ((axis.getMaxAxis() - axis.getMinAxis()) / binSize) + 2;
        }

        double start = axis.getMinAxis();
        if (minimumBinCount < 0)
            start = Math.round(start); // to convert x-axis into integer
        FrequencyDistribution frequency = new FrequencyDistribution(start, binCount, binSize);

        for (int i = 0; i < raw.getCount(); i++) {
            frequency.addValue(((Number) raw.get(i)).doubleValue());
        }

        return frequency;
    }

    /**
     * Get the FrequencyDistribution object
     */
//    protected FrequencyDistribution getDiscreteFrequencyDistribution(Variate data, Map<Integer, Integer> categoryOrderMap) {
//        int min = ((Number) data.getMin()).intValue();
//        int max = ((Number) data.getMax()).intValue();
//
//        if (min == max) {
//            min -= -1.0;
//            max += 1.0;
//        }
//
//        double binSize = 1;
//        int binCount = max - min + 1;
//
//        double start = min;
//
//        FrequencyDistribution frequency = new FrequencyDistribution(start, binCount, binSize);
//
//        for (int i = 0; i < raw.getCount(); i++) {
//            int index = categoryOrderMap.get(((Number)raw.get(i)).intValue());
//            frequency.addValue(index);
//        }
//
//        return frequency;
//    }

    /**
     * Set quantiles to use (0 for none).
     */
    public void setQuantiles(double quantiles) {
        this.quantiles = quantiles;
        hasQuantiles = (quantiles > 0.0);
        hasIntervals = false;
    }

    /**
     * Set arbitrary intervals to use (0 for none).
     */
    public void setIntervals(double upper, double lower) {
        hasQuantiles = false;
        hasIntervals = (upper > 0.0 || lower > 0.0);
        upperInterval = upper;
        lowerInterval = lower;
    }

    /**
     * Use incredible set for tails.
     */
    public void setIncredibleSet(Set<Integer> incredibleSet) {
        this.incredibleSet = incredibleSet;
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

        double x1, y1, x2, y2;
        double lower = 0.0, upper = 0.0;
        int n = xData.getCount();

        if (hasQuantiles) {
            lower = getRawData().getQuantile(quantiles);
            upper = getRawData().getQuantile(1.0 - quantiles);
        } else if (hasIntervals) {
            lower = lowerInterval;
            upper = upperInterval;
        }

        g2.setStroke(lineStroke);
        for (int i = 0; i < n; i += 2) {

            x1 = (Double) xData.get(i);
            y1 = (Double) yData.get(i);
            x2 = (Double) xData.get(i + 1);
            y2 = (Double) yData.get(i + 1);

            if (y1 != y2) {
                if (barPaint != null) {
                    if (hasQuantiles || hasIntervals) {
                        if (x1 < lower) {
                            if (x2 <= lower) {
                                g2.setPaint(quantilePaint);
                                fillRect(g2, x1, y1, x2, y2);
                            } else {
                                g2.setPaint(quantilePaint);
                                fillRect(g2, x1, y1, lower, y2);
                                g2.setPaint(barPaint);
                                fillRect(g2, lower, y1, x2, y2);
                            }
                        } else if (x2 > upper) {
                            if (x1 >= upper) {
                                g2.setPaint(quantilePaint);
                                fillRect(g2, x1, y1, x2, y2);
                            } else {
                                g2.setPaint(barPaint);
                                fillRect(g2, x1, y1, upper, y2);
                                g2.setPaint(quantilePaint);
                                fillRect(g2, upper, y1, x2, y2);
                            }
                        } else {
                            g2.setPaint(barPaint);
                            fillRect(g2, x1, y1, x2, y2);
                        }
                    } else {
                        g2.setPaint(barPaint);
                        fillRect(g2, x1, y1, x2, y2);
                    }
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

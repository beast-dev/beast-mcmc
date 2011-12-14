/*
 * FrequencyPlot.java
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
import dr.inference.trace.TraceFactory;
import dr.stats.Variate;
import dr.util.FrequencyDistribution;

import java.awt.*;
import java.util.List;

public class FrequencyPlot extends Plot.AbstractPlot {

    protected Variate raw = null;

//    protected Paint barPaint = Color.blue;
    protected Paint barPaint = new Color(124, 164, 221);
    private Paint quantilePaint = new Color(232, 114, 103);

    private boolean hasQuantiles = false;
    private double quantiles = 0;

    private boolean hasIntervals = false;
    private double upperInterval = 0.0;
    private double lowerInterval = 0.0;

    private boolean hasIncredibleSet = false;
//    private TraceDistribution.CredibleSet credSet;

    protected TraceDistribution traceDistribution = null;

    protected FrequencyPlot(TraceDistribution traceDistribution) {
        super();
        this.traceDistribution = traceDistribution;
    }

    public FrequencyPlot(Variate.D data, int minimumBinCount) {
        super();
        setData(data, minimumBinCount);
    }

    public FrequencyPlot(List<Double> data, int minimumBinCount) {
        super();
        setData(data, minimumBinCount);
    }

    public FrequencyPlot(List<Double> data, int minimumBinCount, TraceDistribution traceDistribution) {
        this(traceDistribution);
        setData(data, minimumBinCount);
    }

//    public FrequencyPlot(Integer[] data, int minimumBinCount, TraceDistribution traceD) {
//        this(traceD);
//        Double[] doubleData = new Double[data.length];
//        for (int i = 0; i < data.length; i++) {
//            doubleData[i] = data[i].doubleValue();
//        }
//        setData(doubleData, minimumBinCount);
//    }

//    public FrequencyPlot(String[] data, int minimumBinCount, TraceDistribution traceD) {
//        this(traceD);
//        categoryDataMap.clear();
//        double[] doubleData = new double[data.length];
//        for (int i = 0; i < data.length; i++) {
//            doubleData[i] = (double) traceD.getIndex(data[i]);
//            categoryDataMap.put(doubleData[i], data[i]);
//        }
//        setData(doubleData, minimumBinCount);
//    }

    /**
     * Set data
     */
    public void setData(List<Double> data, int minimumBinCount) {
        Variate.D d = new Variate.D(data);
        setData(d, minimumBinCount);
    }

    /**
     * Set data
     */
    public void setData(Variate.D data, int minimumBinCount) {

        this.raw = data;
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

    public Variate getRawData() {
        return raw;
    }

    public void setRawData(Variate raw) {
        this.raw = raw;
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
     * Set arbitrary intervals to use (0 for none).
     */
    public void setInCredibleSet(TraceDistribution traceD) {
        this.traceDistribution = traceD;
        hasIncredibleSet = traceD.inCredibleSet.size() > 0;
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
            lower = raw.getQuantile(quantiles);
            upper = raw.getQuantile(1.0 - quantiles);
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
                    } else if (hasIncredibleSet) {
                        if (traceDistribution.inCredibleSetContains((int) x1) || traceDistribution.inCredibleSetContains((int) x2)) {
                            g2.setPaint(quantilePaint);
                        } else {
                            g2.setPaint(barPaint);
                        }
                        fillRect(g2, x1, y1, x2, y2);
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

    protected void fillRect(Graphics2D g2, double x1, double y1, double x2, double y2) {
        if (traceDistribution != null && traceDistribution.getTraceType() != TraceFactory.TraceType.DOUBLE) {
            super.fillRect(g2, x1-(x2-x1), y1, x2, y2);
        } else {
            super.fillRect(g2, x1, y1, x2, y2);
        }
    }

    protected void drawRect(Graphics2D g2, double x1, double y1, double x2, double y2) {
        if (traceDistribution != null && traceDistribution.getTraceType() != TraceFactory.TraceType.DOUBLE) {
            super.drawRect(g2, x1-(x2-x1), y1, x2, y2);
        } else {
            super.drawRect(g2, x1, y1, x2, y2);
        }
    }

}

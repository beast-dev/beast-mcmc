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

import dr.inference.trace.TraceDistribution;
import dr.math.distributions.GammaKDEDistribution;
import dr.math.distributions.KernelDensityEstimatorDistribution;
import dr.math.distributions.LogTransformedNormalKDEDistribution;
import dr.math.distributions.NormalKDEDistribution;
import dr.stats.Variate;
import dr.util.FrequencyDistribution;

import java.awt.*;
import java.awt.geom.GeneralPath;
import java.awt.geom.Rectangle2D;
import java.util.*;

/**
 * @author Andrew Rambaut
 */
public class ViolinPlot extends KDENumericalDensityPlot { //Plot.AbstractPlot {

    private final boolean isVertical;
    private final double violinWidth;
    private final boolean showTails;
    private final boolean showQuantiles;
    private final double upper, lower;
    private double y1, y2;

    public ViolinPlot(double violinWidth, java.util.List<Double> data, int minimumBinCount) {
        this(true, violinWidth, true, data, minimumBinCount);
    }

    public ViolinPlot(boolean isVertical, double violinWidth, boolean showTails, java.util.List<Double> data) {
        this(true, violinWidth, showTails, data, DEFAULT_KDE_BINS);
    }

    public ViolinPlot(boolean isVertical, double violinWidth, boolean showTails, java.util.List<Double> data, int minimumBinCount) {
        super(data, minimumBinCount);

        this.isVertical = isVertical;
        this.violinWidth = violinWidth;
        this.showQuantiles = showTails;
        this.showTails = showTails;

        // setData again because in super, the width weren't set.
        setData(data, minimumBinCount);

        if (showQuantiles) {
            lower = getQuantile(0.025);
            upper = getQuantile(0.975);
        } else {
            lower = Collections.min(data);
            upper = Collections.max(data);
        }
    }

    public ViolinPlot(boolean isVertical, double violinWidth, double lower, double upper, boolean showTails, java.util.List<Double> data) {
        this(isVertical, violinWidth, lower, upper, showTails, data, DEFAULT_KDE_BINS);
    }

    public ViolinPlot(boolean isVertical, double violinWidth, double lower, double upper, boolean showTails, java.util.List<Double> data, int minimumBinCount) {
        super(data, minimumBinCount);

        this.isVertical = isVertical;
        this.violinWidth = violinWidth;
        this.showQuantiles = true;

        this.lower = lower;
        this.upper = upper;

        this.showTails = showTails;

        // setData again because in super, the width weren't set.
        setData(data, minimumBinCount);

    }

    /**
     * Set data
     * X data is values, y data is density
     */
    public void setData(Variate.N xData, Variate.N yData) {

        double scale = (0.5 * violinWidth) / ((Number)yData.getMax()).doubleValue();
        this.xData = new Variate.D();
        this.yData = new Variate.D();
        for (int i = 0; i < yData.getCount(); i++) {
            this.xData.add(((Number) xData.get(i)).doubleValue());
            this.yData.add(((Number) yData.get(i)).doubleValue() * scale);
        }
        for (int i = yData.getCount() - 1; i >= 0; i--) {
            this.xData.add(((Number) xData.get(i)).doubleValue());
            this.yData.add(-((Number) yData.get(i)).doubleValue() * scale);
        }

        y1 = getDensity(lower) * scale;
        y2 = getDensity(upper) * scale;
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
     * Paint data series
     */
    protected void paintData(Graphics2D g2, Variate.N xData, Variate.N yData) {

        int n = xData.getCount();

        GeneralPath path = new GeneralPath();
        GeneralPath intervalPath = new GeneralPath();

        if (isVertical) {
            // TODO this code is replicated below

            // rotate -90
            float y = (float) transformY(((Number)xData.get(0)).doubleValue());
            float x = (float) transformX(((Number)yData.get(0)).doubleValue() + getPlotNumber() + 1);

            if (showTails) {
                path.moveTo(x, y);
            }

            // create the quantile cropped path
            float v1 = (float) transformY(lower);
            float v2 = (float) transformY(upper);

            y = v1;
            x = (float) transformX(getPlotNumber() + 1);
            intervalPath.moveTo(x, y);

            x = (float) transformX(y1 + getPlotNumber() + 1);
            intervalPath.lineTo(x, y);

            if (!showTails) {
                path.moveTo(x, y);
            }

            boolean crossedAxis = false;

            for (int i = 1; i < n; i++) {
                y = (float) transformY(((Number) xData.get(i)).doubleValue());
                x = (float) transformX(((Number) yData.get(i)).doubleValue() + getPlotNumber() + 1);

                if (showTails) {
                    path.lineTo(x, y);
                }

                if (y < v1) {
                    if (y > v2) {
                        if (!showTails) {
                            path.lineTo(x, y);
                        }

                        intervalPath.lineTo(x, y);
                    } else {
                        if (!crossedAxis) {
                            y = v2;
                            x = (float) transformX(y2 + getPlotNumber() + 1);
                            intervalPath.lineTo(x, y);

                            if (!showTails) {
                                path.lineTo(x, y);
                            }

                            x = (float) transformX(-y2 + getPlotNumber() + 1);
                            intervalPath.lineTo(x, y);

                            if (!showTails) {
                                path.moveTo(x, y);
                            }

                            crossedAxis = true;
                        }
                    }
                }
            }

            // finish the quantile cropped path
            y = v1;
            x = (float) transformX(-y1 + getPlotNumber() + 1);
            intervalPath.lineTo(x, y);

            if (!showTails) {
                path.lineTo(x, y);
            }

            x = (float) transformX(getPlotNumber() + 1);
            intervalPath.lineTo(x, y);

        } else {
            // TODO this code is replicated above

            float x = (float) transformX(((Number)xData.get(0)).doubleValue());
            float y = (float) transformY(((Number)yData.get(0)).doubleValue() + getPlotNumber() + 1);

            path.moveTo(x, y);

            // create the quantile cropped path
            float v1 = (float) transformX(lower);
            float v2 = (float) transformX(upper);

            x = v1;
            y = (float) transformY(getPlotNumber() + 1);
            intervalPath.moveTo(x, y);

            y = (float) transformY(y1 + getPlotNumber() + 1);
            intervalPath.lineTo(x, y);

            boolean crossedAxis = false;

            for (int i = 1; i < n; i++) {
                x = (float) transformX(((Number) xData.get(i)).doubleValue());
                y = (float) transformY(((Number) yData.get(i)).doubleValue() + getPlotNumber() + 1);

                path.lineTo(x, y);

                if (x < v1) {
                    if (x > v2) {
                        intervalPath.lineTo(x, y);
                    } else if (!crossedAxis) {
                        x = v2;
                        y = (float) transformY(y2 + getPlotNumber() + 1);
                        intervalPath.lineTo(x, y);
                        y = (float) transformY(-y2 + getPlotNumber() + 1);
                        intervalPath.lineTo(x, y);
                        crossedAxis = true;
                    }
                }
            }

            // finish the quantile cropped path
            x = v1;
            y = (float) transformY(-y1 + getPlotNumber() + 1);
            intervalPath.lineTo(x, y);

            y = (float) transformY(getPlotNumber() + 1);
            intervalPath.lineTo(x, y);

        }

        if (solid) {
            Paint fillPaint = new Color(
                    ((Color) linePaint).getRed(),
                    ((Color) linePaint).getGreen(),
                    ((Color) linePaint).getBlue(), 32);
            g2.setPaint(fillPaint);

            if (showQuantiles) {
                intervalPath.closePath();
                g2.fill(intervalPath);
            }
        }

        if (showTails) {
            path.closePath();
        }
        
        g2.setStroke(lineStroke);
        g2.setPaint(linePaint);
        g2.draw(path);
    }


}

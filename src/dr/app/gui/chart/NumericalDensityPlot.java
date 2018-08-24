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

import dr.inference.trace.TraceType;
import dr.stats.Variate;
import dr.util.FrequencyDistribution;

import java.awt.*;
import java.awt.geom.GeneralPath;
import java.awt.geom.Rectangle2D;
import java.util.*;

public class NumericalDensityPlot extends HistogramPlot {

    boolean relativeDensity = true;
    boolean pointsOnly = false;
    int minimumBinCount;

    public boolean isSolid() {
        return solid;
    }

    public void setSolid(boolean solid) {
        this.solid = solid;
    }

    boolean solid = true;

    public NumericalDensityPlot(java.util.List<Double> data, int minimumBinCount) {
        super(data, minimumBinCount);
        this.minimumBinCount = minimumBinCount;
    }

    public NumericalDensityPlot(Variate.D data, int minimumBinCount) {
        super(data, minimumBinCount);
        this.minimumBinCount = minimumBinCount;
    }

    public void setRelativeDensity(boolean relative) {
        relativeDensity = relative;
        setData((Variate.D)getRawData(), minimumBinCount);
    }

    public void setPointsOnly(boolean pointsOnly) {
        this.pointsOnly = pointsOnly;
    }

    /**
     * Set data
     */
    public void setData(Variate.D data, int minimumBinCount) {

        setRawData(data);
        FrequencyDistribution frequency = getFrequencyDistribution(data, minimumBinCount);

        Variate.D xData = new Variate.D();
        Variate.D yData = new Variate.D();

        double x = frequency.getLowerBound() - frequency.getBinSize();
        double maxDensity = 0.0;
        for (int i = 0; i < frequency.getBinCount(); i++) {
            double density = frequency.getFrequency(i) / frequency.getBinSize() / data.getCount();
            if (density > maxDensity) maxDensity = density;
        }

        xData.add(x + (frequency.getBinSize() / 2.0));
        yData.add(0.0);
        x += frequency.getBinSize();

        for (int i = 0; i < frequency.getBinCount(); i++) {
            xData.add(x + (frequency.getBinSize() / 2.0));
            double density = frequency.getFrequency(i) / frequency.getBinSize() / data.getCount();
            if (relativeDensity) {
                yData.add(density / maxDensity);
            } else {
                yData.add(density);
            }
            x += frequency.getBinSize();
        }

        xData.add(x + (frequency.getBinSize() / 2.0));
        yData.add(0.0);

        setData(xData, yData);
    }

    /**
     * Set bar fill style. Use a barPaint of null to not fill bar.
     * Bar outline style is set using setLineStyle
     */
    public void setBarFillStyle(Paint barPaint) {
        throw new IllegalArgumentException();
    }

    /**
     * Paint data series
     */
    protected void paintData(Graphics2D g2, Variate.N xData, Variate.N yData) {

        int n = xData.getCount();

        if (pointsOnly) {
            setMarkStyle(Plot.CIRCLE_MARK, 3, new BasicStroke(0.5F), new Color(44, 44, 44), new Color(249, 202, 105));
            Rectangle2D bounds = mark.getBounds2D();
            float w = (float) bounds.getWidth();
            float h = (float) bounds.getHeight();

            for (int i = 0; i < n; i++) {
                float x = (float) transformX(((Number)xData.get(i)).doubleValue());
                float y = (float) transformY(((Number)yData.get(i)).doubleValue());

                x = x - (w / 2);
                y = y - (h / 2);

                g2.translate(x, y);

                if (markFillPaint != null) {
                    g2.setPaint(markFillPaint);
                    g2.fill(mark);
                }

                g2.setPaint(markPaint);
                g2.setStroke(markStroke);
                g2.draw(mark);

                g2.translate(-x, -y);
            }

        } else {
            float x = (float) transformX(((Number)xData.get(0)).doubleValue());
            float y = (float) transformY(((Number)yData.get(0)).doubleValue());

            GeneralPath path = new GeneralPath();
            path.moveTo(x, y);

            for (int i = 1; i < n; i++) {
                x = (float) transformX(((Number)xData.get(i)).doubleValue());
                y = (float) transformY(((Number)yData.get(i)).doubleValue());

                path.lineTo(x, y);
            }

            if (solid) {
                path.closePath();
                Paint fillPaint = new Color(
                        ((Color) linePaint).getRed(),
                        ((Color) linePaint).getGreen(),
                        ((Color) linePaint).getBlue(), 32);
                g2.setPaint(fillPaint);
                g2.fill(path);
            }

            g2.setStroke(lineStroke);
            g2.setPaint(linePaint);
            g2.draw(path);
        }
    }
}

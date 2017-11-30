/*
 * LineDistributionPlot.java
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

import dr.stats.Variate;

import java.awt.*;

public class LineDistributionPlot extends Plot.AbstractPlot {

    private Variate gradients = null;
    private Variate intercepts = null;

    private int lineCount = 0;
    private double x1, x2;
    private double[] y1;
    private double[] y2;

    private boolean isCalibrated = false;

    /**
     * Constructor
     */
    public LineDistributionPlot(Variate gradients, Variate intercepts) {
        this.gradients = gradients;
        this.intercepts = intercepts;
    }

    /**
     * Paint actual plot
     */
    @Override
    public void paintPlot(Graphics2D g2, double xScale, double yScale,
                          double xOffset, double yOffset, int plotNumber, int plotCount) {

        if (gradients == null || intercepts == null) {
            return;
        }

        super.paintPlot(g2, xScale, yScale, xOffset, yOffset, plotNumber, plotCount);

        if (!isCalibrated) {
            x1 = xAxis.getMinAxis();
            x2 = xAxis.getMaxAxis();

            lineCount = gradients.getCount();
            y1 = new double[lineCount];
            y2 = new double[lineCount];

            for (int i = 0; i < lineCount; i++) {
                y1[i] = ( (Double) gradients.get(i) * x1) + (Double) intercepts.get(i);
                y2[i] = ( (Double) gradients.get(i) * x2) + (Double) intercepts.get(i);
            }

            isCalibrated = true;
        }

        Paint paint = new Color((float) 0.1, (float) 0.1, (float) 0.1,
                (float) (50.0 / lineCount));

        g2.setPaint(paint);
        g2.setStroke(lineStroke);

        for (int i = 0; i < lineCount; i++) {
            drawLine(g2, x1, y1[i], x2, y2[i]);
        }
    }

    /**
     * Paint data series
     */
    @Override
    protected void paintData(Graphics2D g2, Variate.N xData, Variate.N yData) {
        // do nothing because paintPlot is overridden
	}
}

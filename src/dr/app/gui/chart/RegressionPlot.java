/*
 * RegressionPlot.java
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

import dr.stats.Regression;
import dr.stats.Variate;

import java.awt.*;
import java.util.List;

public class RegressionPlot extends Plot.AbstractPlot {

    private Regression regression;

    /**
     * Constructor
     */
    public RegressionPlot(Regression regression) {
        super(regression.getXData(), regression.getYData());
        this.regression = regression;
    }

    /**
     * Constructor
     */
    public RegressionPlot(Variate.N xData, Variate.N yData, boolean forceOrigin) {
        super(xData, yData);
        setForceOrigin(forceOrigin);
    }

    /**
     * Constructor
     */
    public RegressionPlot(List<Double> xData, List<Double> yData, boolean forceOrigin) {
        super(xData, yData);
        setForceOrigin(forceOrigin);
    }

    /**
     * Set data
     */
    public void setData(List<Double> xData, List<Double> yData) {
        super.setData(xData, yData);
        regression = new Regression(this.xData, this.yData);
    }

    /**
     * Set data
     */
    public void setData(Variate.N xData, Variate.N yData) {
        super.setData(xData, yData);
        regression = new Regression(this.xData, this.yData);
    }

    public void setForceOrigin(boolean forceOrigin) {
        regression.setForceOrigin(forceOrigin);
    }

    public double getGradient() {
        return regression.getGradient();
    }

    public double getYIntercept() {
        return regression.getYIntercept();
    }

    public double getXIntercept() {
        return regression.getXIntercept();
    }

    public double getResidualMeanSquared() {
        return regression.getResidualMeanSquared();
    }

    public Regression getRegression() {
        return regression;
    }

    public String toString() {
        StringBuffer statString = new StringBuffer("Gradient=");
        statString.append(Double.toString(getGradient()));
        statString.append(", Intercept=");
        statString.append(Double.toString(getYIntercept()));
        statString.append(", RMS=");
        statString.append(Double.toString(getResidualMeanSquared()));

        return statString.toString();
    }

    /**
     * Paint data series
     */
    protected void paintData(Graphics2D g2, Variate.N xData, Variate.N yData) {

        g2.setPaint(linePaint);
        g2.setStroke(lineStroke);

        double gradient = getGradient();
        double intercept = getYIntercept();

        double x1 = xAxis.getMinAxis();
        double y1 = (gradient * x1) + intercept;

        double x2 = xAxis.getMaxAxis();
        double y2 = (gradient * x2) + intercept;

        drawLine(g2, x1, y1, x2, y2);
	}
}

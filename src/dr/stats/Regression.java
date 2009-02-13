/*
 * Regression.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.stats;

import dr.util.Variate;

/**
 * simple regression of two variables
 *
 * @version $Id: Regression.java,v 1.5 2005/05/24 20:26:01 rambaut Exp $
 *
 * @author Andrew Rambaut
 */
public class Regression {
    private Variate xData = null;
    private Variate yData = null;

    private boolean forceOrigin = false;
    private boolean regressionKnown = false;

    private double gradient;
    private double intercept;
    private double covariance;
    private double residualMeanSquared;
    private double correlationCoefficient;

    /**
     * Constructor
     */
    public Regression() { }

    /**
     * Constructor
     */
    public Regression(Variate xData, Variate yData) {
        setData(xData, yData);
    }

    /**
     * Constructor
     */
    public Regression(double[] xData, double[] yData) {
        setData(xData, yData);
    }

    /**
     * Constructor
     */
    public Regression(Variate xData, Variate yData, boolean forceOrigin) {
        setData(xData, yData);
        setForceOrigin(forceOrigin);
    }

    /**
     * Constructor
     */
    public Regression(double[] xData, double[] yData, boolean forceOrigin) {
        setData(xData, yData);
        setForceOrigin(forceOrigin);
    }

    /**
     *	Set data
     */
    public void setData(double[] xData, double[] yData) {
        Variate.Double xd = new Variate.Double();
        Variate.Double yd = new Variate.Double();

        for (int i = 0; i < xData.length; i++) {
            xd.add(xData[i]);
            yd.add(yData[i]);
        }

        this.xData = xd;
        this.yData = yd;

        regressionKnown = false;
    }

    /**
     *	Set data
     */
    public void setData(Variate xData, Variate yData) {
        this.xData = xData;
        this.yData = yData;

        regressionKnown = false;
    }

    public void setForceOrigin(boolean forceOrigin) {
        this.forceOrigin = forceOrigin;

        regressionKnown = false;
    }

    public double getGradient() {
        if (!regressionKnown)
            calculateRegression();
        return gradient;
    }

    public double getIntercept() {
        if (!regressionKnown)
            calculateRegression();
        return intercept;
    }

    public double getYIntercept() {
        return getIntercept();
    }

    public double getXIntercept() {
        return -getIntercept() / getGradient();
    }

    public double getCovariance() {
        if (!regressionKnown)
            calculateRegression();
        return covariance;
    }

    public double getResidualMeanSquared() {
        if (!regressionKnown)
            calculateRegression();
        return residualMeanSquared;
    }

    public double getCorrelationCoefficient  () {
        if (!regressionKnown) {
            calculateRegression();
        }
        return correlationCoefficient;
    }

    public double getRSquared() {
        if (!regressionKnown) {
            calculateRegression();
        }
        return correlationCoefficient * correlationCoefficient;
    }

    public double getResidual(final double x, final double y) {
        return y - ((gradient * x) + intercept);
    }

    public Variate getXData() {
        return xData;
    }

    public Variate getYData() {
        return yData;
    }

    private void calculateRegression() {
        int i, n = xData.getCount();
        double meanX=0.0, meanY=0.0;

        if (!forceOrigin) {
            meanX = xData.getMean();
            meanY = yData.getMean();
        }

        //Calculate sum of products & sum of x squares

        double sumProducts=0.0;
        double sumSquareX=0.0;
        double sumSquareY=0.0;
        double x1, y1;
        for (i = 0; i < n; i++) {
            x1 = xData.get(i) - meanX;
            y1 = yData.get(i) - meanY;
            sumProducts += x1 * y1;
            sumSquareX += x1 * x1;
            sumSquareY += y1 * y1;
        }

        //Calculate gradient and intercept of regression line. Calculate covariance.

        gradient = sumProducts / sumSquareX; 			// Gradient
        intercept = meanY - (gradient * meanX);			// Intercept
        covariance = sumProducts / (n - 1);				// Covariance

        correlationCoefficient = sumProducts / Math.sqrt(sumSquareX * sumSquareY);

        //Calculate residual mean square

        double sumResidualsSquared=0;
        double residual;
        for (i = 0; i < n; i++) {
            residual = yData.get(i) - ((gradient * xData.get(i)) + intercept);
            sumResidualsSquared += residual * residual;
        }

        residualMeanSquared = sumResidualsSquared / (n - 2);// Residual Mean Square

        regressionKnown = true;
    }


}

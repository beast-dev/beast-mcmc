/*
 * Regression.java
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

package dr.stats;

/**
 * simple regression of two variables
 *
 * @author Andrew Rambaut
 * @version $Id: Regression.java,v 1.5 2005/05/24 20:26:01 rambaut Exp $
 */
public class Regression {
    private Variate.D xData = null;
    private Variate.D yData = null;

    private boolean forceOrigin = false;
    private boolean regressionKnown = false;

    private double gradient;
    private double intercept;
    private double covariance;
    private double sumResidualsSquared;
    private double residualMeanSquared;
    private double correlationCoefficient;

    /**
     * Constructor
     */
    public Regression() {
    }

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
     * Set data
     */
    public void setData(double[] xData, double[] yData) {
        Variate.D xd = new Variate.D();
        Variate.D yd = new Variate.D();

        for (int i = 0; i < xData.length; i++) {
            xd.add(xData[i]);
            yd.add(yData[i]);
        }

        this.xData = xd;
        this.yData = yd;

        regressionKnown = false;
    }

    /**
     * Set data
     */
    public void setData(Variate xData, Variate yData) {
        this.xData = (Variate.D) xData;
        this.yData = (Variate.D) yData;

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

    public double getSumResidualsSquared() {
        if (!regressionKnown)
            calculateRegression();
        return sumResidualsSquared;
    }

    public double getCorrelationCoefficient() {
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
        return y - ((getGradient() * x) + getIntercept());
    }

    public double getX(final double y) {
        return (y - getIntercept()) / getGradient();
    }

    public double getY(final double x) {
        return x * getGradient() + getIntercept();
    }

    public Variate.N getXData() {
        return xData;
    }

    public Variate.N getYData() {
        return yData;
    }

    public Variate getYResidualData() {
        Variate.D rd = new Variate.D();

        for (int i = 0; i < xData.getCount(); i++) {
            rd.add(getResidual(xData.get(i), yData.get(i)));
        }

        return rd;
    }

    private void calculateRegression() {
        int i, n = xData.getCount();
        double meanX = 0.0, meanY = 0.0;

        if (!forceOrigin) {
            meanX = xData.getMean();
            meanY = yData.getMean();
        }

        //Calculate sum of products & sum of x squares

        double sumProducts = 0.0;
        double sumSquareX = 0.0;
        double sumSquareY = 0.0;
        double x1, y1;
        for (i = 0; i < n; i++) {
            x1 = xData.get(i) - meanX;
            y1 = yData.get(i) - meanY;
            sumProducts += x1 * y1;
            sumSquareX += x1 * x1;
            sumSquareY += y1 * y1;
        }

        //Calculate gradient and intercept of regression line. Calculate covariance.

        gradient = sumProducts / sumSquareX;             // Gradient
        intercept = meanY - (gradient * meanX);            // Intercept
        covariance = sumProducts / (n - 1);                // Covariance

        correlationCoefficient = sumProducts / Math.sqrt(sumSquareX * sumSquareY);

        //Calculate residual mean square

        sumResidualsSquared = 0;
        double residual;
        for (i = 0; i < n; i++) {
            residual = yData.get(i) - ((gradient * xData.get(i)) + intercept);
            sumResidualsSquared += residual * residual;
        }

        residualMeanSquared = sumResidualsSquared / (n - 2);// Residual Mean Square

        regressionKnown = true;
    }


}

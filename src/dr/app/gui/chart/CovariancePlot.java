/*
 * CovariancePlot.java
 *
 * Copyright (c) 2002-2017 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
import java.util.Random;

/**
 * Description:	A covariance ellipse plot.
 *
 * @author Guy Baele
 */
public class CovariancePlot extends Plot.AbstractPlot {

    private int plotCount;

    public CovariancePlot(java.util.List<Double> data, int minimumBinCount) {
        super(data, data);
        System.out.println("CovariancePlot: " + data.size());
        setName("null");
        this.plotCount = 1;
    }

    public CovariancePlot(java.util.List<Double> xData, java.util.List<Double> yData) {
        super(xData, yData);
        System.out.println("xData: " + xData.size());
        System.out.println("yData: " + yData.size());
        setName("null");
        this.plotCount = 1;
    }

    public CovariancePlot(String name, java.util.List<Double> xData, java.util.List<Double> yData) {
        super(xData, yData);
        setName(name);
        this.plotCount = 1;
    }

    public void setTotalPlotCount(int count) {
        this.plotCount = count;
    }

    /**
     * Paint data series
     */
    protected void paintData(Graphics2D g2, Variate.N xData, Variate.N yData) {

        //System.out.println("CovariancePlot: paintData");
        //System.out.println("PlotNumber = " + plotNumber);
        //System.out.println("TotalPlotCount = " + this.plotCount);

        //TODO replace random coloring with colors based on correlation value
        Random rand = new Random();
        Color fillColor = new Color(rand.nextInt(255), rand.nextInt(255), rand.nextInt(255));

        g2.setColor(fillColor);

        //g2.setPaint(linePaint);
        //g2.setStroke(lineStroke);

        double x1 = (plotNumber/(int)Math.sqrt(plotCount))+0.75;
        double y1 = (plotNumber%(int)Math.sqrt(plotCount))+0.75;
        double x2 = (plotNumber/(int)Math.sqrt(plotCount))+1.25;
        double y2 = (plotNumber%(int)Math.sqrt(plotCount))+1.25;

        //System.out.println("(" + x1 + "," + y1 + ") (" + x2 + "," + y2 + ")");

        //TODO draw actual ellipse plots
        drawRect(g2, x1, y1, x2, y2);
        g2.drawString("plot" + plotNumber, (float)transformX(x1), (float)transformY(y1));

    }

}

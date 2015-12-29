/*
 * AreaPlot.java
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
import java.awt.geom.GeneralPath;
import java.util.List;

/**
 * Description:	An area plot.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: AreaPlot.java,v 1.4 2005/05/24 20:25:59 rambaut Exp $
 */

public class AreaPlot extends Plot.AbstractPlot {

    protected Variate.N xData2 = null;
    protected Variate.N yData2 = null;

    /**
     * Constructor
     */
    public AreaPlot(Variate.N xData, Variate.N yData) {
        super(xData, yData);
    }

    /**
     * Constructor
     */
    public AreaPlot(List<Double> xData, List<Double> yData) {
        super(xData, yData);
    }

    /**
     * Constructor
     */
    public AreaPlot(Variate.N xData1, Variate.N yData1, Variate.N xData2, Variate.N yData2) {
        super(xData1, yData1);
        this.xData2 = xData2;
        this.yData2 = yData2;
    }

    /**
     * Constructor
     */
    public AreaPlot(List<Double> xData1, List<Double> yData1, List<Double> xData2, List<Double> yData2) {
        super(xData1, yData1);
        this.xData2 = new Variate.D(xData2);
        this.yData2 = new Variate.D(yData2);
    }

    /**
     * Set data
     */
    public void setData(List<Double> xData1, List<Double> yData1, List<Double> xData2, List<Double> yData2) {
        setData(xData1, yData1);
        this.xData2 = new Variate.D(xData2);
        this.yData2 = new Variate.D(yData2);
    }

    /**
     * Set data
     */
    public void setData(Variate.N xData1, Variate.N yData1, Variate.N xData2, Variate.N yData2) {
        setData(xData1, yData1);
        this.xData2 = xData2;
        this.yData2 = yData2;
    }

    /**
     * Set axes
     */
    public void setAxes(Axis xAxis, Axis yAxis) {
        super.setAxes(xAxis, yAxis);

        if (xData2 != null && yData2 != null) {
            setupAxis(xAxis, yAxis, xData2, yData2);
        }
    }

    /**
     * Resets axis ranges (if new data has been added)
     */
    public void resetAxes() {
        super.resetAxes();

        if (xData2 != null && yData2 != null) {
            setupAxis(xAxis, yAxis, xData2, yData2);
        }
    }

    /**
     * Paint data series
     */
    protected void paintData(Graphics2D g2, Variate.N xData, Variate.N yData) {

        double x0 = transformX(((Number) xData.get(0)).doubleValue());
        double y0 = transformY(((Number) yData.get(0)).doubleValue());

        GeneralPath path = new GeneralPath();
        path.moveTo((float) x0, (float) y0);

        double x = x0;
        double y = y0;

        for (int i = 1, n = xData.getCount(); i < n; i++) {
            x = transformX(((Number) xData.get(i)).doubleValue());
            y = transformY(((Number) yData.get(i)).doubleValue());
            path.lineTo((float) x, (float) y);
        }

        if (xData2 != null & yData2 != null) {
            for (int i = xData2.getCount() - 1; i >= 0; i--) {
                x = transformX(((Number) xData2.get(i)).doubleValue());
                y = transformY(((Number) yData2.get(i)).doubleValue());
                path.lineTo((float) x, (float) y);
            }

        } else {
            double y1 = transformY(0.0);
            path.lineTo((float) x, (float) y1);
            path.lineTo((float) x0, (float) y1);
            path.lineTo((float) x0, (float) y0);
        }

        path.closePath();

        g2.setPaint(linePaint);
        g2.setStroke(lineStroke);

		g2.fill(path);
	}

}


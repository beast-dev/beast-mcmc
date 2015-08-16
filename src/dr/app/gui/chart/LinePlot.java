/*
 * LinePlot.java
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
 * Description:	A line plot.
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: LinePlot.java,v 1.8 2005/05/24 20:25:59 rambaut Exp $
 */

public class LinePlot extends Plot.AbstractPlot {


    /**
     * Constructor
     */
    public LinePlot(Variate.N xData, Variate.N yData) {
        super(xData, yData);
    }

    /**
     * Constructor
     */
    public LinePlot(List<Double> xData, List<Double> yData) {
        super(xData, yData);
    }

    /**
     * Paint data series
     */
    protected void paintData(Graphics2D g2, Variate.N xData, Variate.N yData) {

        double x = transformX(((Number) xData.get(0)).doubleValue());
        double y = transformY(((Number) yData.get(0)).doubleValue());

        GeneralPath path = new GeneralPath();
        path.moveTo((float) x, (float) y);

        int n = xData.getCount();
        boolean failed = false;
        for (int i = 1; i < n; i++) {
            x = transformX(((Number) xData.get(i)).doubleValue());
            y = transformY(((Number) yData.get(i)).doubleValue());
            if (x == Double.NEGATIVE_INFINITY || y == Double.NEGATIVE_INFINITY ||
                    Double.isNaN(x) || Double.isNaN(y)) {
                failed = true;
            } else if (failed) {
                failed = false;
                path.moveTo((float) x, (float) y);
            } else {
                path.lineTo((float) x, (float) y);
            }
        }

        g2.setPaint(linePaint);
        g2.setStroke(lineStroke);

        g2.draw(path);
	}

}


/*
 * ErrorBarPlot.java
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
 * Description:	An error bar plot (to be put under a scatter plot).
 *
 * @author Andrew Rambaut
 * @version $Id: LinePlot.java,v 1.8 2005/05/24 20:25:59 rambaut Exp $
 */

public class ErrorBarPlot extends Plot.AbstractPlot {
         public enum Orientation {
             HORIZONTAL,
             VERTICAL
         }

    /**
     * Constructor
     */
    public ErrorBarPlot(Orientation orientation, Variate.N xData, Variate.N yData, Variate.N errorData) {
        super(xData, yData);

        this.errorData = errorData;
        this.orientation = orientation;
    }

    /**
     * Constructor
     */
    public ErrorBarPlot(Orientation orientation, List<Double> xData, List<Double> yData, List<Double> errorData) {
        super(xData, yData);

        this.errorData = new Variate.D(errorData);
        this.orientation = orientation;
    }

    /**
     * Paint data series
     */
    protected void paintData(Graphics2D g2, Variate.N xData, Variate.N yData) {

        g2.setPaint(linePaint);
        g2.setStroke(lineStroke);

        int n = xData.getCount();
        for (int i = 0; i < n; i++) {
            GeneralPath path = new GeneralPath();

            double x0 = ((Number) xData.get(i)).doubleValue();
            double y0 = ((Number) yData.get(i)).doubleValue();
            double e = ((Number) errorData.get(i)).doubleValue() / 2;

            float fx = (float)transformX(x0);
            float fy = (float)transformY(y0);

            if (!Double.isInfinite(fx) && !Double.isInfinite(fy) &&
                    !Double.isNaN(fx) && !Double.isNaN(fy)) {
                if (orientation == Orientation.HORIZONTAL) {
                    float fx1 = (float)transformX(x0 - e);
                    float fx2 = (float)transformX(x0 + e);

                    path.moveTo(fx, fy);
                    path.lineTo(fx1, fy);
                    path.moveTo(fx, fy);
                    path.lineTo(fx2, fy);
                } else if (orientation == Orientation.VERTICAL) {
                    float fy1 = (float)transformY(y0 - e);
                    float fy2 = (float)transformY(y0 + e);

                    path.moveTo(fx, fy);
                    path.lineTo(fx, fy1);
                    path.moveTo(fx, fy);
                    path.lineTo(fx, fy2);
                }
            }
            g2.draw(path);
        }


	}

    private Orientation orientation;
    protected Variate.N errorData = null;
}


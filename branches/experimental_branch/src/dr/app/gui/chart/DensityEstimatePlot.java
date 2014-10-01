/*
 * DensityPlot.java
 *
 * Copyright (C) 2002-2009 Alexei Drummond and Andrew Rambaut
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
 * BEAST is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.gui.chart;

import dr.inference.trace.DensityEstimate;
import dr.stats.Variate;

import java.awt.*;
import java.awt.geom.GeneralPath;

public class DensityEstimatePlot extends Plot.AbstractPlot {

    public boolean isSolid() {
        return solid;
    }

    public void setSolid(boolean solid) {
        this.solid = solid;
    }

    boolean solid = true;

    public DensityEstimatePlot(DensityEstimate densityEstimate) {
        setData(densityEstimate.getXCoordinates(), densityEstimate.getYCoordinates());
    }


    public void setRelativeDensity(boolean relative) {
//        relativeDensity = relative;
//        setData(getRawData(), minimumBinCount);
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
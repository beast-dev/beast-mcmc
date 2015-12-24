/*
 * ParentPlot.java
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

package dr.app.tempest;

import dr.app.gui.chart.Plot;
import dr.stats.Variate;

import java.awt.*;
import java.awt.geom.GeneralPath;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;

/**
 * Description:	A line plot.
 *
 * @author Andrew Rambaut
 * @version $Id$
 */

public class ParentPlot extends Plot.AbstractPlot {


    /**
     * Constructor
     */
    public ParentPlot(Variate xData, Variate yData, List<Double> xParentData, List<Double> yParentData) {
        super(xParentData, yParentData);

        this.xTipData = xData;
        this.yTipData = yData;

        this.xParentData = new Variate.D(xParentData);
        this.yParentData = new Variate.D(yParentData);
    }

    /**
     * Paint data series
     */
    protected void paintData(Graphics2D g2, Variate.N xData, Variate.N yData) {

        g2.setPaint(linePaint);
        g2.setStroke(lineStroke);

        if (getSelectedPoints() != null && getSelectedPoints().size() > 0) {
            for (int i : getSelectedPoints()) {

                double x = ((Number) xTipData.get(i)).doubleValue();
                double y = ((Number) yTipData.get(i)).doubleValue();

                double x1 = transformX(x);
                double y1 = transformY(y);

                double x2 = transformX(((Number) xData.get(0)).doubleValue());
                double y2 = transformY(((Number) yData.get(0)).doubleValue());

                GeneralPath path = new GeneralPath();
                path.moveTo((float) x1, (float) y1);
//            path.lineTo((float) x2, (float) y1);
                path.lineTo((float) x2, (float) y2);

                g2.draw(path);
            }
        } else {
        for (int i = 0; i < xData.getCount(); i++) {

            double x1 = transformX(((Number) xTipData.get(i)).doubleValue());
            double y1 = transformY(((Number) yTipData.get(i)).doubleValue());

            double x2 = transformX(((Number) xData.get(i)).doubleValue());
            double y2 = transformY(((Number) yData.get(i)).doubleValue());

            GeneralPath path = new GeneralPath();
            path.moveTo((float) x1, (float) y1);
//            path.lineTo((float) x2, (float) y1);
            path.lineTo((float) x2, (float) y2);

            g2.draw(path);
        }
        }


	}

    private final Variate xTipData;
    private final Variate yTipData;

    private final Variate.N xParentData;
    private final Variate.N yParentData;

    public void setSelectedPoints(Set<Integer> selectedPoints, double mrcaTime, double mrcaDistance) {
        List<Double> x = new ArrayList<Double>();
        x.add(mrcaTime);
        List<Double> y = new ArrayList<Double>();
        y.add(mrcaDistance);
        setData(x, y);
        setSelectedPoints(selectedPoints);
    }

    public void clearSelection() {
        setData(xParentData, yParentData);
        super.clearSelection();
    }

}


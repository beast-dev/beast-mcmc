/*
 * ParallelChart.java
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

import dr.app.gui.chart.Axis;
import dr.app.gui.chart.DiscreteAxis;
import dr.app.gui.chart.JChart;
import dr.app.gui.chart.Plot;

/**
 * A chart that has the primary dimension vertically and parallel plots horizonally (i.e.,
 * for box plots, intervals and violins.
 *
 * @author Andrew Rambaut
 */
public class JParallelChart extends JChart {

    private final boolean horizontal;
    private final DiscreteAxis parallelAxis;

    public JParallelChart(Axis dimensionAxis) {
        this(false, dimensionAxis);
    }

    public JParallelChart(boolean horizontal, Axis dimensionAxis) {
        super(null,null);

        this.horizontal = horizontal;

        parallelAxis = new DiscreteAxis(true, true);
        setXAxis(horizontal ? dimensionAxis : parallelAxis);
        setYAxis(horizontal ? parallelAxis : dimensionAxis);
    }

    @Override
    public void addPlot(Plot plot) {
        super.addPlot(plot);

        // set the range manually to the number of plots
        parallelAxis.setRange(1.0, getPlotCount());
        recalibrate();
        repaint();
    }

    @Override
    protected String getXAxisLabel(double value) {
        if (horizontal) {
            return super.getXAxisLabel(value);
        } else {
            return getPlotName(value);
        }
    }

    @Override
    protected String getYAxisLabel(double value) {
        if (!horizontal) {
            return super.getYAxisLabel(value);
        } else {
            return getPlotName(value);
        }
    }


    private String getPlotName(double value) {
        int index = (int)(value);
        if (index >= 1 && index <= getPlotCount()) {
            Plot plot = getPlot(index - 1);
            return plot.getName();
        } else {
            return "";
        }

    }

}
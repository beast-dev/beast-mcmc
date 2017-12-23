/*
 * JGridChart.java
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

import dr.app.gui.chart.*;

/**
 * An NxK grid of individual plots
 * @author Guy Baele
 */
public class JGridChart extends JChart {

    private final CustomAxis xVariableAxis;
    private final CustomAxis yVariableAxis;

    private int rowCount;
    private int columnCount;

    public JGridChart(double aspectRatio) {
        super(null, null, aspectRatio);

        xVariableAxis = new CustomAxis(1, 2);
        yVariableAxis = new CustomAxis(1, 2);

        //this sets the axes on the JGridChart
        setXAxis(xVariableAxis);
        setYAxis(yVariableAxis);

        // set the row and column counts to the sqrt of the number of plots
        rowCount = 0;
        columnCount = 0;
    }

    public JGridChart(int rowCount, int columnCount, double aspectRatio) {
        super(null, null, aspectRatio);

        this.rowCount = rowCount;
        this.columnCount = columnCount;

        xVariableAxis = new CustomAxis(1, 2);
        yVariableAxis = new CustomAxis(1, 2);

        //this sets the axes on the JGridChart
        setXAxis(xVariableAxis);
        setYAxis(yVariableAxis);

    }

    public void setDimensions(int rowCount, int columnCount) {
        this.rowCount = rowCount;
        this.columnCount = columnCount;

        xVariableAxis.setRange(1.0, columnCount);
        yVariableAxis.setRange(1.0, rowCount);

        recalibrate();
        repaint();
    }


    @Override
    public void addPlot(Plot plot) {

        super.addPlot(plot);

        int rowCount = this.rowCount;
        int columnCount = this.columnCount;

        if (rowCount < 1 || columnCount < 1) {
            // set the range manually to the square root of the number of plots
            int k = (int)Math.sqrt(getPlotCount());
            rowCount = getPlotCount() / k;
            columnCount = getPlotCount() % k;

            if (plot.getXLocation() < 0 || plot.getYLocation() < 0) {
                int x = 0;
                int y = 0;
                for (Plot p : getPlots()) {
                    p.setLocation(x, y);
                    x++;
                    if (x >= columnCount) {
                        x = 0;
                        y++;
                    }
                }
            }
        }

        xVariableAxis.setRange(1.0, columnCount);
        yVariableAxis.setRange(1.0, rowCount);


        recalibrate();
        repaint();
    }

    @Override
    protected String getXAxisLabel(double value) {
        return getPlotName(value);
    }

    @Override
    protected String getYAxisLabel(double value) {
        return getPlotName(value);
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

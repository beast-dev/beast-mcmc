/*
 * MapperPanel.java
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

package dr.app.mapper.application.mapper;

import dr.app.gui.chart.*;
import dr.app.mapper.application.MapperFrameOld;
import dr.inference.trace.TraceList;
import jam.framework.Exportable;

import javax.swing.*;
import java.awt.*;
import java.util.*;

/**
 * @author Andrew Rambaut
 * @version $Id$
 */
public class MapperPanel extends JPanel implements Exportable {

    private int rowCount = 2;
    private int columnCount = 2;

    JChart[][] charts = new JChart[rowCount][columnCount];
    Axis[] rowAxis = new Axis[rowCount];
    Axis[] columnAxis = new Axis[columnCount];

    TraceList[] traceLists;
    java.util.List<Layer> layers;

    public MapperPanel(MapperFrameOld mapperFrame) {
        for (int column = 0; column < columnCount; column++) {
            columnAxis[column] = new LinearAxis();
        }

        JPanel panel1 = new JPanel();
        panel1.setLayout(new GridBagLayout());

        GridBagConstraints gc = new GridBagConstraints();
        gc.weightx = 1.0;
        gc.weighty = 1.0;
        gc.fill = GridBagConstraints.BOTH;

        for (int row = 0; row < rowCount; row++) {
            rowAxis[row] = new LinearAxis();

            for (int column = 0; column < columnCount; column++) {
                charts[row][column] = new JChart(rowAxis[row], columnAxis[column]);

                charts[row][column].setBorder(BorderFactory.createLineBorder(Color.black));
                gc.gridx = column;
                gc.gridy = row;
                panel1.add(charts[row][column], gc);
            }
        }

        setLayout(new BorderLayout());
        add(panel1, BorderLayout.CENTER);
    }

    /**
     * This function takes a multiple statistics in a single log files
     */
    public void setTraces(TraceList[] traceLists) {
        this.traceLists = traceLists;

        setupGraphs();
    }

    /**
     * This function takes a multiple statistics in a single log files
     */
    public void setLayers(java.util.List<Layer> layers) {
        this.layers = layers;

        setupGraphs();
    }


   private void setupGraphs() {
       for (Layer layer : layers) {
           for (int row = 0; row < rowCount; row++) {
               for (int column = 0; column < columnCount; column++) {
                   java.util.List<Double> xData = new ArrayList<Double>();
                   java.util.List<Double> yData = new ArrayList<Double>();

                   Plot plot = new ScatterPlot(xData, yData);
                   charts[row][column].addPlot(plot);
               }
           }
       }
    }

    public void doCopy() {

        java.awt.datatransfer.Clipboard clipboard =
                Toolkit.getDefaultToolkit().getSystemClipboard();

        java.awt.datatransfer.StringSelection selection =
                new java.awt.datatransfer.StringSelection(getExportText());

        clipboard.setContents(selection, selection);

    }

    public String getExportText() {
        return "";
    }


   @Override
    public JComponent getExportableComponent() {
        return this;
    }
}

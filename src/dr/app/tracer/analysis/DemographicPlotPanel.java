/*
 * DemographicPlotPanel.java
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

package dr.app.tracer.analysis;

import dr.app.gui.chart.*;
import dr.stats.Variate;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;


/**
 * A panel that displays demographic plot
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: DemographicPlotPanel.java,v 1.4 2006/11/22 23:49:37 rambaut Exp $
 */
public class DemographicPlotPanel extends JPanel {

    private JChart demoChart = new JChart(new LinearAxis(Axis.AT_DATA, Axis.AT_DATA), new LogAxis());
    private JChartPanel chartPanel = new JChartPanel(demoChart, null, "", "");

    private JComboBox meanMedianComboBox = new JComboBox(new String[]{"Median", "Mean"});
    private JCheckBox solidIntervalCheckBox = new JCheckBox("Solid interval");

    private ChartSetupDialog chartSetupDialog = null;

    private Variate.D xData = null;
    private Variate.D yDataMean = null;
    private Variate.D yDataMedian = null;
    private Variate.D yDataUpper = null;
    private Variate.D yDataLower = null;

    private double timeMedian = -1;
    private double timeMean = -1;
    private double timeUpper = -1;
    private double timeLower = -1;

    /**
     * Creates new DemographicPlotPanel
     *
     * @param frame the parent frame
     */
    public DemographicPlotPanel(final JFrame frame) {
        setMinimumSize(new Dimension(300, 150));
        setLayout(new BorderLayout());

        JToolBar toolBar = new JToolBar();
        toolBar.setLayout(new FlowLayout(FlowLayout.LEFT));
        toolBar.setFloatable(false);

        JLabel label = new JLabel("Show:");
        label.setLabelFor(meanMedianComboBox);
        toolBar.add(label);

        toolBar.add(meanMedianComboBox);

        JButton chartSetupButton = new JButton("Setup Axes");
        toolBar.add(new JToolBar.Separator(new Dimension(8, 8)));
        toolBar.add(chartSetupButton);

        toolBar.add(new JToolBar.Separator(new Dimension(8, 8)));
        toolBar.add(solidIntervalCheckBox);

        toolBar.add(new JToolBar.Separator(new Dimension(8, 8)));

        meanMedianComboBox.setFont(UIManager.getFont("SmallSystemFont"));

        add(chartPanel, BorderLayout.CENTER);
        add(toolBar, BorderLayout.SOUTH);

        meanMedianComboBox.addItemListener(
                new java.awt.event.ItemListener() {
                    public void itemStateChanged(java.awt.event.ItemEvent ev) {
                        updatePlots();
                    }
                }
        );

        chartSetupButton.addActionListener(
                new java.awt.event.ActionListener() {
                    public void actionPerformed(ActionEvent actionEvent) {
                        if (chartSetupDialog == null) {
                            chartSetupDialog = new ChartSetupDialog(frame, false, true,
                                    Axis.AT_DATA, Axis.AT_DATA, Axis.AT_DATA, Axis.AT_DATA);
                        }

                        chartSetupDialog.showDialog(demoChart);
                        validate();
                        repaint();
                    }
                }
        );

        solidIntervalCheckBox.addItemListener(
                new java.awt.event.ItemListener() {
                    public void itemStateChanged(java.awt.event.ItemEvent ev) {
                        updatePlots();
                    }
                }
        );

    }

    public void setupPlot(String title, Variate.D xData,
                          Variate.D yDataMean, Variate.D yDataMedian,
                          Variate.D yDataUpper, Variate.D yDataLower,
                          double timeMean, double timeMedian,
                          double timeUpper, double timeLower) {

        this.xData = xData;
        this.yDataMean = yDataMean;
        this.yDataMedian = yDataMedian;
        this.yDataUpper = yDataUpper;
        this.yDataLower = yDataLower;

        this.timeMean = timeMean;
        this.timeMedian = timeMedian;
        this.timeUpper = timeUpper;
        this.timeLower = timeLower;

        if (xData == null) {
            demoChart.removeAllPlots();
            chartPanel.setTitle("");
            chartPanel.setXAxisTitle("");
            chartPanel.setYAxisTitle("");
            return;
        }

        chartPanel.setTitle(title);

        updatePlots();

        chartPanel.setXAxisTitle("Time");
        chartPanel.setYAxisTitle("Population Size");
    }

    public void updatePlots() {

        demoChart.removeAllPlots();

        if (solidIntervalCheckBox.isSelected()) {

            AreaPlot areaPlot = new AreaPlot(xData, yDataUpper, xData, yDataLower);
            areaPlot.setLineColor(new Color(0x9999FF));
            demoChart.addPlot(areaPlot);
        } else {
            LinePlot plot = new LinePlot(xData, yDataLower);
            plot.setLineStyle(new BasicStroke(1.0F), new Color(0x9999FF));
            demoChart.addPlot(plot);

            plot = new LinePlot(xData, yDataUpper);
            plot.setLineStyle(new BasicStroke(1.0F), new Color(0x9999FF));
            demoChart.addPlot(plot);

        }

        LinePlot linePlot;
        if (meanMedianComboBox.getSelectedItem().equals("Median")) {
            linePlot = new LinePlot(xData, yDataMedian);
        } else {
            linePlot = new LinePlot(xData, yDataMean);
        }
        linePlot.setLineStyle(new BasicStroke(2.0F), Color.black);
        demoChart.addPlot(linePlot);

        Variate.D y1 = new Variate.D();
        y1.add(demoChart.getYAxis().getMinAxis());
        y1.add(demoChart.getYAxis().getMaxAxis());

        if (timeMean > 0.0 && timeMedian > 0.0) {
            Variate.D x1 = new Variate.D();
            if (meanMedianComboBox.getSelectedItem().equals("Median")) {
                x1.add(timeMedian);
                x1.add(timeMedian);

            } else {
                x1.add(timeMean);
                x1.add(timeMean);
            }
            LinePlot linePlot2 = new LinePlot(x1, y1);
            linePlot2.setLineStyle(new BasicStroke(2F, BasicStroke.CAP_ROUND, BasicStroke.JOIN_BEVEL, 0.0F,
                    new float[]{0.5F, 3.0F}, 0.0F), Color.black);
            demoChart.addPlot(linePlot2);
        }

        if (timeLower > 0.0) {
            Variate.D x2 = new Variate.D();
            x2.add(timeLower);
            x2.add(timeLower);

            LinePlot linePlot3 = new LinePlot(x2, y1);
            linePlot3.setLineStyle(new BasicStroke(1.0F, BasicStroke.CAP_ROUND, BasicStroke.JOIN_BEVEL, 0.0F,
                    new float[]{0.5F, 2.0F}, 0.0F), Color.black);
            demoChart.addPlot(linePlot3);
        }

        if (timeUpper > 0.0) {
            Variate.D x3 = new Variate.D();
            x3.add(timeUpper);
            x3.add(timeUpper);

            LinePlot linePlot4 = new LinePlot(x3, y1);
            linePlot4.setLineStyle(new BasicStroke(1.0F, BasicStroke.CAP_ROUND, BasicStroke.JOIN_BEVEL, 0.0F,
                    new float[]{0.5F, 2.0F}, 0.0F), Color.black);
            demoChart.addPlot(linePlot4);
        }

        validate();
        repaint();
    }

    public JComponent getExportableComponent() {
        return chartPanel;
    }

}

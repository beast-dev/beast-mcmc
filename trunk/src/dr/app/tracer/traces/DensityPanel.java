/*
 * DensityPanel.java
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

package dr.app.tracer.traces;

import dr.gui.chart.*;
import dr.inference.trace.*;
import dr.stats.Variate;
import org.virion.jam.framework.Exportable;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.util.HashMap;
import java.util.Map;

/**
 * A panel that displays density plots of traces
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: DensityPanel.java,v 1.3 2006/11/29 09:54:30 rambaut Exp $
 */
public class DensityPanel extends JPanel implements Exportable {
    public static int COLOUR_BY_TRACE = 0;
    public static int COLOUR_BY_FILE = 1;
    public static int COLOUR_BY_ALL = 2;

    private static final Paint[] paints = new Paint[]{
            Color.BLACK,
            new Color(64, 35, 225),
            new Color(229, 35, 60),
            new Color(255, 174, 34),
            new Color(86, 255, 34),
            new Color(35, 141, 148),
            new Color(146, 35, 142),
            new Color(255, 90, 34),
            new Color(239, 255, 34),
            Color.DARK_GRAY
    };

    private ChartSetupDialog chartSetupDialog = null;

//    private JChart traceChart = new JChart(new LinearAxis(Axis.AT_MAJOR_TICK, Axis.AT_MAJOR_TICK), new LinearAxis());
    private DiscreteJChart traceChart = new DiscreteJChart(new LinearAxis(Axis.AT_MAJOR_TICK, Axis.AT_MAJOR_TICK), new LinearAxis());

    private JChartPanel chartPanel = new JChartPanel(traceChart, null, "", "");

    private int minimumBins = 100;
    private JLabel labelBins;
    private JComboBox binsCombo = new JComboBox(
            new Integer[]{10, 20, 50, 100, 200, 500, 1000});

    private JCheckBox relativeDensityCheckBox = new JCheckBox("Relative density");
    private JCheckBox solidCheckBox = new JCheckBox("Fill plot");
    private JComboBox legendCombo = new JComboBox(
            new String[]{"None", "Top-Left", "Top", "Top-Right", "Left",
                    "Right", "Bottom-Left", "Bottom", "Bottom-Right"}
    );
    private JComboBox colourByCombo = new JComboBox(
            new String[]{"Trace", "Trace File", "All"}
    );
    private JLabel messageLabel = new JLabel("No data loaded");

    private int colourBy = COLOUR_BY_TRACE;

    private final JFrame frame;

    /**
     * Creates new FrequencyPanel
     */
    public DensityPanel(final JFrame frame) {
        this.frame = frame;

        setOpaque(false);

        setMinimumSize(new Dimension(300, 150));
        setLayout(new BorderLayout());

        JToolBar toolBar = new JToolBar();
        toolBar.setOpaque(false);
        toolBar.setLayout(new FlowLayout(FlowLayout.LEFT));
        toolBar.setFloatable(false);

        JButton chartSetupButton = new JButton("Axes...");
        chartSetupButton.putClientProperty(
                "Quaqua.Button.style", "placard"
        );
        chartSetupButton.setFont(UIManager.getFont("SmallSystemFont"));
        toolBar.add(chartSetupButton);

        binsCombo.setFont(UIManager.getFont("SmallSystemFont"));
        binsCombo.setOpaque(false);
        binsCombo.setSelectedItem(minimumBins);
        labelBins = new JLabel("Bins:");
        labelBins.setFont(UIManager.getFont("SmallSystemFont"));
        labelBins.setLabelFor(binsCombo);
        toolBar.add(labelBins);
        toolBar.add(binsCombo);

        relativeDensityCheckBox.setOpaque(false);
        relativeDensityCheckBox.setFont(UIManager.getFont("SmallSystemFont"));
        toolBar.add(relativeDensityCheckBox);

        // Probably don't need this as an option - takes up space and
        // solid (translucent) plots look cool...
//		solidCheckBox.setOpaque(false);
//		solidCheckBox.setFont(UIManager.getFont("SmallSystemFont"));
//		solidCheckBox.setSelected(true);
//		toolBar.add(solidCheckBox);

        toolBar.add(new JToolBar.Separator(new Dimension(8, 8)));
        JLabel label = new JLabel("Legend:");
        label.setFont(UIManager.getFont("SmallSystemFont"));
        label.setLabelFor(legendCombo);
        toolBar.add(label);
        legendCombo.setFont(UIManager.getFont("SmallSystemFont"));
        legendCombo.setOpaque(false);
        toolBar.add(legendCombo);

        toolBar.add(new JToolBar.Separator(new Dimension(8, 8)));
        label = new JLabel("Colour by:");
        label.setFont(UIManager.getFont("SmallSystemFont"));
        label.setLabelFor(colourByCombo);
        toolBar.add(label);
        colourByCombo.setFont(UIManager.getFont("SmallSystemFont"));
        colourByCombo.setOpaque(false);
        toolBar.add(colourByCombo);

        toolBar.add(new JToolBar.Separator(new Dimension(8, 8)));

        add(messageLabel, BorderLayout.NORTH);
        add(toolBar, BorderLayout.SOUTH);
        add(chartPanel, BorderLayout.CENTER);

        chartSetupButton.addActionListener(
                new java.awt.event.ActionListener() {
                    public void actionPerformed(ActionEvent actionEvent) {
                        if (chartSetupDialog == null) {
                            chartSetupDialog = new ChartSetupDialog(frame, true, false,
                                    Axis.AT_MAJOR_TICK, Axis.AT_MAJOR_TICK, Axis.AT_ZERO, Axis.AT_MAJOR_TICK);
                        }

                        chartSetupDialog.showDialog(traceChart);
                        validate();
                        repaint();
                    }
                }
        );

        binsCombo.addItemListener(
                new java.awt.event.ItemListener() {
                    public void itemStateChanged(java.awt.event.ItemEvent ev) {
                        minimumBins = (Integer) binsCombo.getSelectedItem();
                        setupTraces();
                    }
                }
        );

        relativeDensityCheckBox.addItemListener(
                new java.awt.event.ItemListener() {
                    public void itemStateChanged(java.awt.event.ItemEvent ev) {
                        for (int i = 0; i < traceChart.getPlotCount(); i++) {
                            ((NumericalDensityPlot) traceChart.getPlot(i)).setRelativeDensity(relativeDensityCheckBox.isSelected());
                        }
                        traceChart.recalibrate();
                        validate();
                        repaint();
                    }
                }
        );

        solidCheckBox.addItemListener(
                new java.awt.event.ItemListener() {
                    public void itemStateChanged(java.awt.event.ItemEvent ev) {
                        for (int i = 0; i < traceChart.getPlotCount(); i++) {
                            ((NumericalDensityPlot) traceChart.getPlot(i)).setSolid(solidCheckBox.isSelected());
                        }
                        traceChart.recalibrate();
                        validate();
                        repaint();
                    }
                }
        );

        legendCombo.addItemListener(
                new java.awt.event.ItemListener() {
                    public void itemStateChanged(java.awt.event.ItemEvent ev) {
                        switch (legendCombo.getSelectedIndex()) {
                            case 0:
                                break;
                            case 1:
                                traceChart.setLegendAlignment(SwingConstants.NORTH_WEST);
                                break;
                            case 2:
                                traceChart.setLegendAlignment(SwingConstants.NORTH);
                                break;
                            case 3:
                                traceChart.setLegendAlignment(SwingConstants.NORTH_EAST);
                                break;
                            case 4:
                                traceChart.setLegendAlignment(SwingConstants.WEST);
                                break;
                            case 5:
                                traceChart.setLegendAlignment(SwingConstants.EAST);
                                break;
                            case 6:
                                traceChart.setLegendAlignment(SwingConstants.SOUTH_WEST);
                                break;
                            case 7:
                                traceChart.setLegendAlignment(SwingConstants.SOUTH);
                                break;
                            case 8:
                                traceChart.setLegendAlignment(SwingConstants.SOUTH_EAST);
                                break;
                        }
                        traceChart.setShowLegend(legendCombo.getSelectedIndex() != 0);
                        validate();
                        repaint();
                    }
                }
        );

        colourByCombo.addItemListener(
                new java.awt.event.ItemListener() {
                    public void itemStateChanged(java.awt.event.ItemEvent ev) {
                        colourBy = colourByCombo.getSelectedIndex();
                        setupTraces();
                    }
                }
        );

    }

    private TraceList[] traceLists = null;
    private java.util.List<String> traceNames = null;

    public void setTraces(TraceList[] traceLists, java.util.List<String> traceNames) {
        this.traceLists = traceLists;
        this.traceNames = traceNames;
        setupTraces();
    }


    private void setupTraces() {
        traceChart.removeAllPlots();

        if (traceLists == null || traceNames == null || traceNames.size() == 0) {
            chartPanel.setXAxisTitle("");
            chartPanel.setYAxisTitle("");
            messageLabel.setText("No traces selected");
            add(messageLabel, BorderLayout.NORTH);
            return;
        }

        remove(messageLabel);

        Class iniTraceType = null;
        int i = 0;
//        Map<Integer, String> categoryDataMap = new HashMap<Integer, String>();
        int numOfBarsInt = 0;
        int numOfBarsCat = 0;
        for (TraceList tl : traceLists) {
            for (String traceName : traceNames) {
                int traceIndex = tl.getTraceIndex(traceName);
                Trace trace = tl.getTrace(traceIndex);
                if (trace != null) {
                    if (trace.getTraceType() == Integer.class) {
                        numOfBarsInt++;
                    } else if (trace.getTraceType() == String.class) {
                        numOfBarsCat++;
                    }
                }
            }
        }
        int barIntId = 0; // start from 0
        int barCatId = 0;
        for (TraceList tl : traceLists) {
            int n = tl.getStateCount();

            for (String traceName : traceNames) {
                int traceIndex = tl.getTraceIndex(traceName);
                Trace trace = tl.getTrace(traceIndex);
                TraceCorrelation td = tl.getCorrelationStatistics(traceIndex);
                FrequencyPlot plot = null;

                if (trace != null) {
                    if (iniTraceType == null) iniTraceType = trace.getTraceType();
                    if (iniTraceType == trace.getTraceType()) {
                        Map<Integer, String> categoryDataMap = new HashMap<Integer, String>();
                        if (trace.getTraceType() == Double.class) {
                            Double values[] = new Double[tl.getStateCount()];
                            tl.getValues(traceIndex, values);
                            boolean[] selected = new boolean[tl.getStateCount()];
                            tl.getSelected(traceIndex, selected);

                            if (td != null) {
                                plot = new NumericalDensityPlot(Trace.arrayConvert(values, selected), minimumBins, td);
                            } else {
                                plot = new NumericalDensityPlot(Trace.arrayConvert(values), minimumBins, td);
                            }

                            traceChart.setXAxis(false, new HashMap<Integer, String>());// make HashMap empty
                            chartPanel.setYAxisTitle("Density");

                            relativeDensityCheckBox.setVisible(true);
                            labelBins.setVisible(true);
                            binsCombo.setVisible(true);

                        } else if (trace.getTraceType() == Integer.class) {
                            Integer values[] = new Integer[tl.getStateCount()];
                            tl.getValues(traceIndex, values);
                            boolean[] selected = new boolean[tl.getStateCount()];
                            tl.getSelected(traceIndex, selected);

                            if (td != null) {
                                plot = new CategoryDensityPlot(Trace.arrayConvert(values, selected), -1, td, numOfBarsInt, barIntId);
                            } else {
                                plot = new CategoryDensityPlot(Trace.arrayConvert(values), -1, td, numOfBarsInt, barIntId);
                            }

                            barIntId++;
                            traceChart.setXAxis(true, new HashMap<Integer, String>());
                            chartPanel.setYAxisTitle("Probability");

                            relativeDensityCheckBox.setVisible(false);
                            labelBins.setVisible(false);
                            binsCombo.setVisible(false);

                        } else if (trace.getTraceType() == String.class) {
                            String initValues[] = new String[tl.getStateCount()];
                            tl.getValues(traceIndex, initValues);
                            boolean[] selected = new boolean[tl.getStateCount()];
                            tl.getSelected(traceIndex, selected);

                            String[] values;
                            if (td != null) {
                                values = Trace.arrayConvert(initValues, selected);
                            } else {
                                values = Trace.arrayConvert(initValues);
                            }

                            int[] intData = new int[values.length];
                            for (int v = 0; v < values.length; v++) {
                                intData[v] = td.credSet.getIndex(values[v]);
                                categoryDataMap.put(intData[v], values[v]);
                            }

                            plot = new CategoryDensityPlot(intData, -1, td, numOfBarsCat, barCatId);
                            barCatId++;
                            traceChart.setXAxis(false, categoryDataMap);
                            chartPanel.setYAxisTitle("Probability");

                            relativeDensityCheckBox.setVisible(false);
                            labelBins.setVisible(false);
                            binsCombo.setVisible(false);

                        } else {
                            throw new RuntimeException("Trace type is not recognized: " + trace.getTraceType());
                        }

                        String name = tl.getTraceName(traceIndex);
                        if (traceLists.length > 1) {
                            name = tl.getName() + " - " + name;
                        }

                        plot.setName(name);
                        if (tl instanceof CombinedTraces) {
                            plot.setLineStyle(new BasicStroke(2.0f), paints[i]);
                        } else {
                            plot.setLineStyle(new BasicStroke(1.0f), paints[i]);
                        }

                        traceChart.addPlot(plot);

                        if (colourBy == COLOUR_BY_TRACE || colourBy == COLOUR_BY_ALL) {
                            i++;
                        }
                        if (i == paints.length) i = 0;
                    } else {
//                         JOptionPane.showMessageDialog(frame,
//                                 "Selected traces contain different trace types\rso that the plot is displayed improperly.",
//                                 "Incompatible trace type",
//                                 JOptionPane.ERROR_MESSAGE);
                    }
                }
            }
            if (colourBy == COLOUR_BY_FILE) {
                i++;
            } else if (colourBy == COLOUR_BY_TRACE) {
                i = 0;
            }
            if (i == paints.length) i = 0;
        }

        if (traceLists.length == 1) {
            chartPanel.setXAxisTitle(traceLists[0].getName());
        } else if (traceNames.size() == 1) {
            chartPanel.setXAxisTitle(traceNames.get(0));
        } else {
            chartPanel.setXAxisTitle("Multiple Traces");
        }

        validate();
        repaint();
    }


    public JComponent getExportableComponent() {
        return chartPanel;
    }

    public String toString() {
        if (traceChart.getPlotCount() == 0) {
            return "no plot available";
        }

        StringBuffer buffer = new StringBuffer();

        Plot plot = traceChart.getPlot(0);
        Variate xData = plot.getXData();

        buffer.append(chartPanel.getXAxisTitle());
        for (int i = 0; i < traceChart.getPlotCount(); i++) {
            plot = traceChart.getPlot(i);
            buffer.append("\t");
            buffer.append(plot.getName());
        }
        buffer.append("\n");

        for (int i = 0; i < xData.getCount(); i++) {
            buffer.append(String.valueOf(xData.get(i)));
            for (int j = 0; j < traceChart.getPlotCount(); j++) {
                plot = traceChart.getPlot(j);
                Variate yData = plot.getYData();
                buffer.append("\t");
                buffer.append(String.valueOf(yData.get(i)));
            }
            buffer.append("\n");
        }

        return buffer.toString();
    }

}

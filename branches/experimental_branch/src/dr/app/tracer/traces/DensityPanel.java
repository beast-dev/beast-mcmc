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

import dr.app.gui.chart.*;
import dr.inference.trace.Trace;
import dr.inference.trace.TraceCorrelation;
import dr.inference.trace.TraceFactory;
import dr.inference.trace.TraceList;
import dr.stats.Variate;
import jam.framework.Exportable;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * A panel that displays density plots of traces
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: DensityPanel.java,v 1.3 2006/11/29 09:54:30 rambaut Exp $
 */
public class DensityPanel extends JPanel implements Exportable {
    private static final int DEFAULT_KDE_BINS = 5000;

    public static enum ColourByOptions {
        COLOUR_BY_TRACE,
        COLOUR_BY_FILE,
        COLOUR_BY_ALL
    };

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

    private class Settings {
        ChartSetupDialog chartSetupDialog = null;
        KDESetupDialog kdeSetupDialog = null;
        int minimumBins = 100;
        boolean showKDE = true;
        boolean showHistogram = false;
        boolean drawSolid = true;
        boolean relativeDensity = false;
        int barCount = 0;
        int legendAlignment = 0;
        ColourByOptions colourBy = ColourByOptions.COLOUR_BY_TRACE;
    }

    private Settings currentSettings = new Settings();
    private Map<String, Settings> settingsMap = new HashMap<String, Settings>();

    private JChart densityChart = new JChart(new LinearAxis(Axis.AT_MAJOR_TICK, Axis.AT_MAJOR_TICK), new LinearAxis());

    // as far as I can see DiscreteJChart is superfluous (discrete stats use the CategoryDensityPlot):
//    protected DiscreteJChart densityChart = new DiscreteJChart(new LinearAxis(Axis.AT_MAJOR_TICK_PLUS, Axis.AT_MAJOR_TICK), new LinearAxis());

    protected JChartPanel chartPanel = new JChartPanel(densityChart, null, "", "");

    protected JLabel labelBins;
    protected JComboBox binsCombo = new JComboBox(
            new Integer[]{10, 20, 50, 100, 200, 500, 1000});

    private JComboBox displayCombo = new JComboBox(
            new String[]{"KDE", "Histogram", "Both"}
    );

//    private JCheckBox kdeCheckBox = new JCheckBox("KDE");
//    private JButton kdeSetupButton = new JButton("Settings...");

    protected JCheckBox relativeDensityCheckBox = new JCheckBox("Relative density");
    private JCheckBox solidCheckBox = new JCheckBox("Fill plot");
    private JComboBox legendCombo = new JComboBox(
            new String[]{"None", "Top-Left", "Top", "Top-Right", "Left",
                    "Right", "Bottom-Left", "Bottom", "Bottom-Right"}
    );
    private JComboBox colourByCombo = new JComboBox(
            new String[]{"Trace", "Trace File", "All"}
    );
    private JButton chartSetupButton = new JButton("Axes...");
    private JLabel messageLabel = new JLabel("No data loaded");

    private TraceFactory.TraceType traceType = null;

    private final JFrame frame;

    /**
     * Creates new FrequencyPanel
     */
    public DensityPanel(final JFrame frame) {
        this.frame = frame;

        setOpaque(false);

        setMinimumSize(new Dimension(300, 150));
        setLayout(new BorderLayout());

        JToolBar toolBar = setupToolBar(frame);

        add(messageLabel, BorderLayout.NORTH);
        add(toolBar, BorderLayout.SOUTH);
        add(chartPanel, BorderLayout.CENTER);
    }

    protected JToolBar setupToolBar(final JFrame frame) {
        JToolBar toolBar = new JToolBar();
        toolBar.setOpaque(false);
        toolBar.setLayout(new FlowLayout(FlowLayout.LEFT));
        toolBar.setFloatable(false);

        chartSetupButton.putClientProperty(
                "Quaqua.Button.style", "placard"
        );
        chartSetupButton.setFont(UIManager.getFont("SmallSystemFont"));
        toolBar.add(chartSetupButton);

        JLabel label = new JLabel("Display:");
        label.setFont(UIManager.getFont("SmallSystemFont"));
        label.setLabelFor(displayCombo);
        toolBar.add(label);
        displayCombo.setFont(UIManager.getFont("SmallSystemFont"));
        displayCombo.setOpaque(false);
        toolBar.add(displayCombo);

        binsCombo.setFont(UIManager.getFont("SmallSystemFont"));
        binsCombo.setOpaque(false);
        binsCombo.setSelectedItem(currentSettings.minimumBins);
        labelBins = new JLabel("Bins:");
        labelBins.setFont(UIManager.getFont("SmallSystemFont"));
        labelBins.setLabelFor(binsCombo);
        toolBar.add(labelBins);
        toolBar.add(binsCombo);

        // KDE's don' do this at present so just taking up space on the toolbar...
//        relativeDensityCheckBox.setOpaque(false);
//        relativeDensityCheckBox.setFont(UIManager.getFont("SmallSystemFont"));
//        toolBar.add(relativeDensityCheckBox);

        // Probably don't need this as an option - takes up space and
        // solid (translucent) plots look cool...
//		solidCheckBox.setOpaque(false);
//		solidCheckBox.setFont(UIManager.getFont("SmallSystemFont"));
//		solidCheckBox.setSelected(true);
//		toolBar.add(solidCheckBox);

        toolBar.add(new JToolBar.Separator(new Dimension(8, 8)));
        label = new JLabel("Legend:");
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

//        kdeCheckBox.setFont(UIManager.getFont("SmallSystemFont"));
//        toolBar.add(kdeCheckBox);
//
//        kdeSetupButton.putClientProperty(
//                "Quaqua.Button.style", "placard"
//        );
//        kdeSetupButton.setFont(UIManager.getFont("SmallSystemFont"));
//        toolBar.add(kdeSetupButton);
//
//        kdeSetupButton.setEnabled(kdeCheckBox.isSelected());

        chartSetupButton.addActionListener(
                new java.awt.event.ActionListener() {
                    public void actionPerformed(ActionEvent actionEvent) {
                        if (currentSettings.chartSetupDialog == null) {
                            currentSettings.chartSetupDialog = new ChartSetupDialog(frame, true, false,
                                    Axis.AT_MAJOR_TICK, Axis.AT_MAJOR_TICK, Axis.AT_ZERO, Axis.AT_MAJOR_TICK);
                        }

                        currentSettings.chartSetupDialog.showDialog(densityChart);
                        validate();
                        repaint();
                    }
                }
        );

        displayCombo.addItemListener(
                new java.awt.event.ItemListener() {
                    public void itemStateChanged(java.awt.event.ItemEvent ev) {
                        currentSettings.showHistogram = displayCombo.getSelectedIndex() >= 1;
                        currentSettings.showKDE = displayCombo.getSelectedIndex() != 1;

                        binsCombo.setEnabled(currentSettings.showHistogram);
                        setupTraces();
                    }
                }
        );

        binsCombo.addItemListener(
                new java.awt.event.ItemListener() {
                    public void itemStateChanged(java.awt.event.ItemEvent ev) {
                        currentSettings.minimumBins = (Integer) binsCombo.getSelectedItem();
                        setupTraces();
                    }
                }
        );

        relativeDensityCheckBox.addItemListener(
                new java.awt.event.ItemListener() {
                    public void itemStateChanged(java.awt.event.ItemEvent ev) {
                        currentSettings.relativeDensity = relativeDensityCheckBox.isSelected();
                        setupTraces();
                    }
                }
        );

        solidCheckBox.addItemListener(
                new java.awt.event.ItemListener() {
                    public void itemStateChanged(java.awt.event.ItemEvent ev) {
                        currentSettings.drawSolid = solidCheckBox.isSelected();
                        setupTraces();
                    }
                }
        );

        legendCombo.addItemListener(
                new java.awt.event.ItemListener() {
                    public void itemStateChanged(java.awt.event.ItemEvent ev) {
                        currentSettings.legendAlignment = legendCombo.getSelectedIndex();
                        setupTraces();
                    }
                }
        );

        colourByCombo.addItemListener(
                new java.awt.event.ItemListener() {
                    public void itemStateChanged(java.awt.event.ItemEvent ev) {
                        currentSettings.colourBy = ColourByOptions.values()[colourByCombo.getSelectedIndex()];
                        setupTraces();
                    }
                }
        );

//        kdeCheckBox.addItemListener(
//                new java.awt.event.ItemListener() {
//                    public void itemStateChanged(java.awt.event.ItemEvent ev) {
//                        currentSettings.showKDE = kdeCheckBox.isSelected();
//                        kdeSetupButton.setEnabled(currentSettings.showKDE);
//                        setupTraces();
//                    }
//                }
//        );
//        kdeSetupButton.addActionListener(
//                new java.awt.event.ActionListener() {
//                    public void actionPerformed(ActionEvent actionEvent) {
//                        if (currentSettings.kdeSetupDialog == null) {
//                            currentSettings.kdeSetupDialog = new KDESetupDialog(frame, true, false,
//                                    Axis.AT_MAJOR_TICK, Axis.AT_MAJOR_TICK, Axis.AT_ZERO, Axis.AT_MAJOR_TICK);
//                        }
//
//                        currentSettings.kdeSetupDialog.showDialog(densityChart);
//                        setupTraces();
//                    }
//                }
//        );

        return toolBar;
    }

    private TraceList[] traceLists = null;
    private java.util.List<String> traceNames = null;

    public void setTraces(TraceList[] traceLists, java.util.List<String> traceNames) {
        this.traceLists = traceLists;
        this.traceNames = traceNames;

        if (traceNames.size() > 0) {
            // find the first settings for the one of the selected traces...
            Settings settings = null;

            for (String name : traceNames) {
                settings = settingsMap.get(name);
                if (settings != null) {
                    break;
                }
            }
            if (settings == null) {
                // if none of the traces have settings yet, create and store one for the
                // first selected trace
                settings = new Settings();
                settingsMap.put(traceNames.get(0), settings);
            }
            currentSettings = settings;
        }

        displayCombo.setSelectedIndex(currentSettings.showHistogram && currentSettings.showKDE ? 2 : (currentSettings.showKDE ? 0 : 1));
        binsCombo.setEnabled(currentSettings.showHistogram);

        binsCombo.setSelectedItem(currentSettings.minimumBins);
        relativeDensityCheckBox.setSelected(currentSettings.relativeDensity);
        legendCombo.setSelectedIndex(currentSettings.legendAlignment);
        colourByCombo.setSelectedIndex(currentSettings.colourBy.ordinal());
//        kdeCheckBox.setSelected(currentSettings.showKDE);
//        kdeSetupButton.setEnabled(currentSettings.showKDE);

        traceType = null;

//        barCount = 0;
        for (TraceList tl : traceLists) {
            for (String traceName : traceNames) {
                int traceIndex = tl.getTraceIndex(traceName);
                Trace trace = tl.getTrace(traceIndex);
                if (trace != null) {
                    if (traceType == null) {
                        traceType = trace.getTraceType();
                    }
                    if (trace.getTraceType() != traceType) {
                        densityChart.removeAllPlots();

                        chartPanel.setXAxisTitle("");
                        chartPanel.setYAxisTitle("");
                        messageLabel.setText("Unable to display a mixture statistics types.");
                        return;
                    }
                }
            }
        }

        // only enable controls relevent to continuous densities...
        displayCombo.setEnabled(traceType == TraceFactory.TraceType.DOUBLE);
        relativeDensityCheckBox.setEnabled(traceType == TraceFactory.TraceType.DOUBLE);
        labelBins.setEnabled(traceType == TraceFactory.TraceType.DOUBLE);
        binsCombo.setEnabled(traceType == TraceFactory.TraceType.DOUBLE);
//        kdeCheckBox.setEnabled(traceType == TraceFactory.TraceType.DOUBLE);
//        kdeSetupButton.setEnabled(traceType == TraceFactory.TraceType.DOUBLE);

        setupTraces();
    }

    protected Plot setupDensityPlot(TraceList tl, int traceIndex, TraceCorrelation td) {
        List values = tl.getValues(traceIndex);
        NumericalDensityPlot plot = new NumericalDensityPlot(values, currentSettings.minimumBins, td);
        return plot;
    }

    protected Plot setupKDEPlot(TraceList tl, int traceIndex, TraceCorrelation td) {
        List values = tl.getValues(traceIndex);
        Plot plot = new KDENumericalDensityPlot(values, DEFAULT_KDE_BINS, td);
        return plot;
    }

    protected Plot setupIntegerPlot(TraceList tl, int traceIndex, TraceCorrelation td, int barCount, int barId) {
        List values = tl.getValues(traceIndex);
        CategoryDensityPlot plot = new CategoryDensityPlot(values, -1, td, barCount, barId);
        return plot;
    }

    protected Plot setupCategoryPlot(TraceList tl, int traceIndex, TraceCorrelation td, Map<Integer, String> categoryDataMap, int barCount, int barId) {
        List values = tl.getValues(traceIndex);

        List<Double> intData = new ArrayList<Double>();
        for (int v = 0; v < values.size(); v++) {
            int index = td.getIndex(values.get(v).toString());
            intData.add(v, (double) index);
            categoryDataMap.put(index, values.get(v).toString());
        }

        CategoryDensityPlot plot = new CategoryDensityPlot(intData, -1, td, barCount, barId);

        return plot;
    }


    private void setupTraces() {
        densityChart.removeAllPlots();

        if (traceLists == null || traceNames == null || traceNames.size() == 0) {
            chartPanel.setXAxisTitle("");
            chartPanel.setYAxisTitle("");
            messageLabel.setText("No traces selected");
            add(messageLabel, BorderLayout.NORTH);
            return;
        }

        remove(messageLabel);

        int barId = 0;
        int i = 0;
        for (TraceList tl : traceLists) {
            int n = tl.getStateCount();

            for (String traceName : traceNames) {
                int traceIndex = tl.getTraceIndex(traceName);
                Trace trace = tl.getTrace(traceIndex);
                TraceCorrelation td = tl.getCorrelationStatistics(traceIndex);
                Plot plot = null;

                if (trace != null) {
                    String name = tl.getTraceName(traceIndex);
                    if (traceLists.length > 1) {
                        name = tl.getName() + " - " + name;
                    }

                    Map<Integer, String> categoryDataMap = new HashMap<Integer, String>();
                    if (trace.getTraceType() == TraceFactory.TraceType.DOUBLE) {
                        if (currentSettings.showHistogram) {
                            plot = setupDensityPlot(tl, traceIndex, td);
                            ((NumericalDensityPlot)plot).setRelativeDensity(currentSettings.relativeDensity);
                            ((NumericalDensityPlot)plot).setPointsOnly(currentSettings.showKDE);
                        } else {
                            plot = null;
                        }

                        if (currentSettings.showKDE) {
                            if (plot != null) {
                                ((NumericalDensityPlot)plot).setSolid(false);
                                ((NumericalDensityPlot)plot).setLineStroke(null);
                                ((NumericalDensityPlot)plot).setMarkStyle(Plot.POINT_MARK, 1, new BasicStroke(0.5f),
                                        Color.black, Color.black);
                            }

                            Plot plot2 = setupKDEPlot(tl, traceIndex, td);
                            plot2.setName(name + " KDE");
                            if (tl instanceof CombinedTraces) {
                                plot2.setLineStyle(new BasicStroke(2.0f), paints[i]);
                            } else {
                                plot2.setLineStyle(new BasicStroke(1.0f), paints[i]);
                            }
                            densityChart.addPlot(plot2);
                        }

                    } else if (trace.getTraceType() == TraceFactory.TraceType.INTEGER) {

                        plot = setupIntegerPlot(tl, traceIndex, td, currentSettings.barCount, barId);
                        barId++;

                    } else if (trace.getTraceType() == TraceFactory.TraceType.STRING) {

                        plot = setupCategoryPlot(tl, traceIndex, td, categoryDataMap, currentSettings.barCount, barId);
                        barId++;

                    } else {
                        throw new RuntimeException("Trace type is not recognized: " + trace.getTraceType());
                    }

                    if (plot != null) {
                        plot.setName(name);
                        if (tl instanceof CombinedTraces) {
                            plot.setLineStyle(new BasicStroke(2.0f), paints[i]);
                        } else {
                            plot.setLineStyle(new BasicStroke(1.0f), paints[i]);
                        }

                        densityChart.addPlot(plot);
                    }
                    if (currentSettings.colourBy == ColourByOptions.COLOUR_BY_TRACE || currentSettings.colourBy == ColourByOptions.COLOUR_BY_ALL) {
                        i++;
                    }
                    if (i == paints.length) i = 0;
                }
            }
            if (currentSettings.colourBy == ColourByOptions.COLOUR_BY_FILE) {
                i++;
            } else if (currentSettings.colourBy == ColourByOptions.COLOUR_BY_TRACE) {
                i = 0;
            }
            if (i == paints.length) i = 0;
        }

        switch (currentSettings.legendAlignment) {
            case 0:
                break;
            case 1:
                densityChart.setLegendAlignment(SwingConstants.NORTH_WEST);
                break;
            case 2:
                densityChart.setLegendAlignment(SwingConstants.NORTH);
                break;
            case 3:
                densityChart.setLegendAlignment(SwingConstants.NORTH_EAST);
                break;
            case 4:
                densityChart.setLegendAlignment(SwingConstants.WEST);
                break;
            case 5:
                densityChart.setLegendAlignment(SwingConstants.EAST);
                break;
            case 6:
                densityChart.setLegendAlignment(SwingConstants.SOUTH_WEST);
                break;
            case 7:
                densityChart.setLegendAlignment(SwingConstants.SOUTH);
                break;
            case 8:
                densityChart.setLegendAlignment(SwingConstants.SOUTH_EAST);
                break;
        }
        densityChart.setShowLegend(currentSettings.legendAlignment != 0);

        if (currentSettings.chartSetupDialog != null) {
            currentSettings.chartSetupDialog.applySettings(densityChart);
        }

        if (traceLists.length == 1) {
            chartPanel.setXAxisTitle(traceLists[0].getName());
        } else if (traceNames.size() == 1) {
            chartPanel.setXAxisTitle(traceNames.get(0));
        } else {
            chartPanel.setXAxisTitle("Multiple Traces");
        }

        if (traceType == TraceFactory.TraceType.DOUBLE) {
            chartPanel.setYAxisTitle("Density");
//            densityChart.setXAxis(false, new HashMap<Integer, String>());// make HashMap empty
        } else {
            chartPanel.setYAxisTitle("Probability");
        }

        validate();
        repaint();
    }


    public JComponent getExportableComponent() {
        return chartPanel;
    }

    public String toString() {
        if (densityChart.getPlotCount() == 0) {
            return "no plot available";
        }

        StringBuffer buffer = new StringBuffer();

        Plot plot = densityChart.getPlot(0);
        Variate xData = plot.getXData();

        buffer.append(chartPanel.getXAxisTitle());
        for (int i = 0; i < densityChart.getPlotCount(); i++) {
            plot = densityChart.getPlot(i);
            buffer.append("\t");
            buffer.append(plot.getName());
        }
        buffer.append("\n");

        for (int i = 0; i < xData.getCount(); i++) {
            buffer.append(String.valueOf(xData.get(i)));
            for (int j = 0; j < densityChart.getPlotCount(); j++) {
                plot = densityChart.getPlot(j);
                Variate yData = plot.getYData();
                buffer.append("\t");
                buffer.append(String.valueOf(yData.get(i)));
            }
            buffer.append("\n");
        }

        return buffer.toString();
    }

}

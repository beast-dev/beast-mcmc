/*
 * CorrelationPanel.java
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
import dr.inference.trace.TraceDistribution;
import dr.inference.trace.TraceFactory;
import dr.inference.trace.TraceList;
import dr.stats.Variate;
import jam.framework.Exportable;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * A panel that displays correlation plots of 2 traces
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: CorrelationPanel.java,v 1.1.1.2 2006/04/25 23:00:09 rambaut Exp $
 */
public class JointDensityPanel extends JPanel implements Exportable {

    private ChartSetupDialog chartSetupDialog = null;

    private JIntervalsChart correlationChart = new JIntervalsChart(new LinearAxis(), new LinearAxis());
    private JChartPanel chartPanel = new JChartPanel(correlationChart, null, "", "");
    private TableScrollPane tableScrollPane = new TableScrollPane();

    private JComboBox cateTableProbTypeCombo = new JComboBox(CateTableProbType.values());
    private JCheckBox defaultNumberFormatCheckBox = new JCheckBox("Use default number format");

    private JLabel messageLabel = new JLabel("No data loaded");

    private JCheckBox sampleCheckBox = new JCheckBox("Sample only");
    private JCheckBox pointsCheckBox = new JCheckBox("Draw as points");
    private JCheckBox translucencyCheckBox = new JCheckBox("Use translucency");

    private TraceList tl1 = null;
    private TraceList tl2 = null;
    private int traceIndex1 = -1;
    private int traceIndex2 = -1;

    private String name1;
    private String name2;

    public enum CateTableProbType {
        JOINT_PRO("Joint Probability"), COND_PRO_X("Conditional Prob (?|row)"),
        COND_PRO_Y("Conditional Prob (?|column)"), COUNT("Count");

        CateTableProbType(String name) {
            this.name = name;
        }

        public String toString() {
            return name;
        }

        private final String name;
    }

    /**
     * Creates new CorrelationPanel
     */
    public JointDensityPanel(final JFrame frame) {

        setOpaque(false);
        setMinimumSize(new Dimension(300, 150));
        setLayout(new BorderLayout());

//        add(messageLabel, BorderLayout.NORTH);
//        add(chartPanel, BorderLayout.CENTER);

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

        sampleCheckBox.setOpaque(false);
        sampleCheckBox.setFont(UIManager.getFont("SmallSystemFont"));
        sampleCheckBox.setSelected(true);
        toolBar.add(sampleCheckBox);

        pointsCheckBox.setOpaque(false);
        pointsCheckBox.setFont(UIManager.getFont("SmallSystemFont"));
        toolBar.add(pointsCheckBox);

        translucencyCheckBox.setOpaque(false);
        translucencyCheckBox.setFont(UIManager.getFont("SmallSystemFont"));
        toolBar.add(translucencyCheckBox);

        cateTableProbTypeCombo.setOpaque(false);
        cateTableProbTypeCombo.setFont(UIManager.getFont("SmallSystemFont"));
        toolBar.add(cateTableProbTypeCombo);

        defaultNumberFormatCheckBox.setOpaque(false);
        defaultNumberFormatCheckBox.setFont(UIManager.getFont("SmallSystemFont"));
        defaultNumberFormatCheckBox.setSelected(true);
        toolBar.add(defaultNumberFormatCheckBox);

        toolBar.add(new JToolBar.Separator(new Dimension(8, 8)));

        add(messageLabel, BorderLayout.NORTH);
        add(toolBar, BorderLayout.SOUTH);
        add(chartPanel, BorderLayout.CENTER);

        chartSetupButton.addActionListener(
                new java.awt.event.ActionListener() {
                    public void actionPerformed(ActionEvent actionEvent) {
                        if (chartSetupDialog == null) {
                            chartSetupDialog = new ChartSetupDialog(frame, true, true,
                                    Axis.AT_MAJOR_TICK, Axis.AT_MAJOR_TICK, Axis.AT_MAJOR_TICK, Axis.AT_MAJOR_TICK);
                        }

                        chartSetupDialog.showDialog(correlationChart);
                        validate();
                        repaint();
                    }
                }
        );

        ActionListener listener = new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent ev) {
                setupChartOrTable();
            }
        };
        sampleCheckBox.addActionListener(listener);
        pointsCheckBox.addActionListener(listener);
        translucencyCheckBox.addActionListener(listener);
        cateTableProbTypeCombo.addActionListener(listener);
        defaultNumberFormatCheckBox.addActionListener(listener);
    }

    public void setCombinedTraces() {
        chartPanel.setXAxisTitle("");
        chartPanel.setYAxisTitle("");
        messageLabel.setText("Can't show correlation of combined traces");
    }

    public void setTraces(TraceList[] traceLists, java.util.List<String> traceNames) {

//        correlationChart.removeAllPlots();

        if (traceLists != null && traceNames != null && traceLists.length == 2 && traceNames.size() == 1) {
            tl1 = traceLists[0];
            name1 = tl1.getName();
            tl2 = traceLists[1];
            name2 = tl2.getName();
            traceIndex1 = tl1.getTraceIndex(traceNames.get(0));
            traceIndex2 = tl2.getTraceIndex(traceNames.get(0));
            name1 = name1 + " - " + tl1.getTraceName(traceIndex1);
            name2 = name2 + " - " + tl2.getTraceName(traceIndex2);
        } else if (traceLists != null && traceNames != null && traceLists.length == 1 && traceNames.size() == 2) {
            tl1 = traceLists[0];
            tl2 = traceLists[0];
            traceIndex1 = tl1.getTraceIndex(traceNames.get(0));
            traceIndex2 = tl2.getTraceIndex(traceNames.get(1));
            name1 = tl1.getTraceName(traceIndex1);
            name2 = tl2.getTraceName(traceIndex2);
        } else {
            tl1 = null;
            tl2 = null;
        }

        setupChartOrTable();
    }

    private void setupChartOrTable() {
        correlationChart.removeAllIntervals();

        if (tl1 == null || tl2 == null) {
//            correlationChart.removeAllPlots();
            chartPanel.remove(tableScrollPane);

            chartPanel.setXAxisTitle("");
            chartPanel.setYAxisTitle("");
            messageLabel.setText("Select two statistics or traces from the table to view their correlation");
            return;
        }

        TraceDistribution td1 = tl1.getDistributionStatistics(traceIndex1);
        TraceDistribution td2 = tl2.getDistributionStatistics(traceIndex2);
        if (td1 == null || td2 == null) {
//            correlationChart.removeAllPlots();
            chartPanel.remove(tableScrollPane);

            chartPanel.setXAxisTitle("");
            chartPanel.setYAxisTitle("");
            messageLabel.setText("Waiting for analysis to complete");
            return;
        }

        messageLabel.setText("");

        if (td1.getTraceType() != TraceFactory.TraceType.DOUBLE && td2.getTraceType() != TraceFactory.TraceType.DOUBLE) {
            chartPanel.remove(correlationChart);
            chartPanel.add(tableScrollPane, "Table");

            sampleCheckBox.setVisible(false);
            pointsCheckBox.setVisible(false);
            translucencyCheckBox.setVisible(false);
            cateTableProbTypeCombo.setVisible(true);
            defaultNumberFormatCheckBox.setVisible(true);

            Object[] rowNames = td1.getRange().toArray();
            Object[] colNames = td2.getRange().toArray();
            double[][] data = categoricalPlot(td1, td2);

            tableScrollPane.setTable(rowNames, colNames, data, defaultNumberFormatCheckBox.isSelected());

        } else {
            chartPanel.remove(tableScrollPane);
            chartPanel.add(correlationChart, "Chart");
//            correlationChart.removeAllPlots();
            cateTableProbTypeCombo.setVisible(false);
            defaultNumberFormatCheckBox.setVisible(false);

            if (td1.getTraceType() == TraceFactory.TraceType.STRING) {
                mixedCategoricalPlot(td1, false); // isFirstTraceListNumerical

                sampleCheckBox.setVisible(false);
                pointsCheckBox.setVisible(false);
                translucencyCheckBox.setVisible(false);


            } else if (td2.getTraceType() == TraceFactory.TraceType.STRING) {
                mixedCategoricalPlot(td2, true); // isFirstTraceListNumerical

                sampleCheckBox.setVisible(false);
                pointsCheckBox.setVisible(false);
                translucencyCheckBox.setVisible(false);

                String swapName = name1;
                name1 = name2;
                name2 = swapName;

            } else {
                numericalPlot(td1, td2);

                sampleCheckBox.setVisible(true);
                pointsCheckBox.setVisible(true);
                translucencyCheckBox.setVisible(true);
            }
        }
        chartPanel.setXAxisTitle(name1);
        chartPanel.setYAxisTitle(name2);

        validate();
        repaint();
    }

    private void mixedCategoricalPlot(TraceDistribution td, boolean isFirstTraceListNumerical) {
        correlationChart.setXAxis(new DiscreteAxis(true, true));
        List<String> categoryValues = td.getRange();
        Map<String, TraceDistribution> categoryTdMap = new HashMap<String, TraceDistribution>();

        if (categoryValues == null || categoryValues.size() < 1) return;

        int maxCount = Math.max(tl1.getStateCount(), tl2.getStateCount());
        int minCount = Math.min(tl1.getStateCount(), tl2.getStateCount());

        int sampleSize = minCount;

        double samples1[] = new double[sampleSize];
        int k = 0;

        List values;
        if (isFirstTraceListNumerical) {
            values = tl1.getValues(traceIndex1);
        } else {
            values = tl2.getValues(traceIndex2);
        }

        for (int i = 0; i < sampleSize; i++) {
            samples1[i] = ((Number) values.get(k)).doubleValue();
            k += minCount / sampleSize;
        }


        String samples2[] = new String[sampleSize];
        k = 0;

        List values2;
        if (isFirstTraceListNumerical) {
            values2 = tl2.getValues(traceIndex2);
        } else {
            values2 = tl1.getValues(traceIndex1);
        }
        for (int i = 0; i < sampleSize; i++) {
            samples2[i] = values2.get(k).toString();
            k += minCount / sampleSize;
        }

        // separate samples into categoryTdMap
        ArrayList[] sepValues = new ArrayList[categoryValues.size()];
        for (int i = 0; i < categoryValues.size(); i++) {
            sepValues[i] = new ArrayList<Double>();
            for (int j = 0; j < samples2.length; j++) {
                if (categoryValues.get(i).equals(samples2[j])) {
                    sepValues[i].add(samples1[j]);
                }
            }

            TraceDistribution categoryTd = new TraceDistribution(sepValues[i], TraceFactory.TraceType.DOUBLE); // todo ?
            categoryTdMap.put(categoryValues.get(i), categoryTd);
        }

        for (String categoryValue : categoryValues) {
            TraceDistribution categoryTd = categoryTdMap.get(categoryValue);
            correlationChart.addIntervals(categoryValue, categoryTd.getMean(), categoryTd.getUpperHPD(), categoryTd.getLowerHPD(), false);
        }
    }

    private double[][] categoricalPlot(TraceDistribution td1, TraceDistribution td2) {
        List<String> rowNames = td1.getRange();
        List<String> colNames = td2.getRange();

        double[][] data = new double[rowNames.size()][colNames.size()];

        int maxCount = Math.max(tl1.getStateCount(), tl2.getStateCount());
        int minCount = Math.min(tl1.getStateCount(), tl2.getStateCount());

        int sampleSize = minCount;

        if (sampleSize <= 0) System.err.println("sampleSize cannot be 0. sampleSize = " + sampleSize);

        String samples1[] = new String[sampleSize];
        int k = 0;

        List values = tl1.getValues(traceIndex1);
        TraceFactory.TraceType type = tl1.getTrace(traceIndex1).getTraceType();
        for (int i = 0; i < sampleSize; i++) {
            if (type == TraceFactory.TraceType.INTEGER) { // as Integer is stored as Double in Trace
                samples1[i] = Integer.toString( ((Number) values.get(k)).intValue() );
            } else {
                samples1[i] = values.get(k).toString();
            }
            k += minCount / sampleSize; // = 1 for non-continous vs non-continous
        }

        String samples2[] = new String[sampleSize];
        k = 0;

        values = tl2.getValues(traceIndex2);
        type = tl2.getTrace(traceIndex2).getTraceType();
        for (int i = 0; i < sampleSize; i++) {
            if (type == TraceFactory.TraceType.INTEGER) { // as Integer is stored as Double in Trace
                samples2[i] = Integer.toString( ((Number) values.get(k)).intValue() );
            } else {
                samples2[i] = values.get(k).toString();
            }
            k += minCount / sampleSize;
        }

        // calculate count
        for (int i = 0; i < sampleSize; i++) {
            if (rowNames.contains(samples1[i]) && colNames.contains(samples2[i])) {
                data[rowNames.indexOf(samples1[i])][colNames.indexOf(samples2[i])] += 1;
            } else {
//                System.err.println("Not find row or column name. i = " + i);
            }
        }

        if (cateTableProbTypeCombo.getSelectedItem() == CateTableProbType.JOINT_PRO) {
            for (int r = 0; r < data.length; r++) {
                for (int c = 0; c < data[0].length; c++) {
                    data[r][c] = data[r][c] / sampleSize;
                }
            }
        } else if (cateTableProbTypeCombo.getSelectedItem() == CateTableProbType.COND_PRO_X) {
            for (int r = 0; r < data.length; r++) {
                double count = 0;
                for (int c = 0; c < data[0].length; c++) {
                    count = count + data[r][c];
                }
                for (int c = 0; c < data[0].length; c++) {
                    if (count != 0)
                        data[r][c] = data[r][c] / count;
                }
            }

        } else if (cateTableProbTypeCombo.getSelectedItem() == CateTableProbType.COND_PRO_Y) {
            for (int c = 0; c < data[0].length; c++) {
                double count = 0;
                for (int r = 0; r < data.length; r++) {
                    count = count + data[r][c];
                }
                for (int r = 0; r < data.length; r++) {
                    if (count != 0)
                        data[r][c] = data[r][c] / count;
                }
            }

        }
        // else COUNT

        return data;
    }

    private void numericalPlot(TraceDistribution td1, TraceDistribution td2) {
        int maxCount = Math.max(tl1.getStateCount(), tl2.getStateCount());
        int minCount = Math.min(tl1.getStateCount(), tl2.getStateCount());

        int sampleSize = minCount;

        if (sampleCheckBox.isSelected()) {
            if (td1.getESS() < td2.getESS()) {
                sampleSize = (int) td1.getESS();
            } else {
                sampleSize = (int) td2.getESS();
            }
            if (sampleSize < 20) {
                sampleSize = 20;
                messageLabel.setText("One of the traces has an ESS < 20 so a sample size of 20 will be used");
            }
            if (sampleSize > 500) {
                messageLabel.setText("This plot has been sampled down to 500 points");
                sampleSize = 500;
            }
        }

        int k = 0;
        if (td1.getTraceType() == TraceFactory.TraceType.INTEGER) {
            correlationChart.setXAxis(new DiscreteAxis(true, true));
        } else {
            correlationChart.setXAxis(new LinearAxis());
        }
        List values = tl1.getValues(traceIndex1);

        List<Double> samples1 = new ArrayList<Double>();
        for (int i = 0; i < sampleSize; i++) {
            samples1.add(i, ((Number) values.get(k)).doubleValue());
            k += minCount / sampleSize;
        }

        k = 0;
        if (td2.getTraceType() == TraceFactory.TraceType.INTEGER) {
            correlationChart.setYAxis(new DiscreteAxis(true, true));
        } else {
            correlationChart.setYAxis(new LinearAxis());
        }
        values = tl2.getValues(traceIndex2);

        List<Double> samples2 = new ArrayList<Double>();
        for (int i = 0; i < sampleSize; i++) {
            samples2.add(i, ((Number) values.get(k)).doubleValue());
            k += minCount / sampleSize;
        }

        ScatterPlot plot = new ScatterPlot(samples1, samples2);
        plot.setMarkStyle(pointsCheckBox.isSelected() ? Plot.POINT_MARK : Plot.CIRCLE_MARK, pointsCheckBox.isSelected() ? 1.0 : 3.0,
                new BasicStroke(2.0f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_MITER),
                new Color(16, 16, 64, translucencyCheckBox.isSelected() ? 32 : 255),
                new Color(16, 16, 64, translucencyCheckBox.isSelected() ? 32 : 255));
        correlationChart.addPlot(plot);
    }

//    private double[] removeNaN(double[] sample) {
//        List<Double> selectedValuesList = new ArrayList<Double>();
//
//        for (int i = 0; i < sample.length; i++) {
//            if (sample[i] != Double.NaN) {
//                selectedValuesList.add(sample[i]);
//            }
//        }
//
//        double[] dest = new double[selectedValuesList.size()];
//        for (int i = 0; i < dest.length; i++) {
//            dest[i] = selectedValuesList.get(i).doubleValue();
//        }
//
//        return dest;
//    }

    public JComponent getExportableComponent() {
        return chartPanel;
    }

    public String toString() {
        if (correlationChart.getPlotCount() == 0) {
            return "no plot available";
        }

        StringBuffer buffer = new StringBuffer();

        Plot plot = correlationChart.getPlot(0);
        Variate xData = plot.getXData();
        Variate yData = plot.getYData();

        buffer.append(chartPanel.getXAxisTitle());
        buffer.append("\t");
        buffer.append(chartPanel.getYAxisTitle());
        buffer.append("\n");

        for (int i = 0; i < xData.getCount(); i++) {
            buffer.append(String.valueOf(xData.get(i)));
            buffer.append("\t");
            buffer.append(String.valueOf(yData.get(i)));
            buffer.append("\n");
        }

        return buffer.toString();
    }
}

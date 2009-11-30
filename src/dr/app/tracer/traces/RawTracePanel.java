package dr.app.tracer.traces;

import dr.gui.chart.*;
import dr.inference.trace.TraceList;
import org.virion.jam.framework.Exportable;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;

/**
 * A panel that displays information about traces
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: RawTracePanel.java,v 1.2 2006/11/30 17:39:29 rambaut Exp $
 */
public class RawTracePanel extends JPanel implements Exportable {
    public static int COLOUR_BY_TRACE = 0;
    public static int COLOUR_BY_FILE = 1;
    public static int COLOUR_BY_ALL = 2;

    private static final Color[] paints = new Color[]{
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

    private JTraceChart traceChart = new JTraceChart(
            new LinearAxis(Axis.AT_ZERO, Axis.AT_DATA),
            new LinearAxis());
    private JChartPanel chartPanel = new JChartPanel(traceChart, null, "", "");

    private JCheckBox burninCheckBox = new JCheckBox("Show Burn-in");
    private JCheckBox sampleCheckBox = new JCheckBox("Sample only");
    private JCheckBox linePlotCheckBox = new JCheckBox("Draw line plot");
    private JComboBox legendCombo = new JComboBox(
            new String[]{"None", "Top-Left", "Top", "Top-Right", "Left",
                    "Right", "Bottom-Left", "Bottom", "Bottom-Right"}
    );
    private JComboBox colourByCombo = new JComboBox(
            new String[]{"Trace", "Trace File", "All"}
    );
    private JLabel messageLabel = new JLabel("No data loaded");

    private int colourBy = COLOUR_BY_TRACE;
    //private boolean sampleOnly = true;

    /**
     * Creates new RawTracePanel
     */
    public RawTracePanel(final JFrame frame) {

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

        burninCheckBox.setSelected(true);
        burninCheckBox.setFont(UIManager.getFont("SmallSystemFont"));
        burninCheckBox.setOpaque(false);
        toolBar.add(burninCheckBox);

        sampleCheckBox.setSelected(true);
        sampleCheckBox.setFont(UIManager.getFont("SmallSystemFont"));
        sampleCheckBox.setOpaque(false);
        toolBar.add(sampleCheckBox);

        toolBar.add(new JToolBar.Separator(new Dimension(8, 8)));
        linePlotCheckBox.setSelected(true);
        linePlotCheckBox.setFont(UIManager.getFont("SmallSystemFont"));
        linePlotCheckBox.setOpaque(false);
        toolBar.add(linePlotCheckBox);

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
                            chartSetupDialog = new ChartSetupDialog(frame, false, true,
                                    Axis.AT_ZERO, Axis.AT_MAJOR_TICK,
                                    Axis.AT_MAJOR_TICK, Axis.AT_MAJOR_TICK);
                        }

                        chartSetupDialog.showDialog(traceChart);
                        validate();
                        repaint();
                    }
                }
        );

        burninCheckBox.addActionListener(
                new java.awt.event.ActionListener() {
                    public void actionPerformed(java.awt.event.ActionEvent ev) {
                        setupTraces();
                    }
                }
        );

        sampleCheckBox.addActionListener(
                new java.awt.event.ActionListener() {
                    public void actionPerformed(java.awt.event.ActionEvent ev) {
                        traceChart.setUseSample(sampleCheckBox.isSelected());
                        validate();
                        repaint();
                    }
                }
        );

        linePlotCheckBox.addActionListener(
                new java.awt.event.ActionListener() {
                    public void actionPerformed(java.awt.event.ActionEvent ev) {
                        traceChart.setIsLinePlot(linePlotCheckBox.isSelected());
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
                        validate();
                        repaint();
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

        traceChart.removeAllTraces();


        if (traceLists == null || traceNames == null || traceNames.size() == 0) {
            chartPanel.setXAxisTitle("");
            chartPanel.setYAxisTitle("");
            messageLabel.setText("No traces selected");
            add(messageLabel, BorderLayout.NORTH);
            return;
        }

        remove(messageLabel);

        int i = 0;
        for (TraceList tl : traceLists) {
            int stateStart = tl.getBurnIn();
            int stateStep = tl.getStepSize();

            for (String traceName : traceNames) {
                double values[] = new double[tl.getStateCount()];
                int traceIndex = tl.getTraceIndex(traceName);
                tl.getValues(traceIndex, values);
                String name = tl.getTraceName(traceIndex);
                if (traceLists.length > 1) {
                    name = tl.getName() + " - " + name;
                }
                traceChart.addTrace(name, stateStart, stateStep, values, paints[i]);

                if (burninCheckBox.isSelected()) {
                    double burninValues[] = new double[tl.getBurninStateCount()];
                    tl.getBurninValues(traceIndex, burninValues);

                    Color burninColor = new Color(paints[i].getRed(), paints[i].getGreen(), paints[i].getBlue(), 96);
                    traceChart.addBurnin(name, stateStep, burninValues, burninColor, false);
                }
                
                if (colourBy == COLOUR_BY_TRACE || colourBy == COLOUR_BY_ALL) {
                    i++;
                }
                if (i == paints.length) i = 0;
            }
            if (colourBy == COLOUR_BY_FILE) {
                i++;
            } else if (colourBy == COLOUR_BY_TRACE) {
                i = 0;
            }
            if (i == paints.length) i = 0;
        }

        chartPanel.setXAxisTitle("State");
        if (traceLists.length == 1) {
            chartPanel.setYAxisTitle(traceLists[0].getName());
        } else if (traceNames.size() == 1) {
            chartPanel.setYAxisTitle(traceNames.get(0));
        } else {
            chartPanel.setYAxisTitle("Multiple Traces");
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

        //Plot plot = traceChart.getPlot(0);

        double[][] traceStates = new double[traceChart.getPlotCount()][];
        double[][] traceValues = new double[traceChart.getPlotCount()][];
        int maxLength = 0;

        for (int i = 0; i < traceChart.getPlotCount(); i++) {
            Plot plot = traceChart.getPlot(i);
            if (i > 0) {
                buffer.append("\t");
            }
            buffer.append("state");
            buffer.append("\t");
            buffer.append(plot.getName());

            traceStates[i] = traceChart.getTraceStates(i);
            traceValues[i] = traceChart.getTraceValues(i);
            if (traceStates[i].length > maxLength) {
                maxLength = traceStates[i].length;
            }
        }
        buffer.append("\n");

        for (int i = 0; i < maxLength; i++) {
            if (traceStates[0].length > i) {
                buffer.append(Integer.toString((int)traceStates[0][i]));
                buffer.append("\t");
                buffer.append(String.valueOf(traceValues[0][i]));
            } else {
                buffer.append("\t");
            }
            for (int j = 1; j < traceStates.length; j++) {
                if (traceStates[j].length > i) {
                    buffer.append("\t");
                    buffer.append(Integer.toString((int)traceStates[j][i]));
                    buffer.append("\t");
                    buffer.append(String.valueOf(traceValues[j][i]));
                } else {
                    buffer.append("\t\t");
                }
            }
            buffer.append("\n");
        }

        return buffer.toString();
    }

}

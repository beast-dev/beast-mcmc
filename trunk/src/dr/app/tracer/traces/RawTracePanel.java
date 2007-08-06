package dr.app.tracer.traces;

import dr.gui.chart.Axis;
import dr.gui.chart.ChartSetupDialog;
import dr.gui.chart.JChartPanel;
import dr.gui.chart.LinearAxis;
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

    private JTraceChart traceChart = new JTraceChart(
            new LinearAxis(Axis.AT_ZERO, Axis.AT_DATA),
            new LinearAxis());
    private JChartPanel chartPanel = new JChartPanel(traceChart, null, "", "");

    private JCheckBox sampleCheckBox = new JCheckBox("Sample only");
    private JCheckBox linePlotCheckBox = new JCheckBox("Draw line plot");
    private JComboBox legendCombo = new JComboBox(
            new String[]{"None", "Top-Left", "Top", "Top-Right", "Left",
                    "Right", "Bottom-Left", "Bottom", "Bottom-Right"}
    );
    private JComboBox colourByCombo = new JComboBox(
            new String[]{"TraceOld", "TraceOld File", "All"}
    );
    private JLabel messageLabel = new JLabel("No data loaded");

    private int colourBy = COLOUR_BY_TRACE;
    private boolean sampleOnly = true;

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
                        setupTraces(traceLists, traceIndices);
                        validate();
                        repaint();
                    }
                }
        );

    }

    private TraceList[] traceLists = null;
    private int[] traceIndices = null;

    public void setTraces(TraceList[] traceLists, int[] traceIndices) {
        this.traceLists = traceLists;
        this.traceIndices = traceIndices;
        setupTraces(traceLists, traceIndices);
    }


    private void setupTraces(TraceList[] traceLists, int[] traceIndices) {

        traceChart.removeAllTraces();


        if (traceLists == null || traceIndices == null || traceIndices.length == 0) {
            chartPanel.setXAxisTitle("");
            chartPanel.setYAxisTitle("");
            messageLabel.setText("No traces selected");
            add(messageLabel, BorderLayout.NORTH);
            return;
        }

        remove(messageLabel);

        int i = 0;
        for (TraceList tl : traceLists) {
            int n = tl.getStateCount();

            int stateStart = tl.getBurnIn();
            int stateStep = tl.getStepSize();

            for (int j = 0; j < traceIndices.length; j++) {
                double values[] = new double[n];
                tl.getValues(traceIndices[j], values);
                String name = tl.getTraceName(traceIndices[j]);
                if (traceLists.length > 1) {
                    name = tl.getName() + " - " + name;
                }
                traceChart.addTrace(name, stateStart, stateStep, values, paints[i]);

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
        } else if (traceIndices.length == 1) {
            chartPanel.setYAxisTitle(traceLists[0].getTraceName(traceIndices[0]));
        } else {
            chartPanel.setYAxisTitle("Multiple Traces");
        }


        validate();
        repaint();
    }

    public JComponent getExportableComponent() {
        return chartPanel;
    }

    public void copyToClipboard() {
        java.awt.datatransfer.Clipboard clipboard =
                Toolkit.getDefaultToolkit().getSystemClipboard();

        java.awt.datatransfer.StringSelection selection =
                new java.awt.datatransfer.StringSelection(this.toString());

        clipboard.setContents(selection, selection);
    }

    public String toString() {
        return "";
    }

}

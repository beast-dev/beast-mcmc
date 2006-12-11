package dr.app.tracer.traces;

import dr.gui.chart.*;
import org.virion.jam.framework.Exportable;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionListener;
import java.awt.event.ActionEvent;

/**
 * A panel that displays correlation plots of 2 traces
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: CorrelationPanel.java,v 1.1.1.2 2006/04/25 23:00:09 rambaut Exp $
 */
public class CorrelationPanel extends JPanel implements Exportable {

    private ChartSetupDialog chartSetupDialog = null;

    private JChart correlationChart = new JChart(new LinearAxis(), new LinearAxis());
    private JChartPanel chartPanel = new JChartPanel(correlationChart, null, "", "");
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

    /** Creates new CorrelationPanel */
    public CorrelationPanel(final JFrame frame) {

        setOpaque(false);
        setMinimumSize(new Dimension(300,150));
        setLayout(new BorderLayout());

        add(messageLabel, BorderLayout.NORTH);
        add(chartPanel, BorderLayout.CENTER);

        JToolBar toolBar = new JToolBar();
        toolBar.setOpaque(false);
        toolBar.setLayout(new FlowLayout(FlowLayout.LEFT));
        toolBar.setFloatable(false);

        JButton chartSetupButton = new JButton("Setup Axes");
        toolBar.add(new JToolBar.Separator(new Dimension(8,8)));
        toolBar.add(chartSetupButton);

        sampleCheckBox.setOpaque(false);
        sampleCheckBox.setSelected(true);
        toolBar.add(sampleCheckBox);

        pointsCheckBox.setOpaque(false);
        toolBar.add(pointsCheckBox);

        translucencyCheckBox.setOpaque(false);
        toolBar.add(translucencyCheckBox);

        toolBar.add(new JToolBar.Separator(new Dimension(8,8)));

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
                setupChart();
            }
        };
        sampleCheckBox.addActionListener(listener);
        pointsCheckBox.addActionListener(listener);
        translucencyCheckBox.addActionListener(listener);
    }

    public void setCombinedTraces() {
        chartPanel.setXAxisTitle("");
        chartPanel.setYAxisTitle("");
        messageLabel.setText("Can't show correlation of combined traces");
    }

    public void setTraces(TraceList[] traceLists, int[] traceIndices) {

        correlationChart.removeAllPlots();

        if (traceLists != null && traceIndices != null && traceLists.length == 2 && traceIndices.length == 1) {
            tl1 = traceLists[0];
            name1 = tl1.getName();
            tl2 = traceLists[1];
            name2 = tl2.getName();
            traceIndex1 = traceIndices[0];
            traceIndex2 = traceIndices[0];
            name1 = name1 + " - " + tl1.getTraceName(traceIndex1);
            name2 = name2 + " - " + tl2.getTraceName(traceIndex2);
        } else if (traceLists != null && traceIndices != null && traceLists.length == 1 && traceIndices.length == 2) {
            tl1 = traceLists[0];
            tl2 = traceLists[0];
            traceIndex1 = traceIndices[0];
            traceIndex2 = traceIndices[1];
            name1 = tl1.getTraceName(traceIndex1);
            name2 = tl2.getTraceName(traceIndex2);
        } else {
            tl1 = null;
            tl2 = null;
        }

        setupChart();
    }

    private void setupChart() {

        correlationChart.removeAllPlots();

        if (tl1 == null || tl2 == null) {
            chartPanel.setXAxisTitle("");
            chartPanel.setYAxisTitle("");
            messageLabel.setText("Select two statistics or traces from the table to view their correlation");
            return;
        }

        TraceDistribution td1 = tl1.getDistributionStatistics(traceIndex1);
        TraceDistribution td2 = tl2.getDistributionStatistics(traceIndex2);
        if (td1 == null || td2 == null) {
            chartPanel.setXAxisTitle("");
            chartPanel.setYAxisTitle("");
            messageLabel.setText("Waiting for analysis to complete");
            return;
        }

        messageLabel.setText("");

        int count = Math.min(tl1.getStateCount(), tl2.getStateCount());


        int sampleSize = count;

        if (sampleCheckBox.isSelected()) {
            if (td1.getESS() < td2.getESS()) {
                sampleSize = (int)td1.getESS();
            } else {
                sampleSize = (int)td2.getESS();
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

        double values[] = new double[count];

        tl1.getValues(traceIndex1, values);

        double samples1[] = new double[sampleSize];
        int k = 0;
        for (int i = 0; i < sampleSize; i++) {
            samples1[i] = values[k];
            k += count / sampleSize;
        }

        tl2.getValues(traceIndex2, values);

        double samples2[] = new double[sampleSize];
        k = 0;
        for (int i = 0; i < sampleSize; i++) {
            samples2[i] = values[k];
            k += count / sampleSize;
        }

        ScatterPlot plot = new ScatterPlot(samples1, samples2);
        plot.setMarkStyle(pointsCheckBox.isSelected() ? Plot.POINT_MARK : Plot.CIRCLE_MARK, pointsCheckBox.isSelected() ? 1.0 : 3.0,
                new BasicStroke(2.0f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_MITER),
                new Color(16,16,64,translucencyCheckBox.isSelected() ? 32 : 255), new Color(16,16,64,translucencyCheckBox.isSelected() ? 32 : 255));
        correlationChart.addPlot(plot);

        chartPanel.setXAxisTitle(name1);
        chartPanel.setYAxisTitle(name2);

        validate();
        repaint();
    }

    public void copyToClipboard() {
        /*
          java.awt.datatransfer.Clipboard clipboard =
              Toolkit.getDefaultToolkit().getSystemClipboard();

          java.awt.datatransfer.StringSelection selection =
              new java.awt.datatransfer.StringSelection(statisticsTable.toString());

          clipboard.setContents(selection, selection);
      */
    }

    public JComponent getExportableComponent() {
        return chartPanel;
    }
}

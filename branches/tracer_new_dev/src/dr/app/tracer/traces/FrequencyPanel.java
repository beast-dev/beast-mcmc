package dr.app.tracer.traces;

import dr.gui.chart.*;
import dr.inference.trace.Trace;
import dr.inference.trace.TraceDistribution;
import dr.inference.trace.TraceList;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;

import org.virion.jam.framework.Exportable;


/**
 * A panel that displays frequency distributions of traces
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: FrequencyPanel.java,v 1.1.1.2 2006/04/25 23:00:09 rambaut Exp $
 */
public class FrequencyPanel extends JPanel implements Exportable {

    private ChartSetupDialog chartSetupDialog = null;

    private TraceList traceList = null;
    private String traceName = null;

    private int minimumBins = 50;

    private JComboBox binsCombo = new JComboBox(
            new Integer[]{10, 20, 50, 100, 200, 500, 1000});


    private JChart traceChart = new JChart(new LinearAxis(Axis.AT_MAJOR_TICK_PLUS, Axis.AT_MAJOR_TICK_PLUS), new LinearAxis());
    private JChartPanel chartPanel = new JChartPanel(traceChart, null, "", "Frequency");

    /**
     * Creates new FrequencyPanel
     */
    public FrequencyPanel(final JFrame frame) {
        setOpaque(false);

        setMinimumSize(new Dimension(300, 150));
        setLayout(new BorderLayout());
        add(chartPanel, BorderLayout.CENTER);

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

        binsCombo.setOpaque(false);
        binsCombo.setFont(UIManager.getFont("SmallSystemFont"));
        binsCombo.setSelectedItem(100);
        JLabel label = new JLabel("Bins:");
        label.setFont(UIManager.getFont("SmallSystemFont"));
        label.setLabelFor(binsCombo);
        toolBar.add(label);
        toolBar.add(binsCombo);

        add(toolBar, BorderLayout.SOUTH);

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
                        setupTrace();
                    }
                }
        );
    }

    public void setTrace(TraceList traceList, String traceName) {
        this.traceList = traceList;
        this.traceName = traceName;
        setupTrace();
    }


    private void setupTrace() {

        traceChart.removeAllPlots();

        if (traceList == null || traceName == null) {
            chartPanel.setXAxisTitle("");
            chartPanel.setYAxisTitle("");
            return;
        }

        double values[] = new double[traceList.getStateCount()];
        int traceIndex = traceList.getTraceIndex(traceName);

        traceList.getValues(traceIndex, Trace.arrayCopy(values));        

        FrequencyPlot plot = new FrequencyPlot(values, minimumBins);

        TraceDistribution td = traceList.getDistributionStatistics(traceIndex);
        if (td != null) {
            plot.setIntervals(td.getUpperHPD(), td.getLowerHPD());
        }

        traceChart.addPlot(plot);

        chartPanel.setXAxisTitle(traceList.getTraceName(traceIndex));
        chartPanel.setYAxisTitle("Frequency");
    }

    public JComponent getExportableComponent() {
        return chartPanel;
    }

}

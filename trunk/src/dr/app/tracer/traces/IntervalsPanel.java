package dr.app.tracer.traces;

import dr.gui.chart.JChartPanel;
import dr.gui.chart.LinearAxis;
import dr.inference.trace.TraceDistribution;
import dr.inference.trace.TraceList;

import javax.swing.*;
import java.awt.*;


/**
 * A panel that displays frequency distributions of traces
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: IntervalsPanel.java,v 1.1.1.2 2006/04/25 23:00:09 rambaut Exp $
 */
public class IntervalsPanel extends JPanel {

    private JIntervalsChart intervalsChart = new JIntervalsChart(new LinearAxis());
    private JChartPanel chartPanel = new JChartPanel(intervalsChart, null, "", "");

    /**
     * Creates new IntervalsPanel
     */
    public IntervalsPanel() {
        setOpaque(false);
        setMinimumSize(new Dimension(300, 150));
        setLayout(new BorderLayout());
        add(chartPanel, BorderLayout.CENTER);
    }

    public void setTraces(TraceList[] traceLists, int[] traceIndices) {

        intervalsChart.removeAllIntervals();

        if (traceLists == null || traceIndices == null) {
            chartPanel.setXAxisTitle("");
            chartPanel.setYAxisTitle("");
            return;
        }

        for (int i = 0; i < traceLists.length; i++) {
            for (int j = 0; j < traceIndices.length; j++) {
                TraceDistribution td = traceLists[i].getDistributionStatistics(traceIndices[j]);
                if (td != null) {
                    String name = traceLists[i].getTraceName(traceIndices[j]);
                    if (traceLists.length > 1) {
                        name = traceLists[i].getName() + " - " + name;
                    }
                    intervalsChart.addIntervals(name, td.getMean(), td.getUpperHPD(), td.getLowerHPD(), false);
                }
            }
        }

        chartPanel.setXAxisTitle("");
        if (traceLists.length == 1) {
            chartPanel.setYAxisTitle(traceLists[0].getName());
        } else if (traceIndices.length == 1) {
            chartPanel.setYAxisTitle(traceLists[0].getTraceName(traceIndices[0]));
        } else {
            chartPanel.setYAxisTitle("Multiple Traces");
        }
        add(chartPanel, BorderLayout.CENTER);

        validate();
        repaint();
    }

    public JComponent getExportableComponent() {
        return chartPanel;
    }

}

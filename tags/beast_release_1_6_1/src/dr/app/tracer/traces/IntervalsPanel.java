package dr.app.tracer.traces;

import dr.app.gui.chart.JChartPanel;
import dr.app.gui.chart.LinearAxis;
import dr.inference.trace.TraceDistribution;
import dr.inference.trace.TraceList;

import javax.swing.*;
import java.awt.*;

import jam.framework.Exportable;


/**
 * A panel that displays frequency distributions of traces
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: IntervalsPanel.java,v 1.1.1.2 2006/04/25 23:00:09 rambaut Exp $
 */
public class IntervalsPanel extends JPanel implements Exportable {

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

    public void setTraces(TraceList[] traceLists, java.util.List<String> traceNames) {

        intervalsChart.removeAllIntervals();

        if (traceLists == null || traceNames == null) {
            chartPanel.setXAxisTitle("");
            chartPanel.setYAxisTitle("");
            return;
        }

        for (TraceList traceList : traceLists) {
            for (String traceName : traceNames) {
                int index = traceList.getTraceIndex(traceName);
                TraceDistribution td = traceList.getDistributionStatistics(index);
                if (td != null) {
                    String name = "";
                    if (traceLists.length > 1) {
                        name = traceList.getName();
                        if (traceNames.size() > 1) {
                            name += ": ";
                        }
                    }
                    name += traceName;
                    intervalsChart.addIntervals(name, td.getMean(), td.getUpperHPD(), td.getLowerHPD(), false);
                }
            }
        }

        chartPanel.setXAxisTitle("");
        if (traceLists.length == 1) {
            chartPanel.setYAxisTitle(traceLists[0].getName());
        } else if (traceNames.size() == 1) {
            chartPanel.setYAxisTitle(traceNames.get(0));
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

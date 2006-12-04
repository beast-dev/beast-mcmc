package dr.app.tracer.traces;

import dr.gui.chart.*;

import javax.swing.*;
import java.awt.*;


/**
 * A panel that displays frequency distributions of traces
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: FrequencyPanel.java,v 1.1.1.2 2006/04/25 23:00:09 rambaut Exp $
 */
public class FrequencyPanel extends JPanel {

	TraceList traceList = null;
	int traceIndex = -1;

	private JChart traceChart = new JChart(new LinearAxis(Axis.AT_MAJOR_TICK_PLUS, Axis.AT_MAJOR_TICK_PLUS), new LinearAxis());
	private JChartPanel chartPanel = new JChartPanel(traceChart, null, "", "Frequency");

	/** Creates new FrequencyPanel */
	public FrequencyPanel() {
		setOpaque(false);

		setMinimumSize(new Dimension(300,150));
		setLayout(new BorderLayout());
		add(chartPanel, BorderLayout.CENTER);
	}

	public void setTrace(TraceList traceList, int traceIndex) {

		this.traceList = traceList;
		this.traceIndex = traceIndex;

		traceChart.removeAllPlots();

		if (traceList == null || traceIndex == -1) {
			chartPanel.setXAxisTitle("");
			chartPanel.setYAxisTitle("");
			return;
		}

		double values[] = new double[traceList.getStateCount()];
		traceList.getValues(traceIndex, values);
		FrequencyPlot plot = new FrequencyPlot(values, 50);

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

/*
 * FrequencyPanel.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.app.tracer;

import dr.gui.chart.*;

import javax.swing.*;
import java.awt.*;


/**
 * A panel that displays frequency distributions of traces
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: FrequencyPanel.java,v 1.14 2005/07/11 14:07:26 rambaut Exp $
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

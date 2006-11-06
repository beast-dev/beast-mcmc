/*
 * IntervalsPanel.java
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

import dr.gui.chart.JChartPanel;
import dr.gui.chart.LinearAxis;

import javax.swing.*;
import java.awt.*;


/**
 * A panel that displays frequency distributions of traces
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: IntervalsPanel.java,v 1.4 2005/07/11 14:07:26 rambaut Exp $
 */
public class IntervalsPanel extends JPanel {

	/**
	 * 
	 */
	private static final long serialVersionUID = 6835575541003219985L;
	private JIntervalsChart intervalsChart = new JIntervalsChart(new LinearAxis());
	private JChartPanel chartPanel = new JChartPanel(intervalsChart, null, "", "");
	
	/** Creates new IntervalsPanel */
	public IntervalsPanel() {		
		setOpaque(false);
		setMinimumSize(new Dimension(300,150));
		setLayout(new BorderLayout());
		add(chartPanel, BorderLayout.CENTER);
	}
	
	public void setTraces(TraceList traceList, int[] traceIndices) {

		intervalsChart.removeAllIntervals();

		if (traceList == null || traceIndices == null || traceIndices.length < 2) {
			chartPanel.setXAxisTitle("");
			chartPanel.setYAxisTitle("");
			return; 
		}

		for (int i = 0; i < traceIndices.length; i++) {
			TraceDistribution td = traceList.getDistributionStatistics(traceIndices[i]);
			if (td != null) {
				intervalsChart.addIntervals(traceList.getTraceName(traceIndices[i]),
							td.getMean(), td.getUpperHPD(), td.getLowerHPD(), false);
			}
		}
		
		chartPanel.setXAxisTitle("");
		chartPanel.setYAxisTitle(traceList.getName());
		add(chartPanel, BorderLayout.CENTER);
		
		validate();
		repaint();
	}
	
    public JComponent getExportableComponent() {
		return chartPanel;
	} 	
      
}

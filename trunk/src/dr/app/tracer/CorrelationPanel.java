/*
 * CorrelationPanel.java
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
import dr.util.Variate;
import org.virion.jam.framework.Exportable;

import javax.swing.*;
import java.awt.*;

/**
 * A panel that displays correlation plots of 2 traces
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: CorrelationPanel.java,v 1.9 2005/10/03 16:27:36 rambaut Exp $
 */
public class CorrelationPanel extends JPanel implements Exportable {

	private JChart correlationChart = new JChart(new LinearAxis(), new LinearAxis());
	private JChartPanel chartPanel = new JChartPanel(correlationChart, null, "", "");
	private JLabel messageLabel = new JLabel("No data loaded");
	
    private String xAxisTitle;
    private String yAxisTitle;

	/** Creates new CorrelationPanel */
	public CorrelationPanel() {

		setOpaque(false);
		setMinimumSize(new Dimension(300,150));
		setLayout(new BorderLayout());

		add(messageLabel, BorderLayout.NORTH);
		add(chartPanel, BorderLayout.CENTER);
	}
	
	public void setCombinedTraces() {
		chartPanel.setXAxisTitle("");
		chartPanel.setYAxisTitle("");
		messageLabel.setText("Can't show joint-marginal distribution of combined traces");
	}
		
	public void setTraces(TraceList traceList, int traceIndex1, int traceIndex2) {
		
		correlationChart.removeAllPlots();
		
		if (traceList == null || traceIndex1 == -1 || traceIndex2 == -1) {
            xAxisTitle = "";
            yAxisTitle = "";
			chartPanel.setXAxisTitle("");
			chartPanel.setYAxisTitle("");
			messageLabel.setText("Select two statistics from the table to view their joint-marginal distribution");
			return;
		}

		TraceCorrelation tc1 = traceList.getCorrelationStatistics(traceIndex1);
		TraceCorrelation tc2 = traceList.getCorrelationStatistics(traceIndex2);
		if (tc1 == null || tc2 == null) {
            xAxisTitle = "";
            yAxisTitle = "";
			chartPanel.setXAxisTitle("");
			chartPanel.setYAxisTitle("");
			messageLabel.setText("Waiting for analysis to complete");
			return;
		}

		messageLabel.setText("");
		
		int sampleSize;
		if (tc1.getESS() < tc2.getESS()) {
			sampleSize = (int)tc1.getESS();
		} else {
			sampleSize = (int)tc2.getESS();
		}
		if (sampleSize < 20) {
			sampleSize = 20;
			messageLabel.setText("One of the traces has an ESS < 20 so a sample size of 20 will be used");
		}
		if (sampleSize > 10000) {
			messageLabel.setText("This plot has been sampled down to 10,000 points");
			sampleSize = 10000;
		}
		
		int count = traceList.getStateCount();
		double values[] = new double[count];
		
		traceList.getValues(traceIndex1, values);
		
		double samples1[] = new double[sampleSize];
		int k = 0;
		for (int i = 0; i < sampleSize; i++) {
			samples1[i] = values[k];
			k += count / sampleSize;
		}
		
		traceList.getValues(traceIndex2, values);
		
		double samples2[] = new double[sampleSize];
		k = 0;
		for (int i = 0; i < sampleSize; i++) {
			samples2[i] = values[k];
			k += count / sampleSize;
		}
		
		ScatterPlot plot = new ScatterPlot(samples1, samples2);
		plot.setMarkStyle(Plot.POINT_MARK, 1.0, 
			new BasicStroke(2.0f, BasicStroke.CAP_ROUND, BasicStroke.JOIN_MITER), 
			Color.black, Color.black);
		correlationChart.addPlot(plot);
		
        xAxisTitle = traceList.getTraceName(traceIndex1);
        yAxisTitle = traceList.getTraceName(traceIndex2);

        chartPanel.setXAxisTitle(xAxisTitle);
		chartPanel.setYAxisTitle(yAxisTitle);
		
		validate();
		repaint();
	}
	
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

        buffer.append(xAxisTitle);
        buffer.append("\t");
        buffer.append(yAxisTitle);
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

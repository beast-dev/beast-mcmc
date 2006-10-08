/*
 * DensityPanel.java
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
 * A panel that displays density plots of traces
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: DensityPanel.java,v 1.20 2005/10/03 16:27:36 rambaut Exp $
 */
public class DensityPanel extends JPanel implements Exportable {

	private JChart traceChart = new JChart(new LinearAxis(Axis.AT_MAJOR_TICK_PLUS, Axis.AT_MAJOR_TICK_PLUS), new LinearAxis());
	private JChartPanel chartPanel = new JChartPanel(traceChart, null, "", "");
	
	private JCheckBox relativeDensityCheckBox = new JCheckBox("Relative density");
	private JComboBox legendCombo = new JComboBox(
		new String[] { "None", "Top-Left", "Top", "Top-Right", "Left", 
						"Right", "Bottom-Left", "Bottom", "Bottom-Right" }
	);
	private JCheckBox showIndividualCheckBox = new JCheckBox("Show individual traces");
	private JLabel messageLabel = new JLabel("No data loaded");
	
	private CombinedTraces combinedTraces = null;
	private int[] traceIndices = null;
    private String xAxisTitle = "";
	
	/** Creates new FrequencyPanel */
	public DensityPanel() {

		setOpaque(false);

		setMinimumSize(new Dimension(300,150));
		setLayout(new BorderLayout());

		showIndividualCheckBox.setSelected(true);
		showIndividualCheckBox.setEnabled(false);

		JToolBar toolBar = new JToolBar();
		toolBar.setOpaque(false);
		toolBar.setLayout(new FlowLayout(FlowLayout.LEFT));
		toolBar.setFloatable(false);
		relativeDensityCheckBox.setOpaque(false);
		toolBar.add(relativeDensityCheckBox);
		
		toolBar.add(new JToolBar.Separator(new Dimension(8,8)));
		JLabel label = new JLabel("Legend:");
		label.setFont(relativeDensityCheckBox.getFont());
		label.setLabelFor(legendCombo);
		toolBar.add(label);
		legendCombo.setFont(relativeDensityCheckBox.getFont());
		legendCombo.setOpaque(false);
		toolBar.add(legendCombo);

		toolBar.add(new JToolBar.Separator(new Dimension(8,8)));
		showIndividualCheckBox.setOpaque(false);
		toolBar.add(showIndividualCheckBox);

		add(messageLabel, BorderLayout.NORTH);
		add(toolBar, BorderLayout.SOUTH);
		add(chartPanel, BorderLayout.CENTER);
		
		relativeDensityCheckBox.addItemListener(
			new java.awt.event.ItemListener() {
				public void itemStateChanged(java.awt.event.ItemEvent ev) {
					for (int i = 0; i < traceChart.getPlotCount(); i++) {
						((DensityPlot)traceChart.getPlot(i)).setRelativeDensity(relativeDensityCheckBox.isSelected());
					}
					traceChart.recalibrate();
					validate();
					repaint();
				}
			}
		);
		
		legendCombo.addItemListener(
			new java.awt.event.ItemListener() {
				public void itemStateChanged(java.awt.event.ItemEvent ev) {
					switch (legendCombo.getSelectedIndex()) {
						case 0: break;
						case 1: traceChart.setLegendAlignment(SwingConstants.NORTH_WEST); break;
						case 2: traceChart.setLegendAlignment(SwingConstants.NORTH); break;
						case 3: traceChart.setLegendAlignment(SwingConstants.NORTH_EAST); break;
						case 4: traceChart.setLegendAlignment(SwingConstants.WEST); break;
						case 5: traceChart.setLegendAlignment(SwingConstants.EAST); break;
						case 6: traceChart.setLegendAlignment(SwingConstants.SOUTH_WEST); break;
						case 7: traceChart.setLegendAlignment(SwingConstants.SOUTH); break;
						case 8: traceChart.setLegendAlignment(SwingConstants.SOUTH_EAST); break;
					}
					traceChart.setShowLegend(legendCombo.getSelectedIndex() != 0);
					validate();
					repaint();
				}
			}
		);
		
		showIndividualCheckBox.addItemListener(
			new java.awt.event.ItemListener() {
				public void itemStateChanged(java.awt.event.ItemEvent ev) {
					if (combinedTraces != null) {
						setCombinedTraces(combinedTraces, traceIndices);
					}
				}
			}
		);
		
	}
	
	Paint[] paints = new Paint[] {
		Color.black,
		Color.blue,
		Color.red,
		Color.green,
		Color.magenta,
		Color.gray,
		Color.pink,
		Color.orange,
		Color.lightGray,
		Color.darkGray
	};
		
	public void setCombinedTraces(CombinedTraces combinedTraces, int[] traceIndices) {
		// This shows the combined trace overlayed on each of the individuals
	
		this.combinedTraces = combinedTraces;
		this.traceIndices = traceIndices;
		showIndividualCheckBox.setEnabled(true);
		
		traceChart.removeAllPlots();
		
		if (combinedTraces == null || traceIndices == null || traceIndices.length == 0) {
            xAxisTitle = "";
			chartPanel.setXAxisTitle("");
			chartPanel.setYAxisTitle("");
			messageLabel.setText("No traces selected");
			return; 
		}

		messageLabel.setText("");

		for (int i = 0; i < traceIndices.length; i++) {
		
			if (showIndividualCheckBox.isSelected() &&
				combinedTraces.getTraceListCount() > 1) {
				for (int j = 0; j < combinedTraces.getTraceListCount(); j++) {
				
					TraceList tl = combinedTraces.getTraceList(j);
					double values[] = new double[tl.getStateCount()];
				
					tl.getValues(traceIndices[i], values);
					DensityPlot plot = new DensityPlot(values, 50);
					plot.setName(null);
					plot.setLineStyle(new BasicStroke(1.0f), paints[i % paints.length]);
					traceChart.addPlot(plot);
				}
			}
						
			double values[] = new double[combinedTraces.getStateCount()];
			combinedTraces.getValues(traceIndices[i], values);
			DensityPlot plot = new DensityPlot(values, 50);
			plot.setName(combinedTraces.getTraceName(traceIndices[i]));
			plot.setLineStyle(new BasicStroke(2.0f), paints[i % paints.length]);
			traceChart.addPlot(plot);
		}
		
		if (traceIndices.length == 1) {
            xAxisTitle = combinedTraces.getTraceName(traceIndices[0]);
		} else {
			xAxisTitle = "Multiple parameters";
		}
        chartPanel.setXAxisTitle(xAxisTitle);
		chartPanel.setYAxisTitle("Density");

		validate();
		repaint();
	}

	public void setTraces(TraceList traceList, int[] traceIndices) {
		
		this.combinedTraces = null;
		this.traceIndices = traceIndices;
		showIndividualCheckBox.setEnabled(false);
		
		traceChart.removeAllPlots();
		
		if (traceList == null || traceIndices == null || traceIndices.length == 0) {
            xAxisTitle = "";
			chartPanel.setXAxisTitle("");
			chartPanel.setYAxisTitle("");
			messageLabel.setText("No traces selected");
			add(messageLabel, BorderLayout.NORTH);
			return; 
		}

		remove(messageLabel);
		
		for (int i = 0; i < traceIndices.length; i++) {
			double values[] = new double[traceList.getStateCount()];
			traceList.getValues(traceIndices[i], values);
			DensityPlot plot = new DensityPlot(values, 50);
			plot.setName(traceList.getTraceName(traceIndices[i]));
			plot.setLineStyle(new BasicStroke(1.0f), paints[i % paints.length]);

			traceChart.addPlot(plot);
		}
		
		if (traceIndices.length == 1) {
            xAxisTitle = traceList.getTraceName(traceIndices[0]);
		} else {
            xAxisTitle = traceList.getName();
		}
        chartPanel.setXAxisTitle(xAxisTitle);
		
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

        Plot plot = traceChart.getPlot(0);
        Variate xData = plot.getXData();

        buffer.append(xAxisTitle);
        for (int i = 0; i < traceChart.getPlotCount(); i++) {
            plot = traceChart.getPlot(i);
            buffer.append("\t");
            buffer.append(plot.getName());
        }
        buffer.append("\n");

        for (int i = 0; i < xData.getCount(); i++) {
            buffer.append(String.valueOf(xData.get(i)));
            for (int j = 0; j < traceChart.getPlotCount(); j++) {
                plot = traceChart.getPlot(j);
                Variate yData = plot.getYData();
                buffer.append("\t");
                buffer.append(String.valueOf(yData.get(i)));
            }
            buffer.append("\n");
	}

        return buffer.toString();
    }
}

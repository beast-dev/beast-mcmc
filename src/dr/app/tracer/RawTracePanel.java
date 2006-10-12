/*
 * RawTracePanel.java
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
import dr.gui.chart.Plot;
import org.virion.jam.framework.Exportable;

import javax.swing.*;
import java.awt.*;

/**
 * A panel that displays information about traces
 *
 * @author Andrew Rambaut
 * @author Alexei Drummond
 * @version $Id: RawTracePanel.java,v 1.20 2005/10/03 16:27:36 rambaut Exp $
 */
public class RawTracePanel extends JPanel implements Exportable {

	/**
	 * 
	 */
	private static final long serialVersionUID = 7158704446243557171L;
	private JTraceChart traceChart = new JTraceChart(new LinearAxis(), new LinearAxis());
	private JChartPanel chartPanel = new JChartPanel(traceChart, null, "", "");
	private JLabel messageLabel = new JLabel("No data loaded");
	
	private JCheckBox sampleCheckBox = new JCheckBox("Sample only");
	private JCheckBox linePlotCheckBox = new JCheckBox("Draw line plot");
	private JComboBox legendCombo = new JComboBox(
		new String[] { "None", "Top-Left", "Top", "Top-Right", "Left", 
						"Right", "Bottom-Left", "Bottom", "Bottom-Right" }
	);
	
	/** Creates new RawTracePanel */
	public RawTracePanel() {

		setOpaque(false);

		setMinimumSize(new Dimension(300,150));
		setLayout(new BorderLayout());
		
		JToolBar toolBar = new JToolBar();
		toolBar.setOpaque(false);
		toolBar.setLayout(new FlowLayout(FlowLayout.LEFT));
		toolBar.setFloatable(false);
		
		sampleCheckBox.setSelected(true);
		sampleCheckBox.setOpaque(false);
		toolBar.add(sampleCheckBox);
		
		toolBar.add(new JToolBar.Separator(new Dimension(8,8)));
		linePlotCheckBox.setSelected(true);
		linePlotCheckBox.setOpaque(false);
		toolBar.add(linePlotCheckBox);
		
		toolBar.add(new JToolBar.Separator(new Dimension(8,8)));
		JLabel label = new JLabel("Legend:");
		label.setFont(linePlotCheckBox.getFont());
		label.setLabelFor(legendCombo);
		toolBar.add(label);
		legendCombo.setFont(linePlotCheckBox.getFont());
		legendCombo.setOpaque(false);
		toolBar.add(legendCombo);
		add(messageLabel, BorderLayout.NORTH);
		add(toolBar, BorderLayout.SOUTH);
		add(chartPanel, BorderLayout.CENTER);

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
		
	}
		
	public void setCombinedTraces(CombinedTraces combinedTraces, int[] traceIndices) {

		traceChart.removeAllTraces();
		
		if (combinedTraces == null || traceIndices == null || traceIndices.length == 0) {
			chartPanel.setXAxisTitle("");
			chartPanel.setYAxisTitle("");
			messageLabel.setText("No traces selected");
			return; 
		}

		if (traceIndices.length > 1) {
			chartPanel.setXAxisTitle("");
			chartPanel.setYAxisTitle("");
			messageLabel.setText("Can't show traces of multiple combined traces");
			return; 
		}

		messageLabel.setText("");

		for (int i = 0; i < combinedTraces.getTraceListCount(); i++) {
			TraceList tl = combinedTraces.getTraceList(i);
			
			int n = tl.getStateCount();
			double values[] = new double[n];
			
			int stateStart = tl.getBurnIn();
			int stateStep = tl.getStepSize();

			tl.getValues(traceIndices[0], values);
			traceChart.addTrace(tl.getName(), stateStart, stateStep, values);
		}
		
		chartPanel.setXAxisTitle("State");
		chartPanel.setYAxisTitle(combinedTraces.getTraceName(traceIndices[0]));

		validate();
		repaint();
	}

	public void setTraces(TraceList traceList, int[] traceIndices) {
		
		remove(chartPanel);

		if (traceList == null || traceIndices == null || traceIndices.length == 0) {
			chartPanel.setXAxisTitle("");
			chartPanel.setYAxisTitle("");
			messageLabel.setText("No traces selected");
			return; 
		}

		messageLabel.setText("");

		int n = traceList.getStateCount();
		
		int stateStart = traceList.getBurnIn();
		int stateStep = traceList.getStepSize();

		traceChart.removeAllTraces();
		
		for (int i = 0; i < traceIndices.length; i++) {
			double values[] = new double[n];
			traceList.getValues(traceIndices[i], values);
			traceChart.addTrace(traceList.getTraceName(traceIndices[i]), stateStart, stateStep, values);
		}
		
		chartPanel.setXAxisTitle("State");
		if (traceIndices.length == 1) {
			chartPanel.setYAxisTitle(traceList.getTraceName(traceIndices[0]));
		} else {
			chartPanel.setYAxisTitle(traceList.getName());
		}
		add(chartPanel, BorderLayout.CENTER);
		
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

        double[][] traceStates = new double[traceChart.getPlotCount()][];
        double[][] traceValues = new double[traceChart.getPlotCount()][];
        int maxLength = 0;

        for (int i = 0; i < traceChart.getPlotCount(); i++) {
            plot = traceChart.getPlot(i);
            if (i > 0) {
                buffer.append("\t");
            }
            buffer.append("state");
            buffer.append("\t");
            buffer.append(plot.getName());

            traceStates[i] = traceChart.getTraceStates(i);
            traceValues[i] = traceChart.getTraceValues(i);
            if (traceStates[i].length > maxLength) {
                maxLength = traceStates[i].length;
            }
        }
        buffer.append("\n");

        for (int i = 0; i < maxLength; i++) {
            if (traceStates[0].length > i) {
                buffer.append(Integer.toString((int)traceStates[0][i]));
                buffer.append("\t");
                buffer.append(String.valueOf(traceValues[0][i]));
            } else {
                buffer.append("\t");
            }
            for (int j = 1; j < traceStates.length; j++) {
                if (traceStates[j].length > i) {
                    buffer.append("\t");
                    buffer.append(Integer.toString((int)traceStates[j][i]));
                    buffer.append("\t");
                    buffer.append(String.valueOf(traceValues[j][i]));
                } else {
                    buffer.append("\t\t");
                }
            }
            buffer.append("\n");
	}
	
        return buffer.toString();
	}

}

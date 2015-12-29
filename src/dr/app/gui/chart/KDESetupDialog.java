/*
 * KDESetupDialog.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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

package dr.app.gui.chart;

import dr.app.gui.components.RealNumberField;
import jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.border.EmptyBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

/**
 * @author Marc A. Suchard
 */
public class KDESetupDialog {


	private final JFrame frame;

	private final JCheckBox logXAxis;
	private final JCheckBox logYAxis;

	private final JCheckBox manualXAxis;
	private final JCheckBox manualYAxis;
	private final RealNumberField minXValue;
	private final RealNumberField maxXValue;
	private final RealNumberField minYValue;
	private final RealNumberField maxYValue;

	private final boolean canLogXAxis;
	private final boolean canLogYAxis;
	private final int defaultMinXAxisFlag;
	private final int defaultMaxXAxisFlag;
	private final int defaultMinYAxisFlag;
	private final int defaultMaxYAxisFlag;

	private OptionsPanel optionPanel;

	public KDESetupDialog(final JFrame frame, boolean canLogXAxis, boolean canLogYAxis,
	                        int defaultMinXAxisFlag, int defaultMaxXAxisFlag,
	                        int defaultMinYAxisFlag, int defaultMaxYAxisFlag) {
		this.frame = frame;

		this.canLogXAxis = canLogXAxis;
		this.canLogYAxis = canLogYAxis;

		this.defaultMinXAxisFlag = defaultMinXAxisFlag;
		this.defaultMaxXAxisFlag = defaultMaxXAxisFlag;
		this.defaultMinYAxisFlag = defaultMinYAxisFlag;
		this.defaultMaxYAxisFlag = defaultMaxYAxisFlag;

		logXAxis = new JCheckBox("Log axis");
		manualXAxis = new JCheckBox("Manual range");
		minXValue = new RealNumberField();
		minXValue.setColumns(12);
		maxXValue = new RealNumberField();
		maxXValue.setColumns(12);

		logYAxis = new JCheckBox("Log axis");
		manualYAxis = new JCheckBox("Manual range");
		minYValue = new RealNumberField();
		minYValue.setColumns(12);
		maxYValue = new RealNumberField();
		maxYValue.setColumns(12);

		optionPanel = new OptionsPanel(12, 12);

		optionPanel.addSpanningComponent(new JLabel("X Axis"));
		if (canLogXAxis) {
			optionPanel.addComponent(logXAxis);
		}
		optionPanel.addComponent(manualXAxis);
		final JLabel minXLabel = new JLabel("Minimum Value:");
		optionPanel.addComponents(minXLabel, minXValue);
		final JLabel maxXLabel = new JLabel("Maximum Value:");
		optionPanel.addComponents(maxXLabel, maxXValue);
		manualXAxis.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent changeEvent) {
                final boolean xaSelected = manualXAxis.isSelected();
                minXLabel.setEnabled(xaSelected);
				minXValue.setEnabled(xaSelected);
				maxXLabel.setEnabled(xaSelected);
				maxXValue.setEnabled(xaSelected);
			}
		});
		manualXAxis.setSelected(true);
		manualXAxis.setSelected(false);

		optionPanel.addSeparator();

		optionPanel.addSpanningComponent(new JLabel("Y Axis"));
		if (canLogYAxis) {
			optionPanel.addComponent(logYAxis);
		}
		optionPanel.addComponent(manualYAxis);
		final JLabel minYLabel = new JLabel("Minimum Value:");
		optionPanel.addComponents(minYLabel, minYValue);
		final JLabel maxYLabel = new JLabel("Maximum Value:");
		optionPanel.addComponents(maxYLabel, maxYValue);
		manualYAxis.addChangeListener(new ChangeListener() {
			public void stateChanged(ChangeEvent changeEvent) {
                final boolean yaSelected = manualYAxis.isSelected();
                minYLabel.setEnabled(yaSelected);
				minYValue.setEnabled(yaSelected);
				maxYLabel.setEnabled(yaSelected);
				maxYValue.setEnabled(yaSelected);
			}
		});
		manualYAxis.setSelected(true);
		manualYAxis.setSelected(false);
	}

	public int showDialog(JChart chart) {

		JOptionPane optionPane = new JOptionPane(optionPanel,
				JOptionPane.QUESTION_MESSAGE,
				JOptionPane.OK_CANCEL_OPTION,
				null,
				null,
				null);
		optionPane.setBorder(new EmptyBorder(12, 12, 12, 12));

        Axis xAxis = chart.getXAxis();
        Axis yAxis = chart.getYAxis();

        if (canLogXAxis) {
			logXAxis.setSelected(xAxis instanceof LogAxis);
		}

        if (canLogYAxis) {
			logYAxis.setSelected(yAxis instanceof LogAxis);
		}

		if (!manualXAxis.isSelected()) {
			minXValue.setValue(xAxis.getMinAxis());
			maxXValue.setValue(xAxis.getMaxAxis());
		}

		if (!manualYAxis.isSelected()) {
			minYValue.setValue(yAxis.getMinAxis());
			maxYValue.setValue(yAxis.getMaxAxis());
		}

		final JDialog dialog = optionPane.createDialog(frame, "Kernel Density Estimator Chart");
		dialog.pack();

		dialog.setVisible(true);


		final Integer value = (Integer)optionPane.getValue();
        final int result = (value != null && value != -1) ? value : JOptionPane.CANCEL_OPTION;

        if (result == JOptionPane.OK_OPTION) {
//			if (canLogXAxis) {
//				if (logXAxis.isSelected()) {
//					xAxis = new LogAxis();
//				} else {
//					xAxis = new LinearAxis();
//				}
//                chart.setXAxis(xAxis);
//			}
//
//            if (manualXAxis.isSelected()) {
//				xAxis.setManualRange(minXValue.getValue(), maxXValue.getValue());
//                xAxis.setAxisFlags(Axis.AT_VALUE, Axis.AT_VALUE);
//			} else {
//				xAxis.setAxisFlags(defaultMinXAxisFlag, defaultMaxXAxisFlag);
//			}
//
//			if (canLogYAxis) {
//				if (logYAxis.isSelected()) {
//					yAxis = new LogAxis();
//				} else {
//					yAxis = new LinearAxis();
//				}
//                chart.setYAxis(yAxis);
//            }
//
//            if (manualYAxis.isSelected()) {
//				yAxis.setManualRange(minYValue.getValue(), maxYValue.getValue());
//                yAxis.setAxisFlags(Axis.AT_VALUE, Axis.AT_VALUE);
//			} else {
//				yAxis.setAxisFlags(defaultMinYAxisFlag, defaultMaxYAxisFlag);
//			}
		}

		return result;
	}

}


package dr.gui.chart;

import org.virion.jam.components.RealNumberField;
import org.virion.jam.panels.OptionsPanel;

import javax.swing.*;
import javax.swing.event.ChangeListener;
import javax.swing.event.ChangeEvent;
import javax.swing.border.EmptyBorder;

public class ChartSetupDialog {

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

	public ChartSetupDialog(final JFrame frame, boolean canLogXAxis, boolean canLogYAxis,
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
				minXLabel.setEnabled(manualXAxis.isSelected());
				minXValue.setEnabled(manualXAxis.isSelected());
				maxXLabel.setEnabled(manualXAxis.isSelected());
				maxXValue.setEnabled(manualXAxis.isSelected());
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
				minYLabel.setEnabled(manualYAxis.isSelected());
				minYValue.setEnabled(manualYAxis.isSelected());
				maxYLabel.setEnabled(manualYAxis.isSelected());
				maxYValue.setEnabled(manualYAxis.isSelected());
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

		if (canLogXAxis) {
			logXAxis.setSelected(chart.getXAxis() instanceof LogAxis);
		}
		if (canLogYAxis) {
			logYAxis.setSelected(chart.getYAxis() instanceof LogAxis);
		}

		if (!manualXAxis.isSelected()) {
			minXValue.setValue(chart.getXAxis().getMinAxis());
			maxXValue.setValue(chart.getXAxis().getMaxAxis());
		}

		if (!manualYAxis.isSelected()) {
			minYValue.setValue(chart.getYAxis().getMinAxis());
			maxYValue.setValue(chart.getYAxis().getMaxAxis());
		}

		final JDialog dialog = optionPane.createDialog(frame, "Setup Chart");
		dialog.pack();

		dialog.setVisible(true);

		int result = JOptionPane.CANCEL_OPTION;
		Integer value = (Integer)optionPane.getValue();
		if (value != null && value.intValue() != -1) {
			result = value.intValue();
		}

		if (result == JOptionPane.OK_OPTION) {
			if (canLogXAxis) {
				if (logXAxis.isSelected()) {
					chart.setXAxis(new LogAxis());
				} else {
					chart.setXAxis(new LinearAxis());
				}
			}
			if (manualXAxis.isSelected()) {
				chart.getXAxis().setManualRange(minXValue.getValue(), maxXValue.getValue());
                chart.getXAxis().setAxisFlags(Axis.AT_VALUE, Axis.AT_VALUE);
			} else {
				chart.getXAxis().setAxisFlags(defaultMinXAxisFlag, defaultMaxXAxisFlag);
			}

			if (canLogYAxis) {
				if (logYAxis.isSelected()) {
					chart.setYAxis(new LogAxis());
				} else {
					chart.setYAxis(new LinearAxis());
				}
			}
			if (manualYAxis.isSelected()) {
				chart.getYAxis().setManualRange(minYValue.getValue(), maxYValue.getValue());
                chart.getYAxis().setAxisFlags(Axis.AT_VALUE, Axis.AT_VALUE);
			} else {
				chart.getYAxis().setAxisFlags(defaultMinYAxisFlag, defaultMaxYAxisFlag);
			}
		}

		return result;
	}

}

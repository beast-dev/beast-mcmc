package dr.app.bss;

import java.awt.BorderLayout;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;

import jam.framework.Exportable;
import jam.panels.OptionsPanel;

import javax.swing.JComboBox;
import javax.swing.JComponent;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.SwingConstants;

import org.virion.jam.components.RealNumberField;

@SuppressWarnings("serial")
public class ClockRateModelPanel extends JPanel implements Exportable {

	private BeagleSequenceSimulatorFrame frame;
	private BeagleSequenceSimulatorData data;
	private OptionsPanel optionPanel;

	private JComboBox clockCombo;
	private RealNumberField[] clockParameterFields = new RealNumberField[BeagleSequenceSimulatorData.clockParameterNames.length];

	public ClockRateModelPanel(final BeagleSequenceSimulatorFrame frame,
			final BeagleSequenceSimulatorData data) {

		this.frame = frame;
		this.data = data;

		setOpaque(false);
		setLayout(new BorderLayout());

		optionPanel = new OptionsPanel(12, 12, SwingConstants.CENTER);
		add(optionPanel, BorderLayout.NORTH);

		clockCombo = new JComboBox();
		clockCombo.setOpaque(false);

		for (String clockModel : BeagleSequenceSimulatorData.clockModels) {
			clockCombo.addItem(clockModel);
		}// END: fill loop

		clockCombo.addItemListener(new ListenClockCombo());

		for (int i = 0; i < BeagleSequenceSimulatorData.clockParameterNames.length; i++) {
			clockParameterFields[i] = new RealNumberField();
			clockParameterFields[i].setColumns(8);
			clockParameterFields[i].setValue(data.clockParameterValues[i]);
		}// END: fill loop

		setClockArguments();

	}// END: Constructor

	private void setClockArguments() {

		optionPanel.removeAll();
		optionPanel.addComponents(new JLabel("Clock rate model:"), clockCombo);
		optionPanel.addSeparator();
		optionPanel.addLabel("Set parameter values:");

		int index = clockCombo.getSelectedIndex();

		for (int i = 0; i < data.clockParameterIndices[index].length; i++) {

			int k = data.clockParameterIndices[index][i];

			JPanel panel = new JPanel(new BorderLayout(6, 6));
			panel.add(clockParameterFields[k], BorderLayout.WEST);
			panel.setOpaque(false);
			optionPanel.addComponentWithLabel(
					BeagleSequenceSimulatorData.clockParameterNames[k] + ":",
					panel);

		}// END: indices loop

		validate();
		repaint();
	}// END: setClockArguments

	private class ListenClockCombo implements ItemListener {
		public void itemStateChanged(ItemEvent ie) {

			setClockArguments();
			frame.fireModelChanged();

		}// END: actionPerformed
	}// END: ListenClockCombo

	public void collectSettings() {

		data.clockModel = clockCombo.getSelectedIndex();
		for (int i = 0; i < BeagleSequenceSimulatorData.clockParameterNames.length; i++) {

			data.clockParameterValues[i] = clockParameterFields[i].getValue();

		}// END: fill loop
	}// END: collectSettings
	
	@Override
	public JComponent getExportableComponent() {
		return this;
	}// END: getExportableComponent

}// END: class
